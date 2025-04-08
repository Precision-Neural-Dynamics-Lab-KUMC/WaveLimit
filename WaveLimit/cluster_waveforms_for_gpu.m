function [test_mean_waveforms,cutoff_distances] = cluster_waveforms_for_gpu(waveforms, timestamps, mean_waveforms, sorting_options, T, total_spikes)
% Clusters waveforms using GPU acceleration.
% Inputs:
% - waveforms: Matrix of waveforms to cluster.
% - timestamps: Timestamps corresponding to waveforms.
% - mean_waveforms: Initial mean waveforms for clustering.
% - sorting_options: Options for clustering process.
% - T: Total recording time.
% - total_spikes: Total number of spikes (waveforms may be a subset for faster clustering).
% Adam Rouse, 4/24/20, v1.3

% Default options if not provided
if nargin<4
    sorting_options = default_options();
end
if nargin<5 || isempty(T)
    T = max(timestamps);  % Use maximum timestamp if total time not provided
end
if nargin<6 || isempty(total_spikes)
    total_spikes = size(waveforms,1);  %If total number spikes specifice, use number of waveforms;  There can be more total spikes if waveforms are only a subset (for speed of clustering)
end
% Transfer data to GPU for faster computation
timestamps = gpuArray(timestamps);
waveforms = gpuArray(waveforms);
[timestamps, sort_i] = sort(timestamps);  %Make sure time stamps are sorted
waveforms = waveforms(:,sort_i);


% Ensure we proceed with clustering if mean waveforms are provided and enough waveforms exist
[num_time_pts,num_waveforms] = size(waveforms);
if exist('mean_waveforms', 'var')  && num_waveforms > sorting_options.min_num_waveforms
    
    % Initialization before clustering
    keep_iterating = true;
    update_mean_flag = ones(size(mean_waveforms,2),1);  % Track which clusters need updating
    test_mean_waveforms = mean_waveforms;              % Initialize test mean waveforms
    old_best_true = gpuArray(zeros(size(mean_waveforms,2),1));   % Best fraction of true waveforms from previous iterations, zero to start so first sorting will now be the best
    num_iter = 0;                                       %Iteration counter
    multiunit_counter = zeros(size(mean_waveforms,2),1); % Track clusters that may represent multiple units
    
    %Iterate through the clustering algorithm , each step the centers are updated, if the estimated
    %fraction of total waveforms improves, keep iterating, once it no
    %longer improves or gets worse, stop and go back to the last
    %cluster centers, the iteration continues until all of the clusters
    %have reached their local maximum
    while keep_iterating

        % Stop iterating when no clusters need further updating
        if (~any(update_mean_flag))
            keep_iterating = false;
        end
        
       
       % Calculate Euclidean distances between waveforms and cluster centers
        distances = sqrt(sum((repmat(permute(waveforms, [3 2 1]), [size(test_mean_waveforms,2),1,1]) - repmat(permute(test_mean_waveforms, [2 3 1]), [1,size(waveforms,2),1])).^2,3));

        % Penalize distances from non-closest clusters so waveform not included in any cluster but the closest
        pen_value = ceil(max(distances(:)));
        %Penalized distances adds pen_value to every distance that is not the
        %closest cluster center, for example if the distances for a spike
        %are .3, .4, .5 from the 3 cluster centers and pen_value = 1, the penalized
        %distances are .3, 1.4, 1.5
        pen_distances = distances;
        [~,min_d] = min(distances,[],1);
        for n = 1:size(distances,1)
            pen_distances(n,min_d~=n) = distances(n,min_d~=n)+pen_value;
        end
        
        % Determine number of waveforms closest to each cluster center
        num_possible_waveforms = histcounts(min_d, [0.5, (1:size(test_mean_waveforms,2))+0.5]);
        num_possible_waveforms(num_possible_waveforms == num_waveforms) = num_waveforms-1;  %To prevent a later error when there's only one cluster
        
        %Calculate distance from a flat signal of zeros, no spikes farther
        %away then this should be included
        max_distances = sqrt(sum((test_mean_waveforms.^2)));
        for n = 1:length(num_possible_waveforms)
            num_possible_waveforms(n) = min(num_possible_waveforms(n), sum(pen_distances(n,:)<max_distances(n)));
        end

        % Set sorting and cluster order
        sorted_distances =  gpuArray.zeros(size(pen_distances));
        sort_order =  gpuArray.zeros(size(pen_distances));
        cluster_order = gpuArray.zeros(size(sorted_distances));
        for n = 1:size(test_mean_waveforms,2)
            [sorted_distances(n,:), sort_order(n,:)] = sort(pen_distances(n,:),2);  %Order waveforms by their distance from the test spike template
            cluster_order(n,sort_order(n,:)) = gpuArray(1:num_waveforms);     
        end     
        
        % Calculate the difference in sorted distances to estimate how many waveforms are in a  cluster using chi^2 statistics
        diff_sorted_distances = diff(sorted_distances,1,2); 
        [b,a] = butter(1,.001,'low');
        %Low-pass to improve estimating the minimum distance rate change
        smooth_diff_sorted_distances = filter(b,a,diff_sorted_distances,[],2); % Low-pass filtering
        smooth_diff_sorted_distances(:,1:101) = gpuArray(pen_value);  % Ignore the first 101 samples
        
        % Find the mode of the chi^2 distribution for each cluster by
        % finding the minimum change in distance from one waveform distance
        % to the next
        min_point = gpuArray.zeros(size(num_possible_waveforms));
        for n = 1:size(smooth_diff_sorted_distances,1)
            if num_possible_waveforms(n)>0
                [~, min_point(n)] = min(smooth_diff_sorted_distances(n,1:num_possible_waveforms(n)));
            else
                min_point(n) = 0;
            end
        end
        min_point = min([num_possible_waveforms./2; min_point]);  %Make sure that it at most includes only half of all possible waveforms (real waveforms + outliers) for this cluster
        
        %Set up an inverse chi^2 function that can be scaled to map a
        %cumulative probility to an expected distance value which
        %should be distributed like the distance values
        chi2fun = @(x,xdata) chi2inv(x(2).*xdata,num_time_pts)./x(1);  %chi2inv(P,v) p = probability, v = degrees of freedom, in our case the number of time points
        options = optimset('Display','off');
        
        chi2_p_at_mode = chi2cdf(num_time_pts-2,num_time_pts); %The mode of the chi^2 distribution is at degrees of freedom - 2 (num_time_pts-2), use that to determine the cumulative probability at the mode
        est_25th = round((0.25/chi2_p_at_mode)*min_point);  %Use the estimated mode value, to determine how many waveforms make up 25th percentile of waveforms
        for n = 1:size(sorted_distances,1)  %For each unit
            % x values go from 0 to 0.25 (percentile) on the inverted cdf
            x_values = linspace(0.25./est_25th(n), 0.25, est_25th(n));
            % y values are the distance values that were predicted to
            % be the 1st to 25th percentile for the given cluster
            y_values = sorted_distances(n,1:est_25th(n));
            %Do least squared curve fitting to scale the distance
            %values expected for each cumulative probability
            if round(length(x_values)/100)>1
                try
                    [x,~] = lsqcurvefit(chi2fun, [chi2inv(.25, num_time_pts)./gather(sorted_distances(n,est_25th(n))), 1], gather(x_values(1:round(length(x_values)/100):end)), gather(y_values(1:round(length(x_values)/100):end)), [0 0], [1e10 1e10], options); %, options);
                    chi2_sf(n) = x(1);
                catch
                    disp('lsqfit failure')
                    chi2_sf(n) = chi2inv(.25, num_time_pts)./gather(sorted_distances(n,est_25th(n)));
                end
            else
                chi2_sf(n) = 1000;
            end
        end
        % Scale distances and compute cumulative probability using chi^2
        chi2_sf = sorting_options.scale_factor*chi2_sf;
        chi2_cdf = chi2cdf(sorted_distances.*repmat(chi2_sf',1,size(sorted_distances,2)),num_time_pts);
        
        %Estimate fraction of true spikes based on ISI violations
        fract_true_ISI = 0.5*gpuArray.ones(size(test_mean_waveforms,2),length(timestamps));  %Preallocate
        ISI_waveforms = gpuArray.false(size(test_mean_waveforms,2),length(timestamps));
        ISI = diff(timestamps);
        ISI_viol = ISI<sorting_options.max_ISI;
        for n = 1:size(test_mean_waveforms,2)
            %Find when each ISI violation occurs in the clusterning order for the current unit
            ISI_cluster_order = [cluster_order(n,[ISI_viol; false]); cluster_order(n,[false; ISI_viol])];
            ISI_cluster_order = max(ISI_cluster_order, [], 1); %Always assumes the spike added later to the cluster is the ISI violation of the short interval pair
            ISI_cluster_order = unique(ISI_cluster_order);
            ISI_waveforms(n,ISI_cluster_order) = 1;
            cumsum_ISI_waveforms = cumsum(ISI_waveforms(n,:));
            tmp = cumsum_ISI_waveforms.*T./(2.*(sorting_options.max_ISI-sorting_options.min_ISI).*(((total_spikes./length(cumsum_ISI_waveforms)).*(1:length(cumsum_ISI_waveforms))).^2)); %Part of fract ISI equation
            fract_true_ISI(n,tmp<0.25) = sqrt(.25 - tmp(tmp<0.25))+.5;
        end
        
        %Estimate fraction of true spikes based on chi^2 fit
        fract_true_chi2 = gpuArray(0.5*ones(size(test_mean_waveforms,2),length(timestamps)));  %Preallocate
        spike_prior = min_point./sum(min_point);  %Estimate prior probability of spike
        for n = 1:size(test_mean_waveforms,2)
            [~, sort_index] = sort(distances(n,:));
            sorted_dist_for_cluster = distances(:,sort_index);
            sorted_dist_for_cluster(sorted_dist_for_cluster==0) = eps;  %Can't have distance be equal to zero for chi2pdf
            chi2_pdf = chi2pdf(sorted_dist_for_cluster.*repmat(chi2_sf',1,size(distances,2)),28);
            chi2_pdf = chi2_pdf.*repmat(spike_prior',[1,size(chi2_pdf,2)]);
            cumsum_chi2_pdf = cumsum(chi2_pdf,2);
            fract_true_chi2(n,:) = cumsum_chi2_pdf(n,:)./sum(cumsum_chi2_pdf,1);
        end
        
        %To be conservative, we're going to estimate the percentage of true waveforms as the minimum of either the percentage based on ISI violations or chi2 overlap (the case of one unit on a channel will only use ISI measure)
        overall_fract_true = min(fract_true_ISI, fract_true_chi2);  
        ratio_true_waveforms = chi2_cdf.*overall_fract_true;  %the phi value in the manuscript that should be maximized
        [new_best_true, cutoffs] = max(ratio_true_waveforms,[],2);
        cutoffs(cutoffs>num_possible_waveforms') = num_possible_waveforms(cutoffs'>num_possible_waveforms)';
        
        %Check to see if each unit has a cutoff such that it includes
        %at least 50% of waveforms based on the chi2 model and that at
        %least 50% of those waveforms are true waveforms based on ISI
        %violations or chi2 overlap
        %If it's not, we can't isolate a single unit and we're going to
        %call it a multi unit
        multiunit_flag = ~any(chi2_cdf> 0.5 & overall_fract_true > 0.5,2);
        multiunit_counter(~multiunit_flag) = 0;
        multiunit_counter = multiunit_counter + multiunit_flag;
        mean_shortISIwaveforms = gpuArray.zeros(size(test_mean_waveforms));
        new_cutoffs = gpuArray.zeros(size(cutoffs));
        SNR = gpuArray.zeros(size(cutoffs));
        for n = 1:size(test_mean_waveforms,2)  %Iterate through each cluster
            if (cutoffs(n)>sorting_options.min_num_in_cluster) %There needs to be at least 500 waveforms in the cluster
                if ~multiunit_flag(n)
                    new_cutoffs(n) = cutoffs(n);
                else
                    %If its a multiunit, we just include 95% of waveforms
                    %using the chi2 model
                    [~, new_cutoffs(n)] = min(abs(chi2_cdf(n,:)- .95));
                    n_spk = new_cutoffs(n);
                    n_tot_spk = (total_spikes./length(timestamps))*n_spk ; %Correct for ISI based on total number of spikes
                    ISI = diff(sort(timestamps(cluster_order(n,1:(n_spk)))));
                    R = sum(ISI<sorting_options.max_ISI);
                    tmp = R*T/(2*(sorting_options.max_ISI-sorting_options.min_ISI)*(n_tot_spk^2));
                    if tmp<0.25
                        overall_fract_true = sqrt(.25 - tmp) + .5;
                    else
                        overall_fract_true = 0.5;
                    end
                    new_best_true(n) = chi2_cdf(n,n_spk).*overall_fract_true;
                    if multiunit_counter(n) >= sorting_options.max_multi_iterations  %Try to move characteristic waveform for only a few iterations if its always been multiunit anyways
                        update_mean_flag(n) = 0;
                    end
                end
                SNR(n) = signal2noise_PtoP2(waveforms(:,cluster_order(n,1:(new_cutoffs(n)))));
                %If we should still try to update the mean spike
                %template
                if (update_mean_flag(n))
                    %If the new template gives a larger amount of true
                    %waveforms, then we should continue to search to find
                    %a better template
                    if (new_best_true(n) > old_best_true(n))
                        %Save this current test mean spike as the mean spike
                        mean_waveforms(:,n) = test_mean_waveforms(:,n);
                        %Determine the mean of the short ISI waveforms
                        shortISIwaveforms = waveforms(:,ISI_waveforms(n,:) & chi2_cdf(n,cluster_order(n,:))<0.99);
                        if (size(shortISIwaveforms,2) == 0)
                            mean_shortISIwaveforms(:,n) = test_mean_waveforms(:,n);
                            update_mean_flag(n) = 0;
                        elseif (size(shortISIwaveforms,2) == 1)
                            mean_shortISIwaveforms(:,n) = shortISIwaveforms;
                        else
                            mean_shortISIwaveforms(:,n) = mean(shortISIwaveforms,2);
                        end
                        
                        %Update the test_mean_wavforms to try to move away from short ISI waveforms
                        if (update_mean_flag(n) == 1)
                            test_mean_waveforms(:,n) = gather(test_mean_waveforms(:,n) + (sorting_options.update_percentage/100)*(test_mean_waveforms(:,n)-mean_shortISIwaveforms(:,n)));
                        end
                        %Save the best score, to compare to the next
                        %iteration
                        old_best_true(n) = new_best_true(n);
                    else
                        %If it didn't improve, then we stop updating and
                        %use the old template mean_waveforms, and save the old
                        %percent true.
                        update_mean_flag(n) = 0;
                        new_best_true(n) = old_best_true(n);
                        test_mean_waveforms(:,n) = mean_waveforms(:,n);
                    end
                end
            else %If this iteration had less than 500 waveforms
                if (old_best_true(n) ~= 0)  %Check to see if a previous iteration was better and use that
                    update_mean_flag(n) = 0;
                    new_best_true(n) = old_best_true(n);
                else  %Otherwise, there's no good cluster for this mean spike
                    new_best_true(n) = 0;
                    SNR(n) = 0;
                    update_mean_flag(n) = 0;
                    new_cutoffs(n) = 0;
                end
            end
        end
        %Display iteration
        num_iter = num_iter+1;
        display(['Number of iterations: ', num2str(num_iter)])
        if (num_iter>sorting_options.max_iterations)
            keep_iterating = false;
            warning(['Max number of iterations (' num2str(sorting_options.max_iterations) ') completed'], 'User:cluster_waveforms')
        end
    end
    %Organize the cutoff_distances to return
    cutoff_distances = zeros(size(new_cutoffs));
    for n = 1:length(new_cutoffs)
        if new_cutoffs(n)>0
            cutoff_distances(n) = gather(sorted_distances(n,new_cutoffs(n)));
        else
            cutoff_distances(n) = 0;
        end
    end
else
    test_mean_waveforms = zeros(size(waveforms,1), 0);
    cutoff_distances = [];
end

end




function SNR = signal2noise_PtoP2(waves)
%Helper function, Peak to Peak signal-to-noise divided by 2
mean_signal = mean(waves,2);
signal = max(mean_signal) - min(mean_signal);
remain_noise = waves-repmat(mean_signal, [1, size(waves,2)]);
noise = std(remain_noise(:));
SNR = (signal./noise)/2;
end

















