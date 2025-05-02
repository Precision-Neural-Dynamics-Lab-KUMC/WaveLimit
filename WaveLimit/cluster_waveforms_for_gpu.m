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
        
        num_units = size(test_mean_waveforms,2);
        test_mean_waveforms = [test_mean_waveforms,zeros(size(mean_waveforms,1),1)]; %Add zero waveform
        % Calculate Euclidean distances between waveforms and cluster centers
        % distances = sqrt(sum((repmat(permute(waveforms, [3 2 1]), [size(test_mean_waveforms,2),1,1]) - repmat(permute(test_mean_waveforms, [2 3 1]), [1,size(waveforms,2),1])).^2,3));
        noise_waveforms = repmat(permute(waveforms, [3 2 1]), [size(test_mean_waveforms,2),1,1]) - repmat(permute(test_mean_waveforms, [2 3 1]), [1,size(waveforms,2),1]);
        distances = sum(noise_waveforms.^2,3);
        test_mean_waveforms = [test_mean_waveforms(:,1:(end-1))]; %Remove zero waveform

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
        sorted_distances =  gpuArray.zeros(num_units,size(pen_distances,2));
        sort_order =  gpuArray.zeros(num_units,size(pen_distances,2));
        cluster_order = gpuArray.zeros(num_units,size(pen_distances,2));
        for n = 1:size(test_mean_waveforms,2)
            [sorted_distances(n,:), sort_order(n,:)] = sort(pen_distances(n,:),2);  %Order waveforms by their distance from the test spike template
            cluster_order(n,sort_order(n,:)) = gpuArray(1:num_waveforms);     
        end     
        
      

        % % Calculate the difference in sorted distances to estimate how many waveforms are in a  cluster using chi^2 statistics
        % diff_sorted_distances = diff(sorted_distances,1,2); 
        % [b,a] = butter(1,.001,'low');
        % %Low-pass to improve estimating the minimum distance rate change
        % smooth_diff_sorted_distances = filter(b,a,diff_sorted_distances,[],2); % Low-pass filtering
        % smooth_diff_sorted_distances(:,1:101) = gpuArray(pen_value);  % Ignore the first 101 samples
        
        % % Find the mode of the chi^2 distribution for each cluster by
        % % finding the minimum change in distance from one waveform distance
        % % to the next
        % min_point = gpuArray.zeros(size(num_possible_waveforms));
        % for n = 1:size(smooth_diff_sorted_distances,1)
        %     if num_possible_waveforms(n)>0
        %         [~, min_point(n)] = min(smooth_diff_sorted_distances(n,1:num_possible_waveforms(n)));
        %     else
        %         min_point(n) = 0;
        %     end
        % end
        % min_point = min([num_possible_waveforms./2; min_point]);  %Make sure that it at most includes only half of all possible waveforms (real waveforms + outliers) for this cluster
        % 
        %Set up an inverse chi^2 function that can be scaled to map a
        %cumulative probility to an expected distance value which
        %should be distributed like the distance values
        % chi2fun = @(x,xdata) chi2inv(x(2).*xdata,num_time_pts)./x(1);  %chi2inv(P,v) p = probability, v = degrees of freedom, in our case the number of time points
        chi2fun = @(x,xdata) chi2inv(x(2).*xdata,x(1))./x(3) + x(4);  %chi2inv(P,v) p = probability, v = degrees of freedom, related to number of time points

        options = optimset('Display','off');
        
        for n = 1:num_units  %For each unit
            curr_limit = sum(sorted_distances(n,1:num_possible_waveforms(n))< 2*median(sorted_distances(n,1:num_possible_waveforms(n))));
            cov_mat = cov(squeeze(noise_waveforms(n,sort_order(n,1:curr_limit),:)));
            [approx_dof,est_scaleFact,est_offset]  = estimateGenChi2(cov_mat);
            est_25th_thresh(n) = est_scaleFact*chi2inv(0.25, approx_dof) + est_offset;  %Use the estimated chi^2 dist to determine how many waveforms make up 25th percentile of waveforms
            idx = find(sorted_distances(n,1:num_possible_waveforms(n)) < est_25th_thresh(n), 1, 'last');
            if num_possible_waveforms(n) < sorting_options.min_num_waveforms || isempty(idx)  % Safety check to avoid errors when the unit has small number of waveforms or no waveforms less than threshold 
                est_25th(n) = round(0.25*num_possible_waveforms(n));
            else
                est_25th(n) = idx;
            end
            est_25th(n) =  min([est_25th(n), round(0.25*num_possible_waveforms(n))]);
            
            % x values go from 0 to 0.25 (percentile) on the inverted cdf
            x_values = linspace(0.25./est_25th(n), 0.25, est_25th(n));
            % y values are the distance values that were predicted to
            % be the 1st to 25th percentile for the given cluster
            y_values = sorted_distances(n,1:est_25th(n));
            %Do least squared curve fitting to scale the distance
            %values expected for each cumulative probability
            if round(length(x_values)/100)>1
                try
                    [x,~] = lsqcurvefit(chi2fun, [gather(approx_dof), 1, gather(chi2inv(.25, approx_dof)./sorted_distances(n,est_25th(n))), gather(est_offset)], gather(x_values(1:round(length(x_values)/100):end)), gather(y_values(1:round(length(x_values)/100):end)), [1, 0.5, eps, 0], [num_time_pts, 3, Inf, Inf], options); %, options);
                    chi2_dof(n) = x(1);
                    prob_sf(n) = x(2);
                    chi2_sf(n) = x(3);
                    chi2_offset(n) = x(4);
                catch
                    disp('lsqfit failure')
                    chi2_dof(n) = approx_dof;
                    prob_sf(n) = chi2inv(.25, approx_dof)./gather(sorted_distances(n,est_25th(n)));  %TODO check this
                    chi2_sf(n) = est_scaleFact;
                    chi2_offset(n) = est_offset;
                end
            else
                chi2_dof(n) = approx_dof;
                prob_sf(n) = chi2inv(.25, approx_dof)./gather(sorted_distances(n,est_25th(n)));
                chi2_sf(n) = est_scaleFact;
                chi2_offset(n) = est_offset;
            end
        end

        %Test plot of fit
        % n = 1;
        % tmp = chi2fun([chi2_dof(n),prob_sf(n),chi2_sf(n), chi2_offset(n)], linspace(0,1,est_25th(n)*4));
        % figure;
        % plot(sorted_distances(n,:))
        % hold on
        % plot(tmp)

        % Scale distances and compute cumulative probability using chi^2
        chi2_sf = sorting_options.scale_factor*chi2_sf;
        for n = 1:num_units
            chi2_cdf(n,:) = chi2cdf(sorted_distances(n,:).*chi2_sf(n)-chi2_offset(n),chi2_dof(n));
        end
        
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
        
        % Estimate fraction of true spikes based on chi^2 fit
        fract_true_chi2 = gpuArray(0.5*ones(num_units,length(timestamps)));  %Preallocate
        spike_prior = est_25th./sum(est_25th);  %Estimate prior probability of spike, Note this uses 25th percentile. Both numerator and denominator would be multiplied by 4 to estimate total number of spikes (but they then cancel)
        for n = 1:num_units
            [~, sort_index] = sort(distances(n,:));
            sorted_dist_for_cluster = distances(:,sort_index);
            sorted_dist_for_cluster(sorted_dist_for_cluster==0) = eps;  %Can't have distance be equal to zero for chi2pdf
            for curr_n = 1:num_units
                chi2_pdf(curr_n,:) = chi2pdf(sorted_dist_for_cluster(curr_n,:).*chi2_sf(curr_n)-chi2_offset(curr_n),chi2_dof(curr_n));
            end
            tmp = sorted_dist_for_cluster(n,:).*chi2_sf(n)-chi2_offset(n);
            chi2_pdf(n,tmp<(chi2_dof(n)-2)) = chi2pdf(chi2_dof(n)-2,chi2_dof(n));
            chi2_pdf = chi2_pdf.*repmat(spike_prior',[1,size(chi2_pdf,2)]);
            cumsum_chi2_pdf = cumsum(chi2_pdf,2);
            fract_true_chi2(n,:) = cumsum_chi2_pdf(n,:)./sum(cumsum_chi2_pdf,1);
        end
        
        % fract_true_chi2 = gpuArray(0.5*ones(num_units,length(timestamps)));
        % chi2_pdf = diff(chi2_cdf,1,2);
        % for n = 1:num_units
        %      [~, sort_index] = sort(distances(n,:));
        %     sorted_dist_for_cluster = distances(:,sort_index);
        %     sorted_dist_for_cluster(sorted_dist_for_cluster==0) = eps;  %Can't have distance be equal to zero for chi2pdf
        %     for curr_n = 1:num_units
        %         chi2_prob(n,:) = chi2cdf(sorted_distances(n,:).*chi2_sf(n)-chi2_offset(n),chi2_dof(n));
        %     end
        % end

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

% [pc_components,wave_pcs,~,~,~,wave_mu] = pca(waveforms');
% decimation_factor = 10;
% example_a = sort_order(1,1:(new_cutoffs(1)))'; 
% example_b = sort_order(2,1:(new_cutoffs(2)))'; 
% example_c = sort_order(3,1:(new_cutoffs(3)))'; 
% example_u = setdiff(1:num_waveforms, [example_a; example_b; example_c])';
% example_a = example_a(1:decimation_factor:end);
% example_b = example_b(1:decimation_factor:end);
% example_c = example_c(1:decimation_factor:end);
% example_u = example_u(1:decimation_factor:end);
% example_all = [example_a; example_b; example_c; example_u];
% throwout_i = any(waveforms(:,example_all)>0.15,1);
% example_all = example_all(~throwout_i);
% throwout_i = any(waveforms(:,example_u)>0.15,1);
% example_u = example_u(~throwout_i);
% %%
% t = 1e6*((1:size(waveforms,1))-9)/30000;
% mean_pcs = bsxfun(@minus,test_mean_waveforms',wave_mu)*pc_components;
% % tmp = bsxfun(@minus,waveforms',wave_mu)*pc_components;
% 
% 
% colors = [ [0.45, 0.65,0.25]; 
%           [0.3333, 0.4196,0.1843]; % Main green
%           [0.5333, 0.5937, 0.4290]; % Lighter green
%           [0.5, 0.5, 0.5]; % grey
%          ];
% 
% fig_2d_right = figure
% hold on
% h4 = scatter(wave_pcs(example_u, 1), wave_pcs(example_u, 2), 10, [0.2 0.2 0.2], 'filled');
% h1 = scatter(wave_pcs(example_a, 1), wave_pcs(example_a, 2), 10, colors(1,:), 'filled');
% h2 = scatter(wave_pcs(example_b, 1), wave_pcs(example_b, 2), 10, colors(2,:), 'filled');
% h3 = scatter(wave_pcs(example_c, 1), wave_pcs(example_c, 2), 10, colors(3,:), 'filled');
% view(0, 90)
% grid off
% 
% % Create dummy scatter plots for the legend
% hold on
% % Create dummy plots for the legend
% h1_legend = plot(nan, nan, 'Marker', '.', 'MarkerSize', 25, 'DisplayName', 'Unit A', 'LineStyle', 'none', 'Color', h1.CData(1,:));
% h2_legend = plot(nan, nan, 'Marker', '.', 'MarkerSize', 25, 'DisplayName', 'Unit B', 'LineStyle', 'none', 'Color', h2.CData(1,:));
% h3_legend = plot(nan, nan, 'Marker', '.', 'MarkerSize', 25, 'DisplayName', 'Unit C', 'LineStyle', 'none', 'Color', h3.CData(1,:));
% h4_legend = plot(nan, nan, 'Marker', '.', 'MarkerSize', 25, 'DisplayName', 'Unsorted', 'LineStyle', 'none', 'Color', h4.CData(1,:));
% 
% % Create legend and set properties
% lgd = legend([h1_legend, h2_legend, h3_legend, h4_legend]);
% set(lgd, 'Box', 'off'); % Remove the box around the legend
% lgd.FontSize = 15; % Set the font size of the legend
% xlabel('PC 1')
% ylabel('PC 2')
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])
% %print(gcf, 'PostSort_PC','-dpng')
% set(gcf, 'color', 'w'); % Set background color to white
% set(gca, 'FontSize', 15); % Set font size to 15
% set(gca, 'LineWidth', 3); % Make the plot lines box thicker
% 
% xlim([ -0.5231, 0.4769])
% ylim([-0.3312, 0.4688])

end




function SNR = signal2noise_PtoP2(waves)
%Helper function, Peak to Peak signal-to-noise divided by 2
mean_signal = mean(waves,2);
signal = max(mean_signal) - min(mean_signal);
remain_noise = waves-repmat(mean_signal, [1, size(waves,2)]);
noise = std(remain_noise(:));
SNR = (signal./noise)/2;
end

function [approx_dof,scaleFact,offset]  = estimateGenChi2(cov_mat)
%This function approximates a generalized chi-squared distribution by using
%the covariance matrix of the noise distribution around 
%It then fits a chi-squared distribution by approximating the degrees of
%freedom, and scaling and adding an offset using a "method of moments"
%approach

%Diagonalize covariance matrix to obtain equivalent sum of independent
%normal random variables
[~,D] = eig(cov_mat);
d_inc = diag(D)>(max(diag(D))/100); %Remove 

var_values = diag(D(d_inc,d_inc));
expected_val = sum(var_values);  %Expected value of sum of squared random variables
expected_var = 2*sum((var_values).^2);  %Expected variance of sum of squared random variables
approx_dof = sum(var_values)./max(var_values); %Approximated DOF of chi-square distribution to use,  
% the mean of the chi^2 distribution = approx_dof, variance = 2*approx_dof

scaleFact = sqrt( expected_var/(2*approx_dof) );
offset = expected_val - scaleFact*approx_dof;
end















