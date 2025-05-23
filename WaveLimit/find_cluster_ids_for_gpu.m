function template_waveforms = find_cluster_ids_for_gpu(waveforms, timestamps, options, tot_spikes)
% find_cluster_ids for WaveLimit, v1.3
% Adam Rouse, 4/23/20
% Identifies cluster waveforms on the GPU using PCA and template matching.
% Inputs:
% - waveforms: Matrix of waveforms to cluster.
% - timestamps: Timestamps corresponding to the waveforms.
% - options: Configuration options for clustering.
% - tot_spikes: Total number of spikes for scaling interspike interval (ISI) calculations.


% Move data to GPU and sort waveforms by ascending timestamps
timestamps = gpuArray(timestamps);
[timestamps,sort_i] = sort(timestamps,'ascend');
waveforms =  gpuArray(waveforms(:,sort_i));
num_waveforms = size(waveforms,2);
noise_level = std(waveforms(1,:)); % Estimate noise level using the first time sample

% Perform clustering if there are more than 1000 waveforms and the noise level is positive
if num_waveforms>1000 && noise_level>0
  
    % Apply a scale factor based on noise level (to address magnitude dependency issue)
    scale_factor = 0.02/noise_level;
    waveforms = scale_factor*waveforms;
    noise_level = std(waveforms(1,:));  % Re-estimate noise level after scaling

     % Perform PCA (Principal Component Analysis) using SVD (Singular Value Decomposition)
    [U,sigma] = svd(bsxfun(@minus,waveforms,nanmean(waveforms,2))', 'econ');
    pca_score =  bsxfun(@times,U,diag(sigma)');

    small_sigma = 1.0*noise_level;
    large_sigma = 2.5*noise_level;

    % Set random seed for reproducibility if specified in options
    if options.repeatable_random_seed~=0
        rng(options.repeatable_random_seed);
    end

    %Here, we're going to identify 1000 waveforms to test if they look like
    %they are near cluster centers. The first 500 are random. Then we'll
    %add parts of the waveform space that look like they haven't been
    %sampled well to fill in 501-750 and 751-1000 to get our 1000
    % Select a random subset of 500 waveforms for sampling
    rand_sample = gpuArray(randperm(num_waveforms,500));
    % Calculate pairwise distances between samples in PCA space
    sample_dist =  repmat(sqrt( sum( (permute(pca_score(rand_sample(1:500),1:5),[1,3,2])-permute(pca_score(:,1:5),[3,1,2])).^2 ,3) ), [1,1,2]);
    % Filter out waveforms with fewer than 50 neighbors for clustering since these can't be cluster centers
    sample_dist(sum(sample_dist(:,:,1)<large_sigma,2)<50,:,:) = NaN;  
    sample_dist(:,:,1) = exp(-0.5 * ((  sample_dist(:,:,1)  )./small_sigma).^2) ./ (sqrt(2*pi) .* small_sigma);
    sample_dist(:,:,2) = exp(-0.5 * ((  sample_dist(:,:,2)  )./large_sigma).^2) ./ (sqrt(2*pi) .* large_sigma);
    % Handle NaN values by replacing them with small values to avoid divide-by-zero errors
    sample_dist(isnan(sample_dist(:,1,1)),:,1) = 0;
    sample_dist(isnan(sample_dist(:,1,2)),:,2) = 0.001; %Make small number rather than 0 to avoid divide by zero

    % Identify spikes that are poorly clustered (missed spikes)
    missed_spikes = find(sum(sample_dist(:,:,1),1)<1);
    if length(missed_spikes)>= 250  %If enough spikes outside concentrated areas, then random selection
        rand_sample(501:750) = missed_spikes(randperm(length(missed_spikes),250));
    else
        [~,sort_spikes_i] = sort(sum(sample_dist(:,:,1),1),'ascend');  %Otherwise, sort by distances from sample
        rand_sample(501:750) = sort_spikes_i(1:250);
    end

    % Calculate pairwise distances between samples in PCA space for new sample set
    sample_dist(501:750,:,:) =  repmat(sqrt( sum( (permute(pca_score(rand_sample(501:750),1:5),[1,3,2])-permute(pca_score(:,1:5),[3,1,2])).^2 ,3) ), [1,1,2]);
    % Filter out waveforms with fewer than 50 neighbors for clustering since these can't be cluster centers
    sample_dist([false(500,1);true(250,1)] & sum(sample_dist(:,:,1)<large_sigma,2)<50,:,:) = NaN;  %If the test waveform has less than 50 neighbors, it's not a cluster center
    sample_dist(501:750,:,1) = exp(-0.5 * ((  sample_dist(501:750,:,1)  )./small_sigma).^2) ./ (sqrt(2*pi) .* small_sigma);
    sample_dist(501:750,:,2) = exp(-0.5 * ((  sample_dist(501:750,:,2)  )./large_sigma).^2) ./ (sqrt(2*pi) .* large_sigma);
    % Handle NaN values by replacing them with small values to avoid divide-by-zero errors
    sample_dist(isnan(sample_dist(:,1,1)),:,1) = 0;
    sample_dist(isnan(sample_dist(:,1,2)),:,2) = 0.001; 

    % Identify spikes that are poorly clustered (missed spikes)
    missed_spikes = find(sum(sample_dist(:,:,1),1)<1);
    if length(missed_spikes)>= 250  %If enough spikes outside concentrated areas, then random selection
        rand_sample(751:1000) = missed_spikes(randperm(length(missed_spikes),250));
    else
        [~,sort_spikes_i] = sort(sum(sample_dist(:,:,1),1),'ascend');  %Otherwise, sort by distances from sample
        rand_sample(751:1000) = sort_spikes_i(1:250);
    end
    % Calculate pairwise distances between samples in PCA space for new sample set
    sample_dist(751:1000,:,:) =  repmat(sqrt( sum( (permute(pca_score(rand_sample(751:1000),1:5),[1,3,2])-permute(pca_score(:,1:5),[3,1,2])).^2 ,3) ), [1,1,2]);
    % Filter out waveforms with fewer than 50 neighbors for clustering since these can't be cluster centers
    sample_dist([false(750,1);true(250,1)] & sum(sample_dist(:,:,1)<large_sigma,2)<50,:,:) = NaN;  %If the test waveform has less than 50 neighbors, it's not a cluster center
    sample_dist(751:1000,:,1) = exp(-0.5 * ((  sample_dist(751:1000,:,1)  )./small_sigma).^2) ./ (sqrt(2*pi) .* small_sigma);
    sample_dist(751:1000,:,2) = exp(-0.5 * ((  sample_dist(751:1000,:,2)  )./large_sigma).^2) ./ (sqrt(2*pi) .* large_sigma);
    % Handle NaN values by replacing them with small values to avoid divide-by-zero errors
    sample_dist(isnan(sample_dist(:,1,1)),:,1) = 0;
    sample_dist(isnan(sample_dist(:,1,2)),:,2) = 0.001; %Make small number rather than 0 to avoid divide by zero

    % Check for duplicate waveforms in the random sample
    [uniq_rand_sample, uniq_rand_sample_i] = unique(rand_sample);  
    if length(uniq_rand_sample)<1000
        rand_sample = rand_sample(uniq_rand_sample_i);
        sample_dist = sample_dist(uniq_rand_sample_i,:,:);
    end

     % Calculate the ratio of small to large distance weights
    sample_dist_ratio = sum(sample_dist(:,:,1),2)./sum(sample_dist(:,:,2),2);

    % Sort by distance ratios to find candidate clusters
    [sort_dist, sort_i] = sort(sample_dist_ratio, 'descend', 'MissingPlacement', 'last');
    sort_sample_dist = squareform(pdist(pca_score(rand_sample(sort_i),1:5)));
    sort_sample_dist = tril(sort_sample_dist);
    sort_sample_dist(sort_sample_dist==0) = Inf;
    sort_sample_dist_score = normpdf(sort_sample_dist,0,large_sigma)>1;

    % Identify candidate cluster units
    candidate_sort_dist = sort_dist.*double(~any(sort_sample_dist_score,2));
    test_template_units = find(candidate_sort_dist>0);

    % Iterate to refine cluster assignments based on overlap
    keep_going = 1;
    recalculate_overlap = 1;
    while keep_going
        if recalculate_overlap
            sample_dist_fract = sample_dist(sort_i(test_template_units),:,2)./sum(sample_dist(sort_i(test_template_units),:,2),1);
            overlap = (sample_dist_fract*sample_dist_fract')./ ...
                sqrt((sum(sample_dist_fract.^2,2)*sum(sample_dist_fract.^2,2)'));
            overlap(overlap>0.9999) = 0;
            overlap(1,1)=0;
            recalculate_overlap = 0;
        end

        % If there is overlap between clusters, resolve it
        if any(overlap(:)>0.05)
            [~,max_overlap_pair] = max(overlap(:));
            [unit1,unit2] = ind2sub(size(overlap), max_overlap_pair);
            [best_sample_dist,best_unit] = max(sample_dist(sort_i(test_template_units),:,2), [], 1);
            best_unit(best_sample_dist<0.1) = 0;
            curr_timestamps = [timestamps(best_unit==unit1); timestamps(best_unit==unit2)];
            curr_units = [ones(sum(best_unit==unit1),1); 2*ones(sum(best_unit==unit2),1)]; %1 for 1st unit, 2 for 2nd unit
            [sort_curr_timestamps, sort_ts_i] = sort(curr_timestamps);
            sort_curr_units = curr_units(sort_ts_i);
            sort_curr_units = sort_curr_units(1:(end-1)) + sort_curr_units(2:end);  %If ISI violation is one from unit 1 and one from unit 2, sum is 3
            
            % Calculate ISI violations
            num_ISI_violation12 = sum(diff(sort_curr_timestamps)<options.max_ISI/1000 & sort_curr_units==3);
            num_ISI_violation11 = sum(diff(sort_curr_timestamps)<options.max_ISI/1000 & sort_curr_units==2);
            num_ISI_violation22 = sum(diff(sort_curr_timestamps)<options.max_ISI/1000 & sort_curr_units==4);

            %Calculate percent true with the two units combined
            number_combined_violations = num_ISI_violation12+num_ISI_violation11+num_ISI_violation22;
            number_combined_waveforms = length(curr_units)*(tot_spikes/num_waveforms);
            if (.25 - number_combined_violations.*max(timestamps)/(2.*(options.max_ISI-options.min_ISI)/1000.*(number_combined_waveforms.^2))) > 0
                percent_same_combined = sqrt(.25 - number_combined_violations.*max(timestamps)./(2.*(options.max_ISI-options.min_ISI)/1000.*(number_combined_waveforms.^2)))+.5;
            else
                percent_same_combined = 0.5;
            end

             %Calculate percent true for unit 1
            number_combined_waveforms1 = sum(curr_units==1)*(tot_spikes/num_waveforms);
            if (.25 - num_ISI_violation11.*max(timestamps)/(2.*(options.max_ISI-options.min_ISI)/1000.*(number_combined_waveforms1.^2))) > 0
                percent_same_1 = sqrt(.25 - num_ISI_violation11.*max(timestamps)./(2.*(options.max_ISI-options.min_ISI)/1000.*(number_combined_waveforms1.^2)))+.5;
            else
                percent_same_1 = 0.5;
            end

             %Calculate percent true for unit 2
            number_combined_waveforms2 = sum(curr_units==2)*(tot_spikes/num_waveforms);
            if (.25 - num_ISI_violation22.*max(timestamps)/(2.*(options.max_ISI-options.min_ISI)/1000.*(number_combined_waveforms2.^2))) > 0
                percent_same_2 = real(sqrt(.25 - num_ISI_violation22.*max(timestamps)./(2.*(options.max_ISI-options.min_ISI)/1000.*(number_combined_waveforms2.^2)))+.5);
            else
                percent_same_2 = 0.5;
            end
            %If percent correct is better for the combined unit, make it a
            %single unit and remove the worse cluster center
            if percent_same_combined > (max([percent_same_1,percent_same_2])-0.15)
                if  sum(curr_units==1) < sum(curr_units==2) %If less unit 1 spikes, remove unit 1
                    test_template_units = test_template_units([1:(unit1-1),(unit1+1):end]);
                else  %Remove unit 2
                    test_template_units = test_template_units([1:(unit2-1),(unit2+1):end]);

                end
                recalculate_overlap = 1;
            else %else keep both units by saving as no overlap
                overlap(unit1,unit2) = 0;
                overlap(unit2,unit1) = 0;
            end
        else
            keep_going = 0;
        end
    end

    % Handle cases where the number of units exceeds the allowed maximum, default max is 2
    if length(test_template_units) > options.max_units_per_ch
        % Keep only the units with the most adjacent waveforms
        dist_cutoff = prctile(reshape(sample_dist(test_template_units,1:10:end,1),[],1), 100*(1-50/num_waveforms));
        spike_counts = gpuArray.zeros(length(test_template_units),1);
        for n = 1:length(test_template_units)
            spike_counts(n) = sum(sample_dist(sort_i(test_template_units(n)),:,1)>dist_cutoff );
        end
        [~,sort_spike_counts_i] = sort(spike_counts, 'descend');
        test_template_units = test_template_units(sort_spike_counts_i(1:options.max_units_per_ch));
    end

    % Generate template waveforms for the final clusters
    template_waveforms = zeros(size(waveforms,1),length(test_template_units));
    for n = 1:length(test_template_units)
        template_waveforms(:,n) = gather(mean(waveforms(:, sample_dist(sort_i(test_template_units(n)),:,1)>prctile(sample_dist(sort_i(test_template_units(n)),:,1), 100*(1-50/num_waveforms))),2));
    end

    %Return back template waveforms without the scale_factor
    template_waveforms = template_waveforms/scale_factor;
else
    % If fewer than 1000 waveforms, return an empty result
    template_waveforms = [];
end


