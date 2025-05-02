function WaveLimit(input_data_file, output_data_file, options,channels_to_sort, strobeInfo)
% Spike sorting function to process .nex, .plx, or .nev files and output
% a sorted .nex file with new sorting assignments.
%
% INPUT:
% input_data_file    - .nex, .plx, or .nev file to be sorted
% output_data_file   - .nex file to write new sorting assignments
%                      (Default: *_sorted.nex appended to input_data_file)
% options            - Spike sorting options (Default: from default_options function)
% channels_to_sort   - List of channels to be sorted (1-indexed)
%                      (Default: all channels)
% strobeInfo         - Trial start/end event codes, required to limit spike
%                      sorting to trial periods
%
% OUTPUT:
% Output .nex file with sorted spike waveforms.
%
% Adam Rouse, 9/30/24, v1.6
% Comments and formatting added by Xavier Scherschligt, 10/9/24

% Check if required function to write .nex files is in path
if ~exist('writeNexFile', 'file') || ~exist('writeNex5File', 'file')
    error('AGRCodeError:pathError', 'No writeNexFile.m found, add HowToReadAndWriteNexAndNex5FilesInMatlab to path')
end

% Detect file type based on file extension
if ~isempty(regexpi(input_data_file, '\.nex','once'))
    plx_file_flag = 0;
    nev_file_flag = 0;
elseif ~isempty(regexpi(input_data_file, '\.plx','once')) || ~isempty(regexpi(input_data_file, '\.pl2','once'))
    plx_file_flag = 1;
    nev_file_flag = 0;
elseif ~isempty(regexpi(input_data_file, '\.nev','once')) 
    plx_file_flag = 0;
    nev_file_flag = 1;
end

tic; % Start timer for performance tracking

% Default output file name if not specified
if ~exist('output_data_file', 'var') || isempty(output_data_file)
    if plx_file_flag
        output_data_file = [input_data_file(1:(regexpi(input_data_file, '\.plx')-1)), '_sorted.nex'];
    elseif nev_file_flag
        output_data_file = [input_data_file(1:(regexpi(input_data_file, '\.nev')-1)), '_sorted.nex'];
    else
        output_data_file = [input_data_file(1:(regexpi(input_data_file, '\.nex')-1)), '_sorted.nex'];
    end
end
% Default options if not provided
if ~exist('options', 'var') || isempty(options)
    options = default_options();
end

% Handle .plx files
if plx_file_flag
    % Check if Plexon SDK is in the path
    if ~exist('plx_info', 'file')
        error('AGRCodeError:pathError', 'No plx_info.m found, add Plexon Matlab Offline Files SDK to path')
    end
    
    % Load .plx file information
    [tscounts, ~, ~, ~] = plx_info(input_data_file, 1);
    unique_chan_numbers = find(sum(tscounts(:,2:end)>0)>1);
    [max_channels,channel_names] = plx_chan_names(input_data_file);
    
    % Load various file properties
    [~, ~,sample_freq, ~, ~, num_time_pts, ~, SpikePeakV, SpikeADResBits, ~, ~, FileDuration] = plx_information(input_data_file);
    [~,ChGains] = plx_chan_gains(input_data_file);
    
    % Load behavioral markers and correct signed integer codes
    [~, markers_ts, markers_sv] = plx_event_ts(input_data_file,257) ;  
    markers_sv_corrected = markers_sv + 32768;                          
    
    % Create a blank Nex file structure for output
    new_nexFileData = create_blank_nex();
    new_nexFileData.tend = FileDuration;
    new_nexFileData.markers{1}.name = 'events';
    new_nexFileData.markers{1}.values{1}.name = 'event_codes';
    new_nexFileData.markers{1}.values{1}.strings = arrayfun(@(x) num2str(x), markers_sv_corrected, 'Uni', false);
    new_nexFileData.markers{1}.timestamps = markers_ts;
    unit_counter = 1;
    total_time = FileDuration;
    
    % Include analog data if specified in options
    if options.keep_analog	
        [~, samplecounts] = plx_adchan_samplecounts(input_data_file);	
        [~, names] = plx_adchan_names(input_data_file);	
        %[n,gains] = plx_adchan_gains(input_data_file);	
        ch_with_data = find(samplecounts>0);	
        for ch = 1:length(ch_with_data)	
            % Read and store analog channel data
            [adfreq, ~, ts, fn, ad] = plx_ad_v(input_data_file, ch_with_data(ch)-1);	
            new_nexFileData.contvars{ch,1}.name = names(ch_with_data(ch),:);	
            new_nexFileData.contvars{ch,1}.varVersion = 200;	
            new_nexFileData.contvars{ch,1}.ADtoMV = 0.0024;  %Not sure how to calculate this from the gains in the plexon file so just use one default value for now	
            new_nexFileData.contvars{ch,1}.MVOfffset = 0;	
            new_nexFileData.contvars{ch,1}.ADFrequency = adfreq;	
            new_nexFileData.contvars{ch,1}.timestamps = ts;	
            new_nexFileData.contvars{ch,1}.fragmentStarts = [1 cumsum(fn(1:(end-1)))];	
            new_nexFileData.contvars{ch,1}.data = ad;	
        end	
    end	

% Handle .nev files
elseif nev_file_flag
    if ~exist('openNEV', 'file')
        error('AGRCodeError:pathError', 'No openNEV.m found, add NPMK-4.5.*.* to path')
    end
    nevFileData = openNEV(input_data_file, 'noread', 'nosave', 'nomat');  %There seem to be issues with event codes and reading from .mat file, for safety for now only read from original file and don't read data just to get events
    sample_freq = nevFileData.MetaTags.SampleRes;
    new_nexFileData = create_blank_nex();
    new_nexFileData.freq =  nevFileData.MetaTags.TimeRes;
    unit_counter = 1;
    new_nexFileData.markers{1}.name = 'events';
    new_nexFileData.markers{1}.varVarsion = 100;
    new_nexFileData.markers{1}.values{1}.name = 'event_codes';
    new_nexFileData.markers{1}.values{1}.strings = arrayfun(@(x) num2str(x), nevFileData.Data.SerialDigitalIO.UnparsedData, 'Uni', false);
    new_nexFileData.markers{1}.timestamps = nevFileData.Data.SerialDigitalIO.TimeStampSec';

    max_channels = max(nevFileData.Data.Spikes.Electrode);
    total_time = nevFileData.MetaTags.DataDurationSec;
    new_nexFileData.tend = total_time;
    new_nexFileData.contvars = {};

     % Reload NEV data for spike analysis
    nevFileData = openNEV(input_data_file, 'nosave', 'nomat'); 

% Handle .nex files
else
    if ~isempty(regexpi(input_data_file, '\.nex5','once'))
        nexFileData = readNex5File(input_data_file);
    else
        nexFileData = readNexFile(input_data_file);
    end
    
    % Nex file is zero-indexed, correct channel numbers
    nex_file_chan_numbers = cellfun(@(x) x.wireNumber, nexFileData.neurons) + 1; %Nex file is zero indexed
    nex_file_chanName_numbers = cellfun(@(x) str2double(x.name(regexp(x.name,'\d'))), nexFileData.neurons);
    nex_file_unit_numbers = cellfun(@(x) x.unitNumber, nexFileData.neurons);
    unique_chan_numbers = unique(nex_file_chan_numbers);
    num_time_pts = size(nexFileData.waves{1}.waveforms,1);
    
    % If strobeInfo not provided, use the full file duration
    if nargin < 5 %If no strobeInfo
        total_time = nexFileData.tend - nexFileData.tbeg;
    else
        % If strobeInfo, then only sort within trials, Find trial start and end samples
        Events = cellfun(@str2num,   nexFileData.markers{1}.values{1}.strings);
        EventTimeSamples = round(nexFileData.freq*nexFileData.markers{1}.timestamps);
        strobeInfo = find_trial_samples(strobeInfo, Events, EventTimeSamples);
        %Check if any spikes happen outside trial start and end
        if any( (sum((nexFileData.freq*nexFileData.waves{1}.timestamps(1:100:end)-strobeInfo.trial_start_samp') > 0, 2) - ...  %Number of trial starts before waveform timestamp (spot check only every 100 spike for time)
                sum((nexFileData.freq*nexFileData.waves{1}.timestamps(1:100:end)-strobeInfo.trial_end_samp') > 0, 2)) ~= 1)    %Number of trial ends before waveform timestamp %Difference between two numbers should be 1
            warning('Matlab:User', 'Spike times in file occur outside trial start and end times')
            total_time = nexFileData.tend - nexFileData.beg;
        else
             % Calculate total time based on trials
            total_samples = sum( strobeInfo.trial_end_samp - strobeInfo.trial_start_samp + 1  - nexFileData.waves{1}.NPointsWave );
            total_time = total_samples/nexFileData.freq;
        end
    end
    % Prepare the new Nex file structure
    new_nexFileData = nexFileData;
    new_nexFileData.neurons = {};
    new_nexFileData.waves = {};
    unit_counter = 1;
    orig_nexFileData_chans = cellfun(@(x) x.wireNumber, nexFileData.neurons)+1;
    max_channels = max(orig_nexFileData_chans);
end

if ~exist('channels_to_sort', 'var') || isempty(channels_to_sort)
    channels_to_sort = 1:max_channels;
else
    channels_to_sort = intersect(channels_to_sort,1:max_channels);
end


% Process each channel and sort waveforms
max_waveforms = 50000; %Only this many waveforms in sorting algorithm (the rest of the waveforms will be just assigned based on determined sort criteria)
for ch = 1:max_channels  % Loop over all channels
    % If dealing with Plexon or NEV files, data needs to be read to write to the .nex output file
    if plx_file_flag 
        curr_units = find(tscounts(:,ch+1)>0)-1;
        unit_counts = sum(tscounts(:,ch+1));
        waveforms = zeros(unit_counts,num_time_pts);
        timestamps = zeros(unit_counts,1);
        cluster_indexes = zeros(unit_counts,1);
        for u = 1:length(curr_units)
            % Get the range of spikes for the current unit
            curr_range = [sum(tscounts(1:curr_units(u),ch+1))+1,sum(tscounts(1:(curr_units(u)+1),ch+1))];
            % Read waveform and timestamp data from the file
            [~, ~, timestamps(curr_range(1):curr_range(2)), waveforms(curr_range(1):curr_range(2),:)] = plx_waves_v(input_data_file, ch, curr_units(u));
            % Assign the current unit index to the cluster
            cluster_indexes(curr_range(1):curr_range(2)) = (u-1)*ones(curr_range(2)-curr_range(1)+1,1);
        end
        % Transpose waveforms to match the format needed for sorting
        waveforms = waveforms';
        % Sort timestamps to ensure chronological order
        [timestamps, sort_i] = sort(timestamps);  
        waveforms = waveforms(:,sort_i); % Sort waveforms by the same order as timestamps
        
        % If we are using existing clusters, reorder the cluster indexes as well
        if options.use_existing_clusters
            cluster_indexes = cluster_indexes(sort_i);
        end
        curr_ch_name = channel_names(ch,:); % Set current channel name from channel names array
    elseif nev_file_flag  % If reading from NEV file
        % Get spikes for the current channel from NEV data
        curr_spikes_i = nevFileData.Data.Spikes.Electrode == ch;  % Get waveforms
        waveforms = double(nevFileData.Data.Spikes.Waveform(:,curr_spikes_i));
        unit_counts = size(waveforms,2); % Number of spikes
        timestamps = double(nevFileData.Data.Spikes.TimeStamp(curr_spikes_i))'/double(nevFileData.MetaTags.SampleRes);  % Get timestamps
        
        % If using existing clusters, get cluster indexes
        if options.use_existing_clusters
            cluster_indexes = nevFileData.Data.Spikes.Unit(curr_spikes_i);
        end
        curr_ch_name = ['sig', num2str(ch,'%03d')];
    else % Nex file handling
        ch_indexes = find(nex_file_chan_numbers == ch); % Get all indexes for current channel in the Nex file
        curr_ch_name = ['sig', num2str(nex_file_chanName_numbers(ch_indexes(1)),'%03d')]; % Set channel name
    end
    
    % Proceed if channel is to be sorted or if previous sorting is to be removed
    if (ismember(ch, channels_to_sort) || ~options.keep_unsorted_ch_assignments)  
        if ~plx_file_flag && ~nev_file_flag   % Handle Nex file waveforms and timestamps
            ch_indexes = find(nex_file_chan_numbers == ch); % Get channel indexes for Nex file
            if diff(nex_file_chanName_numbers(ch_indexes)) ~= 0 
                warning('Warning, the names (sig001U,a,b...) of the neurons with the same wireNumber do not have the same channel in their name!')
            end
            unit_counts = cellfun(@(x) size(x.waveforms,2),nexFileData.waves(ch_indexes));
            waveforms = zeros(num_time_pts,sum(unit_counts));
            timestamps = zeros(sum(unit_counts),1);

            % If using existing clusters, initialize cluster index array
            if options.use_existing_clusters
                cluster_indexes = zeros(sum(unit_counts),1);
            end

            % Load waveforms and timestamps for each unit on the current channel
            for u = 1:length(ch_indexes)
                curr_range = [sum(unit_counts(1:(u-1)))+1,sum(unit_counts(1:u))];
                waveforms(:,curr_range(1):curr_range(2)) = nexFileData.waves{ch_indexes(u)}.waveforms;
                timestamps(curr_range(1):curr_range(2)) = nexFileData.waves{ch_indexes(u)}.timestamps;
                
                 % Assign cluster indexes if using existing clusters
                if options.use_existing_clusters
                    cluster_indexes(curr_range(1):curr_range(2)) = nex_file_unit_numbers(ch_indexes(u))*ones(length(nexFileData.waves{ch_indexes(u)}.timestamps),1);
                end
            end
            
            % Sort timestamps and waveforms
            [timestamps, sort_i] = sort(timestamps);  
            waveforms = waveforms(:,sort_i);
            if options.use_existing_clusters
                cluster_indexes = cluster_indexes(sort_i);
            end
        end
        
        % Limit the number of waveforms to max_waveforms for generating
        % sort criteria (other waveforms will be assigned)
        if sum(unit_counts) > max_waveforms
            if options.repeatable_random_seed~=0
                rng(options.repeatable_random_seed); % Set random seed for reproducibility
            end
            ISI = diff(timestamps); % Calculate inter-spike intervals
            ISI_viol = 1000*ISI<options.max_ISI;  % Identify ISI violations
            ISI_viol = [false; ISI_viol] | [ISI_viol; false];  % Mark the violations
            num_ISI_viol = sum(ISI_viol); % Count the violations
            num_ISI_legit = sum(~ISI_viol);  % Count the legitimate spikes

            % If fewer than half of the max waveforms are ISI violations, include all violations
            % As the ISI violations are more informative about colliding
            % units include them all but need to correct for total number
            % of spikes
            if num_ISI_viol < (max_waveforms/2)  
                rand_i = find(~ISI_viol); % Select legitimate spikes
                rand_i = rand_i(randperm(length(rand_i),max_waveforms-num_ISI_viol)); % Randomly sample the remaining waveforms
                rand_i = sort([rand_i; find(ISI_viol)]); % Combine with violations
                tot_spikes = length(ISI_viol); % Total number of spikes (actual total with both the subsample and those not looked at)
            elseif num_ISI_legit < (max_waveforms/2)  % Handle the rare case where most spikes are ISI violations
                warning('Rare condition when num_ISI_legit<max_waveforms/2, while sum(unit_counts) > max_waveforms')
                ISI_viol_i = find(ISI_viol); % Get indexes of violations 
                ISI_viol_i = ISI_viol_i(randperm(length(ISI_viol_i), max_waveforms-num_ISI_legit)); % Sample from violations
                rand_i = sort([ISI_viol_i; find(~ISI_viol)]); % Combine with legitimate spikes
                tot_spikes = length(ISI_viol); % Total number of spikes (actual total with both the subsample and those not looked at)
            else % Balance the selection between legitimate and ISI violation spikes
                rand_i = find(~ISI_viol);
                if ~isempty(rand_i)  %Prevent error in very rare occurence when all ISIs are less than options.max_ISI
                    rand_i = rand_i(randperm(length(rand_i), floor(max_waveforms/2))); % Random sample of legitimate spikes
                    ISI_viol_i = find(ISI_viol);
                    ISI_viol_i = ISI_viol_i(randperm(length(ISI_viol_i), floor(max_waveforms/2)));
                    rand_i = sort([rand_i; ISI_viol_i]);
                    tot_spikes = floor((max_waveforms/2)*length(ISI_viol)/sum(ISI_viol)); % Maintain ISI violation proportion
                else
                    rand_i = sort(ISI_viol_i);  % In the rare case where all spikes are ISI violations
                end
            end
        else  % If fewer than max_waveforms, take all waveforms
            rand_i = 1:sum(unit_counts);
            tot_spikes = sum(unit_counts);
        end
        
    end

    % Proceed with spike sorting for the current channel if required
    if ismember(ch, channels_to_sort)
        display(['Sorting ch: ' num2str(ch)])
        if sum(unit_counts)>0  % Only proceed if there are waveforms
            % Align waveforms for sorting, either on GPU or CPU depending on options
        if options.align_waveforms
            if options.use_gpu &&  7 * (8*numel(waveforms)) < options.max_gpu_memory  % Prevent GPU memory overload
                aligned_waveforms = perform_spike_alignment_for_gpu(waveforms); % Use GPU for alignment
                aligned_waveforms = gather(aligned_waveforms); % Gather results from GPU
            else
                aligned_waveforms = perform_spike_alignment_for_cpu(waveforms); % Use CPU for alignment
            end
        else
            aligned_waveforms = waveforms; % Use original waveforms if no alignment is needed
        end
        
        % Cluster the waveforms based on existing clusters or create new clusters
        if options.use_existing_clusters
            if plx_file_flag
                unique_cluster_indexes = find(tscounts(:,1+ch)>0)-1; % Find unique cluster indexes
            elseif nev_file_flag
                unique_cluster_indexes = unique(cluster_indexes);
            else
                unique_cluster_indexes = nex_file_unit_numbers(ch_indexes);
            end
            unique_cluster_indexes = unique_cluster_indexes(unique_cluster_indexes~=0);
            mean_waveforms = zeros(size(aligned_waveforms,1),length(unique_cluster_indexes));

            % Calculate mean waveforms for each cluster
            for u = 1:length(unique_cluster_indexes)
                mean_waveforms(:,u) = mean(aligned_waveforms(:,cluster_indexes==unique_cluster_indexes(u)),2);
            end
        else % Perform clustering from scratch using GPU or CPU
            if options.use_gpu
                mean_waveforms = find_cluster_ids_for_gpu(aligned_waveforms(:,rand_i), timestamps(rand_i), options, tot_spikes);
            else
                mean_waveforms = find_cluster_ids_for_cpu(aligned_waveforms(:,rand_i), timestamps(rand_i), options, tot_spikes);
            end
        end
        else
            mean_waveforms = []; % No mean waveforms if no spikes
        end
    else
        mean_waveforms = []; % No sorting if the channel is not in the list to be sorted
    end
    
    % If cluster centers are found, proceed with clustering and post-processing
    if size(mean_waveforms,2)>0  
        if options.use_gpu % Perform clustering on GPU
            [test_mean_waveforms,cutoff_distances] = cluster_waveforms_for_gpu(aligned_waveforms(:,rand_i), 1000*timestamps(rand_i), mean_waveforms, options, 1000*total_time, tot_spikes);  %Call main clustering function
            test_mean_waveforms = gather(test_mean_waveforms);
        else % Perform clustering on CPU
            [test_mean_waveforms,cutoff_distances] = cluster_waveforms_for_cpu(aligned_waveforms(:,rand_i), 1000*timestamps(rand_i), mean_waveforms, options, 1000*total_time, tot_spikes);  %Call main clustering function
        end

        % Post-process clusters, include waveforms below the cutoff distance
        if size(test_mean_waveforms,2)>0 
            if ~options.keep_zero_clusters
                test_mean_waveforms = test_mean_waveforms(:,cutoff_distances>0);
                cutoff_distances = cutoff_distances(cutoff_distances>0);
            end
            
            % Calculate distances between waveforms and cluster centers
            % distances = sqrt(sum((repmat(permute(aligned_waveforms, [3 2 1]), [size(test_mean_waveforms,2),1,1]) - repmat(permute(test_mean_waveforms, [2 3 1]), [1,size(aligned_waveforms,2),1])).^2,3));
            distances = sum((repmat(permute(aligned_waveforms, [3 2 1]), [size(test_mean_waveforms,2),1,1]) - repmat(permute(test_mean_waveforms, [2 3 1]), [1,size(aligned_waveforms,2),1])).^2,3);
            pen_value = ceil(max(distances(:)));  % Penalty value for distance
            %Penalized distances adds pen_value to every distance that is not the
            %closest cluster center, for example if the distances for a spike
            %are .3, .4, .5 from the 3 cluster centers and pen_value = 1, the penalized
            %distances are .3, 1.4, 1.5
            pen_distances = distances;
            [~,min_d] = min(distances,[],1);
            for n = 1:size(distances,1)
                pen_distances(n,min_d~=n) = distances(n,min_d~=n)+pen_value;
            end
            % Sort distances and assign waveforms to clusters
            [sorted_distances, cluster_order] = sort(pen_distances,2); 
            waveform_assignments{ch} = zeros(size(aligned_waveforms,2),1);
            for n = 1:size(distances,1)
                waveform_assignments{ch}(cluster_order(n, sorted_distances(n,:)<=cutoff_distances(n))) = n;
            end
            
            % Perform cleanup by removing units that fail certain criteria (SNR, ISI violations, etc.)
            unique_cluster_indexes = unique(waveform_assignments{ch});
            for u = 1:length(unique_cluster_indexes)
                if unique_cluster_indexes(u)~=0
                    % Check signal-to-noise ratio (SNR)
                    if options.SNR_minimum>0
                        currSNR = signal2noise_PtoP2(aligned_waveforms(:,waveform_assignments{ch}==unique_cluster_indexes(u)));
                        if currSNR<options.SNR_minimum
                            waveform_assignments{ch}(waveform_assignments{ch}==unique_cluster_indexes(u)) = 0;
                            currSNR = -1;
                        end
                    else
                        currSNR = 100;  %Make large so comparison with options.SNR_minimum is always true
                    end

                    % Additional criteria: long intervals and line noise removal
                    if currSNR>=options.SNR_minimum && (options.remove_long_interval_units || options.remove_line_noise_units)
                        curr_diff_timestamps = diff(sort(timestamps(waveform_assignments{ch}==unique_cluster_indexes(u))));
                        % Expected number of short intervals (Poisson process model)
                        expected_not_long_intervals = length(curr_diff_timestamps)*(1-exp((-length(curr_diff_timestamps)/(1000*total_time))*options.long_interval_time));
                        % Remove long-interval units if needed
                        if options.remove_long_interval_units && (sum((1000*curr_diff_timestamps)<options.long_interval_time)/expected_not_long_intervals) < options.not_long_interval_fraction
                            waveform_assignments{ch}(waveform_assignments{ch}==unique_cluster_indexes(u)) = 0;
                        elseif options.remove_line_noise_units % Remove line noise-contaminated units
                            line_noise_amplitude = (abs(sum(exp(-1i*2*pi*options.line_noise_freq*curr_diff_timestamps))).^2)/length(curr_diff_timestamps);  %Rayleigh z statistic, for data with no 60Hz the mean value will be one, larger numbers more 60Hz ISI distribution
                            if line_noise_amplitude > options.line_noise_cutoff
                                waveform_assignments{ch}(waveform_assignments{ch}==unique_cluster_indexes(u)) = 0;
                            end
                        end
                        
                    end
                    
                end
            end

            % Remove clusters with no spike assignments
            if ~options.keep_zero_clusters  
                unique_cluster_indexes = sort(unique(waveform_assignments{ch}));
                unique_cluster_indexes = unique_cluster_indexes(unique_cluster_indexes>0); %Temporarily remove 0-unsorted units
                for u = 1:length(unique_cluster_indexes)
                    waveform_assignments{ch}(waveform_assignments{ch}==unique_cluster_indexes(u)) = u;
                end
                unique_cluster_indexes = unique(waveform_assignments{ch});
            end
            
                % Assign waveforms to units and save them in the output Nex file
                for u = 1:length(unique_cluster_indexes)
                    if unique_cluster_indexes(u) == 0
                        new_nexFileData.neurons{unit_counter,1}.name = [curr_ch_name, 'U'];
                    else
                        new_nexFileData.neurons{unit_counter,1}.name = [curr_ch_name, char(96+unique_cluster_indexes(u))];
                    end
                  
                    new_nexFileData.neurons{unit_counter,1}.varVersion = 101;
                    new_nexFileData.neurons{unit_counter,1}.wireNumber = ch-1;
                    new_nexFileData.neurons{unit_counter,1}.unitNumber = unique_cluster_indexes(u);
                    new_nexFileData.neurons{unit_counter,1}.xPos = 0;
                    new_nexFileData.neurons{unit_counter,1}.yPos = 0;
                    new_nexFileData.neurons{unit_counter,1}.timestamps = timestamps(waveform_assignments{ch}==unique_cluster_indexes(u));
                    
                    
                    new_nexFileData.waves{unit_counter,1}.name = [new_nexFileData.neurons{unit_counter,1}.name, '_wf'];
                    new_nexFileData.waves{unit_counter,1}.varVersion = new_nexFileData.neurons{unit_counter,1}.varVersion;
                    new_nexFileData.waves{unit_counter,1}.NPointsWave = size(waveforms,1);
                    
                    new_nexFileData.waves{unit_counter,1}.wireNumber = new_nexFileData.neurons{unit_counter,1}.wireNumber;
                    new_nexFileData.waves{unit_counter,1}.unitNumber = new_nexFileData.neurons{unit_counter,1}.unitNumber;
                    new_nexFileData.waves{unit_counter,1}.MVOffset = 0;
                    new_nexFileData.waves{unit_counter,1}.timestamps = new_nexFileData.neurons{unit_counter,1}.timestamps;
                    new_nexFileData.waves{unit_counter,1}.waveforms = waveforms(:,waveform_assignments{ch}==unique_cluster_indexes(u));
                    
                    % Set ADtoMV conversion factor based on file type
                    if plx_file_flag
                        new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                        new_nexFileData.waves{unit_counter,1}.ADtoMV = SpikePeakV/((2^(SpikeADResBits-1))*ChGains(ch));
                    elseif nev_file_flag 
                        new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                        new_nexFileData.waves{unit_counter,1}.ADtoMV = 6e-5;
                        new_nexFileData.waves{unit_counter,1}.waveforms = new_nexFileData.waves{unit_counter,1}.waveforms.*new_nexFileData.waves{unit_counter,1}.ADtoMV;
                    else
                        new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{ch_indexes(1),1}.WFrequency;
                        new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{ch_indexes(1),1}.ADtoMV;
                    end
                    unit_counter = unit_counter+1;
                end
            
        else
            % If no waveforms, assign the entire channel to "unsorted"
            new_nexFileData.neurons{unit_counter,1}.name = [curr_ch_name, 'U'];
            new_nexFileData.neurons{unit_counter,1}.varVersion = 101;
            new_nexFileData.neurons{unit_counter,1}.wireNumber = ch-1;
            new_nexFileData.neurons{unit_counter,1}.unitNumber = 0;
            new_nexFileData.neurons{unit_counter,1}.xPos = 0;
            new_nexFileData.neurons{unit_counter,1}.yPos = 0;
            new_nexFileData.neurons{unit_counter,1}.timestamps = timestamps;
            
            % Save unsorted waveforms
            new_nexFileData.waves{unit_counter,1}.name = [new_nexFileData.neurons{unit_counter,1}.name, '_wf'];
            new_nexFileData.waves{unit_counter,1}.varVersion = new_nexFileData.neurons{unit_counter,1}.varVersion;
            new_nexFileData.waves{unit_counter,1}.NPointsWave = size(waveforms,1);
            new_nexFileData.waves{unit_counter,1}.wireNumber = new_nexFileData.neurons{unit_counter,1}.wireNumber;
            new_nexFileData.waves{unit_counter,1}.unitNumber = new_nexFileData.neurons{unit_counter,1}.unitNumber;
            
            new_nexFileData.waves{unit_counter,1}.MVOffset = 0;
            new_nexFileData.waves{unit_counter,1}.timestamps = new_nexFileData.neurons{unit_counter,1}.timestamps;
            new_nexFileData.waves{unit_counter,1}.waveforms = waveforms;
            
            % Set ADtoMV conversion based on file type
            if ~isempty(ch_indexes)
                if plx_file_flag
                    new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = SpikePeakV/((2^(SpikeADResBits-1))*ChGains(ch));
                elseif nev_file_flag
                    new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = 6e-5;
                    new_nexFileData.waves{unit_counter,1}.waveforms = new_nexFileData.waves{unit_counter,1}.waveforms.*new_nexFileData.waves{unit_counter,1}.ADtoMV;
                else
                    new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{ch_indexes(1),1}.WFrequency;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{ch_indexes(1),1}.ADtoMV;
                end
            else %If there's no waveforms on the channel, just use 1st channel sampling frequency and gain
                new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{1,1}.WFrequency;
                new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{1,1}.ADtoMV;
            end
            
            unit_counter = unit_counter+1;
        end
        
        
    else   %if size(mean_waveforms,2)>0, no waveforms, don't assign new sortings
     
        if options.keep_unsorted_ch_assignments
            if plx_file_flag
                curr_orig_units = unique(cluster_indexes);
                waveform_assignments{ch} = cluster_indexes;
                for u = curr_orig_units'
                    if u == 0
                        new_nexFileData.neurons{unit_counter,1}.name = [curr_ch_name, 'U'];
                    else
                        new_nexFileData.neurons{unit_counter,1}.name = [curr_ch_name, char(96+u)];
                    end
                    new_nexFileData.neurons{unit_counter,1}.varVersion = 101;
                    new_nexFileData.neurons{unit_counter,1}.wireNumber = ch-1;
                    new_nexFileData.neurons{unit_counter,1}.unitNumber = u;
                    new_nexFileData.neurons{unit_counter,1}.xPos = 0;
                    new_nexFileData.neurons{unit_counter,1}.yPos = 0;
                    new_nexFileData.neurons{unit_counter,1}.timestamps = timestamps(cluster_indexes==u);
                    
                    new_nexFileData.waves{unit_counter,1}.name = [new_nexFileData.neurons{unit_counter,1}.name, '_wf'];
                    new_nexFileData.waves{unit_counter,1}.varVersion = new_nexFileData.neurons{unit_counter,1}.varVersion;
                    new_nexFileData.waves{unit_counter,1}.NPointsWave = size(waveforms,1);
                    
                    new_nexFileData.waves{unit_counter,1}.wireNumber = new_nexFileData.neurons{unit_counter,1}.wireNumber;
                    new_nexFileData.waves{unit_counter,1}.unitNumber = new_nexFileData.neurons{unit_counter,1}.unitNumber;
                    new_nexFileData.waves{unit_counter,1}.MVOffset = 0;
                    new_nexFileData.waves{unit_counter,1}.timestamps = new_nexFileData.neurons{unit_counter,1}.timestamps;
                    new_nexFileData.waves{unit_counter,1}.waveforms = waveforms(:,waveform_assignments{ch}==u);
                    new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = SpikePeakV/((2^(SpikeADResBits-1))*ChGains(ch));
                    
                    unit_counter = unit_counter+1;
                end
            elseif nev_file_flag
                curr_ch_spikes = nevFileData.Data.Spikes.Electrode==ch;
                curr_orig_units = unique(nevFileData.Data.Spikes.Unit(curr_ch_spikes));
                for u = curr_orig_units
                    curr_spikes_index =  curr_ch_spikes & nevFileData.Data.Spikes.Unit == u;
                    if u == 0
                        new_nexFileData.neurons{unit_counter,1}.name = [curr_ch_name, 'U'];
                    else
                        new_nexFileData.neurons{unit_counter,1}.name = [curr_ch_name, char(96+u)];
                    end
                    new_nexFileData.neurons{unit_counter,1}.varVersion = 101;
                    new_nexFileData.neurons{unit_counter,1}.wireNumber = ch-1;
                    new_nexFileData.neurons{unit_counter,1}.unitNumber = u;
                    new_nexFileData.neurons{unit_counter,1}.xPos = 0;
                    new_nexFileData.neurons{unit_counter,1}.yPos = 0;
                    new_nexFileData.neurons{unit_counter,1}.timestamps = double(nevFileData.Data.Spikes.TimeStamp(curr_spikes_i))'/double(nevFileData.MetaTags.SampleRes);   
                    
                    new_nexFileData.waves{unit_counter,1}.name = [new_nexFileData.neurons{unit_counter,1}.name, '_wf'];
                    new_nexFileData.waves{unit_counter,1}.varVersion = new_nexFileData.neurons{unit_counter,1}.varVersion;
                    new_nexFileData.waves{unit_counter,1}.NPointsWave = size(waveforms,1);
                    
                    new_nexFileData.waves{unit_counter,1}.wireNumber = new_nexFileData.neurons{unit_counter,1}.wireNumber;
                    new_nexFileData.waves{unit_counter,1}.unitNumber = new_nexFileData.neurons{unit_counter,1}.unitNumber;
                    new_nexFileData.waves{unit_counter,1}.MVOffset = 0;
                    new_nexFileData.waves{unit_counter,1}.timestamps = new_nexFileData.neurons{unit_counter,1}.timestamps;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = 6e-05;
                    new_nexFileData.waves{unit_counter,1}.waveforms = double(nevFileData.Data.Spikes.Waveform(:,curr_spikes_i)).*new_nexFileData.waves{unit_counter,1}.ADtoMV;
                    new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                    
                    
                    unit_counter = unit_counter+1;
                end
            else
                curr_orig_units = find(orig_nexFileData_chans==ch);
                for u = curr_orig_units'
                    new_nexFileData.neurons{unit_counter,1} = nexFileData.neurons{u,1};
                    new_nexFileData.waves{unit_counter,1} = nexFileData.waves{u,1};
                    unit_counter = unit_counter+1;
                end
            end
        else  %All waveforms assigned to unsorted
            new_nexFileData.neurons{unit_counter,1}.name = [curr_ch_name, 'U']; 
            new_nexFileData.neurons{unit_counter,1}.varVersion = 101;
            new_nexFileData.neurons{unit_counter,1}.wireNumber = ch-1;  %nex_file_chanName_numbers(curr_orig_units(1))-1;
            new_nexFileData.neurons{unit_counter,1}.unitNumber = 0;
            new_nexFileData.neurons{unit_counter,1}.xPos = 0;
            new_nexFileData.neurons{unit_counter,1}.yPos = 0;
            new_nexFileData.neurons{unit_counter,1}.timestamps = timestamps;
            
            new_nexFileData.waves{unit_counter,1}.name = [new_nexFileData.neurons{unit_counter,1}.name, '_wf'];
            new_nexFileData.waves{unit_counter,1}.varVersion = new_nexFileData.neurons{unit_counter,1}.varVersion;
            new_nexFileData.waves{unit_counter,1}.NPointsWave = size(waveforms,1);
            new_nexFileData.waves{unit_counter,1}.wireNumber = new_nexFileData.neurons{unit_counter,1}.wireNumber;
            new_nexFileData.waves{unit_counter,1}.unitNumber = new_nexFileData.neurons{unit_counter,1}.unitNumber;
            
            new_nexFileData.waves{unit_counter,1}.MVOffset = 0;
            new_nexFileData.waves{unit_counter,1}.timestamps = new_nexFileData.neurons{unit_counter,1}.timestamps;
            new_nexFileData.waves{unit_counter,1}.waveforms = waveforms;
            if plx_file_flag
                new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                new_nexFileData.waves{unit_counter,1}.ADtoMV = SpikePeakV/((2^(SpikeADResBits-1))*ChGains(ch));
            elseif nev_file_flag
                new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                new_nexFileData.waves{unit_counter,1}.ADtoMV = 6e-05;
                new_nexFileData.waves{unit_counter,1}.waveforms = new_nexFileData.waves{unit_counter,1}.waveforms.*new_nexFileData.waves{unit_counter,1}.ADtoMV;
            else
                if ~isempty(ch_indexes)
                    new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{ch_indexes(1),1}.WFrequency;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{ch_indexes(1),1}.ADtoMV;
                    
                else %If there's no waveforms on the channel, just use 1st channel sampling frequency and gain
                    new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{1,1}.WFrequency;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{1,1}.ADtoMV;
                end
            end
            
            unit_counter = unit_counter+1;
        end
        
    end
end


for k = 1:length(new_nexFileData.neurons)
    % Pad names with null characters for Plexon offline sorter compatibility
    new_nexFileData.neurons{k}.name((end+1):64) = 0;  
    new_nexFileData.waves{k}.name((end+1):64) = 0;
end

% Save sorted data to .nex or .nex5 file
if strcmp(output_data_file(end-4:end), '.nex5')  % save to nex5 file if assigned to, ZL 2019-11-24
    writeNex5File(new_nexFileData, output_data_file);
else
    writeNexFile(new_nexFileData, output_data_file);
    s = dir(output_data_file);
    if s.bytes/1024^3 >= 2  % save to .nex5 file if output .nex file >2GB
        warning('output sorted .nex file greater than 2GB, data for last few channels unreliable')
        output_data_file = [output_data_file '5'];
        disp(['Save a .nex5 file: ' output_data_file] )
        writeNex5File(new_nexFileData, output_data_file);
    end
end

toc % End timer
end

% Helper function to calculate signal-to-noise ratio
function SNR = signal2noise_PtoP2(waves)
%Peak to Peak signal-to-noise divided by 2
mean_signal = mean(waves,2);
signal = max(mean_signal) - min(mean_signal);
remain_noise = waves-repmat(mean_signal, [1, size(waves,2)]);
noise = std(remain_noise(:));
SNR = (signal./noise)/2;
end
