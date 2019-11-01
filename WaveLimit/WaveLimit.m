function WaveLimit(input_data_file,output_data_file,options,channels_to_sort, strobeInfo)
%WaveLimit(input_data_file,output_data_file,options,channels_to_sort)
%channels_to_sort channel numbers that should be sorted - NOTE Must be 1
%indexed from channels 1 to n NOT 0 to n-1
%
% input_data_file - .nex or .plx file to be sorted
% output_data_file - .nex file to write new sorting assignments (Default, add *_sorted.nex to input_data_file)
% options - spike sorting options (Default, options from function: default_options)
% channels_to_sort - list of channels to be sorted (1 indexed) (Default, all channels)
% strobeInfo - trial start and end event codes, necessary if only spikes during trials is included to accurately estimate the total time of the file
% Adam Rouse 10/25/19, v1.1


% TODO if only 1 waveform then perform_spike_align_cpu crashes


if ~exist('writeNexFile', 'file')
    error('AGRCodeError:pathError', 'No writeNexFile.m found, add HowToReadAndWriteNexAndNex5FilesInMatlab to path')
end


if ~isempty(regexpi(input_data_file, '\.nex'))
    plx_file_flag = 0;
elseif ~isempty(regexpi(input_data_file, '\.plx')) || ~isempty(regexpi(input_data_file, '\.pl2'))
    plx_file_flag = 1;
end

tic
if ~exist('output_data_file', 'var') || isempty(output_data_file)
    if plx_file_flag
        output_data_file = [input_data_file(1:(regexpi(input_data_file, '\.plx')-1)), '_sorted.nex'];
    else
        output_data_file = [input_data_file(1:(regexpi(input_data_file, '\.nex')-1)), '_sorted.nex'];
    end
end
if ~exist('options', 'var') || isempty(options)
    options = default_options();
end

if plx_file_flag
    if ~exist('plx_info', 'file')
        error('AGRCodeError:pathError', 'No plx_info.m found, add Plexon Matlab Offline Files SDK to path')
    end
    [tscounts, ~, ~, ~] = plx_info(input_data_file, 1);
    unique_chan_numbers = find(sum(tscounts(:,2:end)>0)>1);
    [max_channels,channel_names] = plx_chan_names(input_data_file);
    
    [~, ~,sample_freq, ~, ~, num_time_pts, ~, SpikePeakV, SpikeADResBits, ~, ~, FileDuration] = plx_information(input_data_file);
    [~,ChGains] = plx_chan_gains(input_data_file);
    [~, markers_ts, markers_sv] = plx_event_ts(input_data_file,257) ;  %read behavioral event marker codes into sv, timestamps into ts
    markers_sv_corrected = markers_sv + 32768;                          %correct signed integer to integer
    new_nexFileData = create_blank_nex();
    new_nexFileData.tend = FileDuration;
    new_nexFileData.markers{1}.name = 'events';
    new_nexFileData.markers{1}.values{1}.name = 'event_codes';
    new_nexFileData.markers{1}.values{1}.strings = arrayfun(@(x) num2str(x), markers_sv_corrected, 'Uni', false);
    new_nexFileData.markers{1}.timestamps = markers_ts;
    unit_counter = 1;
    total_time = FileDuration;
else
    nexFileData = readNexFile(input_data_file);
    nex_file_chan_numbers = cellfun(@(x) x.wireNumber, nexFileData.neurons) + 1; %Nex file is zero indexed
    nex_file_unit_numbers = cellfun(@(x) x.unitNumber, nexFileData.neurons);
    unique_chan_numbers = unique(nex_file_chan_numbers);
    num_time_pts = size(nexFileData.waves{1}.waveforms,1);
    if nargin < 5 %If no strobeInfo
        total_time = nexFileData.tend - nexFileData.tbeg;
    else
        %Find trial start and end samples using find_trial_sampes
        Events = cellfun(@str2num,   nexFileData.markers{1}.values{1}.strings);
        EventTimeSamples = round(nexFileData.freq*nexFileData.markers{1}.timestamps);
        strobeInfo = find_trial_samples(strobeInfo, Events, EventTimeSamples);
        %Check if any spikes happen outside trial start and end
        if any( (sum((nexFileData.freq*nexFileData.waves{1}.timestamps(1:100:end)-strobeInfo.trial_start_samp') > 0, 2) - ...  %Number of trial starts before waveform timestamp (spot check only every 100 spike for time)
                sum((nexFileData.freq*nexFileData.waves{1}.timestamps(1:100:end)-strobeInfo.trial_end_samp') > 0, 2)) ~= 1)    %Number of trial ends before waveform timestamp %Difference between two numbers should be 1
            warning('Matlab:User', 'Spike times in file occur outside trial start and end times')
            total_time = nexFileData.tend - nexFileData.beg;
        else
            %If spikes do indeed fall only within trial start and ends, then the total time is the total trial times
            total_samples = sum( strobeInfo.trial_end_samp - strobeInfo.trial_start_samp + 1  - nexFileData.waves{1}.NPointsWave );
            total_time = total_samples/nexFileData.freq;
        end
    end
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
if plx_file_flag
    
else
    %     zero_unit_channels = zeros(size(nex_file_chan_numbers));  %Will label which units had zero waveforms assigned
end

max_waveforms = 100000;
for ch = 1:max_channels   %length(channels_to_sort)
    
    if plx_file_flag  %If we're reading in a plexon file, then we always have to read in the data to write to nex
        curr_units = find(tscounts(:,ch+1)>0)-1;
            unit_counts = sum(tscounts(:,ch+1));
            waveforms = zeros(unit_counts,num_time_pts);
            timestamps = zeros(unit_counts,1);
            cluster_indexes = zeros(unit_counts,1);
            for u = 1:length(curr_units)
                curr_range = [sum(tscounts(1:curr_units(u),ch+1))+1,sum(tscounts(1:(curr_units(u)+1),ch+1))];
                [~, ~, timestamps(curr_range(1):curr_range(2)), waveforms(curr_range(1):curr_range(2),:)] = plx_waves_v(input_data_file, ch, curr_units(u));
                cluster_indexes(curr_range(1):curr_range(2)) = (u-1)*ones(curr_range(2)-curr_range(1)+1,1);
            end
            waveforms = waveforms';
            [timestamps, sort_i] = sort(timestamps);  %Make sure time stamps are sorted
            waveforms = waveforms(:,sort_i);
            if options.use_existing_clusters
                cluster_indexes = cluster_indexes(sort_i);
            end
    end
    
    if (ismember(ch, channels_to_sort) || ~options.keep_unsorted_ch_assignments)  %If it's a channel to sort or remove previous sorting, then read in waveforms
        if ~plx_file_flag  %If it is a nex file, read in data for sorting
            ch_indexes = find(nex_file_chan_numbers == ch);
            unit_counts = cellfun(@(x) size(x.waveforms,2),nexFileData.waves(ch_indexes));
            %         if sum(unit_counts)<= max_waveforms
            waveforms = zeros(num_time_pts,sum(unit_counts));
            timestamps = zeros(sum(unit_counts),1);
            if options.use_existing_clusters
                cluster_indexes = zeros(sum(unit_counts),1);
            end
            for u = 1:length(ch_indexes)
                curr_range = [sum(unit_counts(1:(u-1)))+1,sum(unit_counts(1:u))];
                waveforms(:,curr_range(1):curr_range(2)) = nexFileData.waves{ch_indexes(u)}.waveforms;
                timestamps(curr_range(1):curr_range(2)) = nexFileData.waves{ch_indexes(u)}.timestamps;
                if options.use_existing_clusters
                    cluster_indexes(curr_range(1):curr_range(2)) = nex_file_unit_numbers(ch_indexes(u))*ones(length(nexFileData.waves{ch_indexes(u)}.timestamps),1);
                end
            end
            
            [timestamps, sort_i] = sort(timestamps);  %Make sure time stamps are sorted
            waveforms = waveforms(:,sort_i);
            if options.use_existing_clusters
                cluster_indexes = cluster_indexes(sort_i);
            end
        end
        
        if sum(unit_counts) > max_waveforms
            if options.repeatable_random_seed~=0
                rng(options.repeatable_random_seed);
            end
            ISI = diff(timestamps);
            ISI_viol = 1000*ISI<options.max_ISI;
            ISI_viol = [false; ISI_viol] | [ISI_viol; false];
            num_ISI_viol = sum(ISI_viol);
            if num_ISI_viol < (max_waveforms/2)
                rand_i = find(~ISI_viol);
                rand_i = rand_i(randperm(length(rand_i),max_waveforms-num_ISI_viol));
                rand_i = sort([rand_i; find(ISI_viol)]);
                tot_spikes = length(ISI_viol);
            else
                rand_i = find(~ISI_viol);
                if ~isempty(rand_i)  %Prevent error in very rare occurence when all ISIs are less than options.max_ISI
                rand_i = rand_i(randperm(length(rand_i), floor(max_waveforms/2)));
                ISI_viol_i = find(ISI_viol);
                ISI_viol_i = ISI_viol_i(randperm(length(ISI_viol_i), floor(max_waveforms/2)));
                rand_i = sort([rand_i; ISI_viol_i]);
                tot_spikes = floor((max_waveforms/2)*length(ISI_viol)/sum(ISI_viol)); %Not the actual total number of spikes but keeps ratio of included ISI violations to tot_spikes the same
                else
                    rand_i = sort(ISI_viol_i);
                end
            end
        else
            rand_i = 1:sum(unit_counts);
            tot_spikes = sum(unit_counts);
        end
        
    end
    if ismember(ch, channels_to_sort)
        display(['Sorting ch: ' num2str(ch)])
        if options.align_waveforms
            if options.use_gpu
                aligned_waveforms = perform_spike_alignment_for_gpu(waveforms);
                aligned_waveforms = gather(aligned_waveforms);
            else
                aligned_waveforms = perform_spike_alignment_for_cpu(waveforms);
            end
        else
            aligned_waveforms = waveforms;
        end
        
        if options.use_existing_clusters
            unique_cluster_indexes = nex_file_unit_numbers(ch_indexes);
            unique_cluster_indexes = unique_cluster_indexes(unique_cluster_indexes~=0);
            mean_waveforms = zeros(size(aligned_waveforms,1),length(unique_cluster_indexes));
            for u = 1:length(unique_cluster_indexes)
                mean_waveforms(:,u) = mean(aligned_waveforms(:,cluster_indexes==unique_cluster_indexes(u)),2);
            end
        else
            if options.use_gpu
                mean_waveforms = find_cluster_ids_for_gpu(aligned_waveforms(:,rand_i), timestamps(rand_i), options, tot_spikes);
            else
                mean_waveforms = find_cluster_ids_for_cpu(aligned_waveforms(:,rand_i), timestamps(rand_i), options, tot_spikes);
            end
        end
    else
        mean_waveforms = [];
    end
    
    
    if size(mean_waveforms,2)>0  %If there were cluster centers found
        if options.use_gpu
            [test_mean_waveforms,cutoff_distances] = cluster_waveforms_for_gpu(aligned_waveforms(:,rand_i), 1000*timestamps(rand_i), mean_waveforms, options, 1000*total_time, tot_spikes);  %Call main clustering function
        else
            [test_mean_waveforms,cutoff_distances] = cluster_waveforms_for_cpu(aligned_waveforms(:,rand_i), 1000*timestamps(rand_i), mean_waveforms, options, 1000*total_time, tot_spikes); new_nexFileData %Call main clustering function
        end
        if size(test_mean_waveforms,2)>0  %If there were cluster centers found
            if ~options.keep_zero_clusters
                test_mean_waveforms = test_mean_waveforms(:,cutoff_distances>0);
                cutoff_distances = cutoff_distances(cutoff_distances>0);
            end
            
            distances = sqrt(sum((repmat(permute(aligned_waveforms, [3 2 1]), [size(test_mean_waveforms,2),1,1]) - repmat(permute(test_mean_waveforms, [2 3 1]), [1,size(aligned_waveforms,2),1])).^2,3));
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
            [sorted_distances, cluster_order] = sort(pen_distances,2);  %Sort the pen_distances to order waveforms by their distance from the test spike template
            waveform_assignments{ch} = zeros(size(aligned_waveforms,2),1);
            for n = 1:size(distances,1)
                waveform_assignments{ch}(cluster_order(n, sorted_distances(n,:)<=cutoff_distances(n))) = n;
            end
            
            %Perform clean-up by removing units that violate the various
            %criteria for being considered a neuron and not just noise
            %contaminiation
            unique_cluster_indexes = unique(waveform_assignments{ch});
            for u = 1:length(unique_cluster_indexes)
                if unique_cluster_indexes(u)~=0
                    if options.SNR_minimum>0
                        currSNR = signal2noise_PtoP2(aligned_waveforms(:,waveform_assignments{ch}==unique_cluster_indexes(u)));
                        if currSNR<options.SNR_minimum
                            waveform_assignments{ch}(waveform_assignments{ch}==unique_cluster_indexes(u)) = 0;
                            currSNR = -1;
                        end
                    end
                    if currSNR>=options.SNR_minimum && (options.remove_long_interval_units || options.remove_line_noise_units)
                        curr_diff_timestamps = diff(sort(timestamps(waveform_assignments{ch}==unique_cluster_indexes(u))));
                        %calculate expect number of intervals less than the
                        %long interval time, using exponential distribution
                        %to model Poisson process, for end of trial
                        %reward artifacts and similar
                        expected_not_long_intervals = length(curr_diff_timestamps)*(1-exp((-length(curr_diff_timestamps)/(1000*total_time))*options.long_interval_time));
                        if options.remove_long_interval_units && (sum((1000*curr_diff_timestamps)<options.long_interval_time)/expected_not_long_intervals) < options.not_long_interval_fraction
                            waveform_assignments{ch}(waveform_assignments{ch}==unique_cluster_indexes(u)) = 0;
                            %calculate power-line noise (60Hz in US) to remove
                            %units that have excess 60Hz rhythmicity
                        elseif options.remove_line_noise_units
                            line_noise_amplitude = (abs(sum(exp(-1i*2*pi*options.line_noise_freq*curr_diff_timestamps))).^2)/length(curr_diff_timestamps);  %Rayleigh z statistic, for data with no 60Hz the mean value will be one, larger numbers more 60Hz ISI distribution
                            if line_noise_amplitude > options.line_noise_cutoff
                                waveform_assignments{ch}(waveform_assignments{ch}==unique_cluster_indexes(u)) = 0;
                            end
                        end
                        
                    end
                    
                end
            end
            if ~options.keep_zero_clusters  %Remove units with no spikes assigned
                unique_cluster_indexes = sort(unique(waveform_assignments{ch}));
                unique_cluster_indexes = unique_cluster_indexes(unique_cluster_indexes>0); %Temporarily remove 0-unsorted units
                for u = 1:length(unique_cluster_indexes)
                    waveform_assignments{ch}(waveform_assignments{ch}==unique_cluster_indexes(u)) = u;
                end
                unique_cluster_indexes = unique(waveform_assignments{ch});
            end
            
                
                for u = 1:length(unique_cluster_indexes)
                    if plx_file_flag
                        if unique_cluster_indexes(u) == 0
                            new_nexFileData.neurons{unit_counter,1}.name = [channel_names(ch,:), 'U'];
                        else
                            new_nexFileData.neurons{unit_counter,1}.name = [channel_names(ch,:), char(96+unique_cluster_indexes(u))];
                        end
                    else
                        if unique_cluster_indexes(u) == 0
                            new_nexFileData.neurons{unit_counter,1}.name = ['sig', num2str(ch,'%03d'), 'U'];
                        else
                            new_nexFileData.neurons{unit_counter,1}.name = ['sig', num2str(ch,'%03d'), char(96+unique_cluster_indexes(u))];
                        end
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
                    if plx_file_flag
                        new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                        new_nexFileData.waves{unit_counter,1}.ADtoMV = SpikePeakV/((2^(SpikeADResBits-1))*ChGains(ch));
                    else
                        new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{ch_indexes(1),1}.WFrequency;
                        new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{ch_indexes(1),1}.ADtoMV;
                    end
                    unit_counter = unit_counter+1;
                end
            
        else
            %All waveforms assigned to unsorted
            new_nexFileData.neurons{unit_counter,1}.name = ['sig', num2str(ch,'%03d'), 'U'];
            new_nexFileData.neurons{unit_counter,1}.varVersion = 101;
            new_nexFileData.neurons{unit_counter,1}.wireNumber = ch-1;
            new_nexFileData.neurons{unit_counter,1}.unitNumber = 0;
            new_nexFileData.neurons{unit_counter,1}.xPos = 0;
            new_nexFileData.neurons{unit_counter,1}.yPos = 0;
            new_nexFileData.neurons{unit_counter,1}.timestamps = timestamps;
            
            new_nexFileData.waves{unit_counter,1}.name = [new_nexFileData.neurons{unit_counter,1}.name, '_wf'];
            new_nexFileData.waves{unit_counter,1}.varVersion = new_nexFileData.neurons{unit_counter,1}.varVersion;
            new_nexFileData.waves{unit_counter,1}.NPointsWave = size(waveforms,1);
            if ~isempty(ch_indexes)
                if plx_file_flag
                    new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = SpikePeakV/((2^(SpikeADResBits-1))*ChGains(ch));
                else
                    
                    new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{ch_indexes(1),1}.WFrequency;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{ch_indexes(1),1}.ADtoMV;
                end
            else %If there's no waveforms on the channel, just use 1st channel sampling frequency and gain
                new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{1,1}.WFrequency;
                new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{1,1}.ADtoMV;
            end
            new_nexFileData.waves{unit_counter,1}.wireNumber = new_nexFileData.neurons{unit_counter,1}.wireNumber;
            new_nexFileData.waves{unit_counter,1}.unitNumber = new_nexFileData.neurons{unit_counter,1}.unitNumber;
            
            new_nexFileData.waves{unit_counter,1}.MVOffset = 0;
            new_nexFileData.waves{unit_counter,1}.timestamps = new_nexFileData.neurons{unit_counter,1}.timestamps;
            new_nexFileData.waves{unit_counter,1}.waveforms = waveforms;
            unit_counter = unit_counter+1;
        end
        
        
    else   %if size(mean_waveforms,2)>0, no waveforms, don't assign new sortings
     
        if options.keep_unsorted_ch_assignments  %Use previous sortings
            if plx_file_flag
                curr_orig_units = unique(cluster_indexes);
                waveform_assignments{ch} = cluster_indexes;
                for u = curr_orig_units'
                    if u == 0
                        new_nexFileData.neurons{unit_counter,1}.name = [channel_names(ch,:), 'U'];
                    else
                        new_nexFileData.neurons{unit_counter,1}.name = [channel_names(ch,:), char(96+u)];
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
            else
                curr_orig_units = find(orig_nexFileData_chans==ch);
                for u = curr_orig_units'
                    new_nexFileData.neurons{unit_counter,1} = nexFileData.neurons{u,1};
                    new_nexFileData.waves{unit_counter,1} = nexFileData.waves{u,1};
                    unit_counter = unit_counter+1;
                end
            end
        else  %All waveforms assigned to unsorted
            new_nexFileData.neurons{unit_counter,1}.name = ['sig', num2str(ch,'%03d'), 'U'];
            new_nexFileData.neurons{unit_counter,1}.varVersion = 101;
            new_nexFileData.neurons{unit_counter,1}.wireNumber = ch-1;
            new_nexFileData.neurons{unit_counter,1}.unitNumber = 0;
            new_nexFileData.neurons{unit_counter,1}.xPos = 0;
            new_nexFileData.neurons{unit_counter,1}.yPos = 0;
            new_nexFileData.neurons{unit_counter,1}.timestamps = timestamps;
            
            new_nexFileData.waves{unit_counter,1}.name = [new_nexFileData.neurons{unit_counter,1}.name, '_wf'];
            new_nexFileData.waves{unit_counter,1}.varVersion = new_nexFileData.neurons{unit_counter,1}.varVersion;
            new_nexFileData.waves{unit_counter,1}.NPointsWave = size(waveforms,1);
           
            if plx_file_flag
                new_nexFileData.waves{unit_counter,1}.WFrequency = sample_freq;
                new_nexFileData.waves{unit_counter,1}.ADtoMV = SpikePeakV/((2^(SpikeADResBits-1))*ChGains(ch));
            else
                if ~isempty(ch_indexes)
                    new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{ch_indexes(1),1}.WFrequency;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{ch_indexes(1),1}.ADtoMV;
                    
                else %If there's no waveforms on the channel, just use 1st channel sampling frequency and gain
                    new_nexFileData.waves{unit_counter,1}.WFrequency = nexFileData.waves{1,1}.WFrequency;
                    new_nexFileData.waves{unit_counter,1}.ADtoMV = nexFileData.waves{1,1}.ADtoMV;
                end
            end
            new_nexFileData.waves{unit_counter,1}.wireNumber = new_nexFileData.neurons{unit_counter,1}.wireNumber;
            new_nexFileData.waves{unit_counter,1}.unitNumber = new_nexFileData.neurons{unit_counter,1}.unitNumber;
            
            new_nexFileData.waves{unit_counter,1}.MVOffset = 0;
            new_nexFileData.waves{unit_counter,1}.timestamps = new_nexFileData.neurons{unit_counter,1}.timestamps;
            new_nexFileData.waves{unit_counter,1}.waveforms = waveforms;
            unit_counter = unit_counter+1;
        end
        
    end
end


for k = 1:length(new_nexFileData.neurons)
    new_nexFileData.neurons{k}.name((end+1):64) = 0;  %Add null character to name to make 64 character string padded by nulls, important for reading back into plexon offline sorter
    new_nexFileData.waves{k}.name((end+1):64) = 0;
end
writeNexFile(new_nexFileData, output_data_file);

toc
end


function SNR = signal2noise_PtoP2(waves)
%Peak to Peak signal-to-noise divided by 2
mean_signal = mean(waves,2);
signal = max(mean_signal) - min(mean_signal);
remain_noise = waves-repmat(mean_signal, [1, size(waves,2)]);
noise = std(remain_noise(:));
SNR = (signal./noise)/2;
end
