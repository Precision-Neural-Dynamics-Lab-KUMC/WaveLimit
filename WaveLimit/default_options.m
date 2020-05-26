function sorting_options = default_options()
%Create default sorting options for WaveLimit, v1.1
%Adam Rouse, 4/23/20

sorting_options.max_ISI = 1.5;  %Maximum Interspike Interval time to be considered a violation in ms
sorting_options.include_multiunits = true;  %Keep units even if less than 1/2 waveforms appear to come from a single neuron
sorting_options.SNR_minimum = 1;
sorting_options.remove_long_interval_units = true;  %Test if excessive long interspike intervals occur, useful for trial-related noise like a spike with every reward
sorting_options.long_interval_time = 500;  %Time in ms to use for defining long interspike interval time
sorting_options.not_long_interval_fraction = 0.3;  %It's more reliable to test for the number of spikes that are not long interspike intervals and if it is less than this fraction of the expected value, then throwout this unit
sorting_options.remove_line_noise_units = true;  %Remove units that have interspike intervals that are periodic with power line noise
sorting_options.line_noise_freq = 60;  %Line noise frequency to test, 60Hz in USA, 
sorting_options.line_noise_cutoff = 500;  % Test value for line noise frequency (60Hz), Rayleigh z test statistic,
sorting_options.min_num_waveforms = 1000;  %Minimum number of waveforms need
sorting_options.min_num_in_cluster = 500;  %Minimum number of waveforms needed in a cluster
sorting_options.throwout_nle_outliers = true; %Throw out nonlinear energy outliers if true
sorting_options.update_percentage = 10;  %Update factor for how far to move central waveform from ISI violation causing waveforms 
sorting_options.align_waveforms = true;  %temporally align waveforms around peak or trough of waveform if true 
sorting_options.max_units_per_ch = 20;    %Maximum units that can be identified and sorted on a channel
sorting_options.max_iterations = 30;
sorting_options.max_multi_iterations = 3;
sorting_options.use_existing_clusters = false;
sorting_options.keep_unsorted_ch_assignments = true;
sorting_options.scale_factor = 1; %How much to favor including units versus excluding units: 1 exactly maximize unit information, >1 more waveforms inclued, <1 less waveforms
sorting_options.repeatable_random_seed = 1; %If a number - the random number generator seed is reset to that value so it is repeatable, if 0-then each iteration is unique
sorting_options.use_existing_clusters = false;  %True: Sort with centers defined by  input file file, False: use WaveLimit custer identification
sorting_options.keep_unsorted_ch_assignments = true; %True: Keep sorting assignments in output file of channels not sorted with WaveLimit call, False: make all waveforms in 
sorting_options.keep_zero_clusters = false;  %Keep units in file even if they have zero spikes because didn't mean sort criteria, Mainly for use_existing_clusters = true
sorting_options.repeatable_random_seed = 33;  %if a value: repeatable-random numbers with same seed, if 0: each run is truly random
sorting_options.keep_analog = true;   %Keep analog channel data in .nex output file
if gpuDeviceCount>0
    sorting_options.use_gpu = true;
    try 
        tmp = gpuArray.zeros(1);
    catch
        warning('There is a problem with the GPU device.  Using CPU instead')
        sorting_options.use_gpu = false;
    end
    gpu_info = gpuDevice;
    sorting_options.max_gpu_memory = gpu_info.AvailableMemory;
else
    sorting_options.use_gpu = false;
end







