function sorting_options = default_options()
%Create default sorting options for WaveLimit, v2.0
%Adam Rouse, 8/5/2025

%Setting needed based on input file
sorting_options.min_ISI = 0.675;  %Minimum time that can occur between spikes, needs to be known from hardware and spike thresholding software 

% Input/output file handling
sorting_options.use_existing_clusters = false;  %True: Sort with centers defined by  input file file, False: use WaveLimit custer identification
sorting_options.keep_unsorted_ch_assignments = true; %True: Keep sorting assignments in output file of channels not sorted with WaveLimit call, False: make all waveforms in other channels unsorted
sorting_options.keep_zero_clusters = false;  %Keep units in file even if they have zero spikes because didn't mean sort criteria, Mainly for use_existing_clusters = true so user can see which submitted clusters were then removed
sorting_options.scale_factor = 1; %How much to favor including units versus excluding units by scaling the chi^2 fit for fraction of spikes observed: 1 exactly maximize unit information phi value as in paper, >1 more waveforms included, <1 less waveforms
sorting_options.keep_analog = true;   %Keep analog channel data in .nex output file

% Default sorting options for WaveLimit spike sorting
sorting_options.max_ISI = 1.5;  %Maximum Interspike Interval time to be considered a violation in ms
sorting_options.include_multiunits = true;  %Keep units even if less than 1/2 waveforms appear to come from a single neuron
sorting_options.SNR_minimum = 0; %Keep units only if their SNR is above a certain value



% Long interval interspike interval unit removal, useful for trial-related noise like a spike with every reward or electrical stimulation artifact removal
sorting_options.remove_long_interval_units = true;  %Test if too many long interspike intervals occur
% When compared to a Poisson process, are less than sorting_options.not_long_interval_fraction (30%) of all interspike interval times less than sorting_options.long_interval_time (500 ms) inter-spike intervals?  
% if so then the spikes are occurring in a regular pattern with long time
% between them, this is likely due to regular stimulation or a trial
% artefact like a reward system
sorting_options.long_interval_time = 500;  %Time in ms to use for defining long interspike interval time
sorting_options.not_long_interval_fraction = 0.3; %Fraction of spikes that must be shorter than the long interspike interval time compared to a Poisson spiking unit, i.e. a Poisson process would have a value of 1
 


% Line noise removal settings, remove unit if too many spikes occur in a regular pattern relative to 60 Hz
sorting_options.remove_line_noise_units = true;  %Remove units that have interspike intervals that are periodic with power line noise
sorting_options.line_noise_freq = 60;  %Line noise frequency to test, 60Hz in USA 
% Test value for line noise frequency (60Hz), Rayleigh z test statistic to see if spikes are occurring in a regular pattern (non-uniform) within the 60Hz cycle, 
% high number of test statistic indicates alignment of spikes with 60Hz, thus setting this cutoff smaller is stricter, higher is more tolerant units with some line noise 
sorting_options.line_noise_cutoff = 500;   

% Random seed settings
sorting_options.repeatable_random_seed = 33;  %if a value: repeatable-random numbers with same seed, so sorting the same data will return the same result if 0: each run is truly random

% Cluster size settings
sorting_options.min_num_waveforms = 1000;  %Minimum number of waveforms needed to sort a channel
sorting_options.min_num_in_cluster = 500;  %Minimum number of waveforms needed in a cluster to save the cluster

% Sorting iteration settings
sorting_options.update_percentage = 10;  %Update factor for how far to move characteristic waveform from ISI violation causing waveforms 
sorting_options.align_waveforms = true;  %temporally align waveforms around peak or trough of waveform if true 
sorting_options.max_units_per_ch = 20;    %Maximum units that can be identified and sorted on a channel
sorting_options.max_iterations = 30; %Maximum iterations per channel of updating the characteristic waveform for each unit
sorting_options.max_multi_iterations = 3; %Maximum iterations per unit for adjusting characteristic waveform if it's always multiunit in each iteration


% GPU usage settings
if gpuDeviceCount>0
    sorting_options.use_gpu = true;
    try 
        % Check if the GPU is functioning properly
        tmp = gpuArray.zeros(1);
    catch
        warning('There is a problem with the GPU device.  Using CPU instead')
        sorting_options.use_gpu = false; % Fallback to CPU if GPU test fails
    end
    gpu_info = gpuDevice;  % Retrieve GPU device info
    sorting_options.max_gpu_memory = gpu_info.AvailableMemory;  % Set maximum GPU memory for sorting
else
    sorting_options.use_gpu = false;  % Use CPU if no GPU is available
end







