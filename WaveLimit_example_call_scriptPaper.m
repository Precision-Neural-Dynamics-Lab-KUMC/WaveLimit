% Example script to call WaveLimit for spike sorting
% Author: Adam Rouse, Date: 5/30/19

% Define paths to the data files
data_path = 'C:\DataFiles\';
% input_data_file = [data_path 'data_processed\monk_p\20160504_COT_precision\P_PlexonData_Centers\P_20170630_GHIJKLA+A_centers.nex'];
% output_data_file = [data_path 'data_processed\monk_p\20160504_COT_precision\P_PlexonData_NewAutosort\P_20170630_GHIJKLA+A-Autosort_test2.nex'];

% input_data_file = [data_path 'data_processed\Nomad test\Trellis_Test\BB_to_Spikes\Q_20190627_arrayH_test.nex'];
% output_data_file = [data_path 'data_processed\Nomad test\Trellis_Test\BB_to_Spikes\Q_20190627_arrayH_test_auto.nex'];

input_data_file = 'R:\SOM RSCH\RouseLab\DataFiles\ArchivedProjects\SchieberLab\data_processed\monk_p\20160504_COT_precision\BB_to_Spikes\P_20170705_GHIJxxxx_BBout.nex';
output_data_file = 'R:\SOM RSCH\RouseLab\DataFiles\ArchivedProjects\SchieberLab\data_processed\monk_p\SpikeSortingPaper\P_20170705_GHIJxxxx_BB1-64_auto.nex';
% input_data_file = 'R:\SOM RSCH\RouseLab\DataFiles\Processed_Data\COTPerturb20210713\monk_A\BBtoSpikes\A_COTPerturb_Ped12_20210802_out.nex';
% output_data_file = 'R:\SOM RSCH\RouseLab\DataFiles\Processed_Data\COTPerturb20210713\monk_A\Autosort\A_COTPerturb_Ped12_20210802_out_auto.nex';

% Define paths to additional MATLAB tools
added_tools_path = 'R:\SOM RSCH\RouseLab\LabWork\DataProcessing\DataProcessingToolboxes\';
%Add Plexon offline SDK to path 
addpath(genpath([added_tools_path 'Matlab Offline Files SDK\'])) 
%Add Nex file reading/writing toolbox
addpath(genpath([added_tools_path 'HowToReadAndWriteNexAndNex5FilesInMatlab\'])) 

WaveLimit_path = '.\';
addpath([WaveLimit_path 'WaveLimit\']) 


% Add WaveLimit toolbox to MATLAB path
WaveLimit_path = '.\';
addpath([WaveLimit_path 'WaveLimit\']);

% Set WaveLimit options using default configuration
options = defaultWaveLimitOptions;

% Modify specific options for this analysis
options.SNR_minimum = 0;  % Set to 0 for no minimum SNR threshold
options.include_multiunits = true;  % Save all units, including multi-units with ISI violations

% Channels to sort: use an empty array to sort all available channels
ch_to_sort =  14; %[];

% Sorting method options
options.use_existing_clusters = false;  % False: Use WaveLimit's cluster identification, True: use input file clusters
options.keep_unsorted_ch_assignments = true;  % Keep channel assignments for unsorted channels in the output file
options.keep_zero_clusters = true;  % Keep clusters with zero spikes (useful if use_existing_clusters is true)
options.remove_line_noise_units = true;  % Remove units with excessive 60Hz noise in their ISI pattern
options.remove_long_interval_units = true;  % Remove units with prolonged ISI intervals

% Call WaveLimit function to perform the sorting
WaveLimit(input_data_file, output_data_file, options, ch_to_sort);
