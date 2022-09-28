%Example script to call WaveLimit
%Adam Rouse, 5/30/19

% data_path = 'C:\DataFiles\';
% input_data_file = [data_path 'data_processed\monk_p\20160504_COT_precision\P_PlexonData_Centers\P_20170630_GHIJKLA+A_centers.nex'];
% output_data_file = [data_path 'data_processed\monk_p\20160504_COT_precision\P_PlexonData_NewAutosort\P_20170630_GHIJKLA+A-Autosort_test2.nex'];

% input_data_file = [data_path 'data_processed\Nomad test\Trellis_Test\BB_to_Spikes\Q_20190627_arrayH_test.nex'];
% output_data_file = [data_path 'data_processed\Nomad test\Trellis_Test\BB_to_Spikes\Q_20190627_arrayH_test_auto.nex'];

input_data_file = 'R:\SOM RSCH\RouseLab\DataFiles\Processed_Data\COTPerturb20210713\monk_A\BBtoSpikes\A_COTPerturb_Ped12_20210802_out.nex';
output_data_file = 'R:\SOM RSCH\RouseLab\DataFiles\Processed_Data\COTPerturb20210713\monk_A\Autosort\A_COTPerturb_Ped12_20210802_out_auto.nex';

added_tools_path = 'C:\Users\arouse\Box\added_matlab_tools\';
%Add Plexon offline SDK to path 
addpath(genpath([added_tools_path 'Matlab Offline Files SDK\'])) 
%Add Nex file reading/writing toolbox
addpath(genpath([added_tools_path 'HowToReadAndWriteNexAndNex5FilesInMatlab\'])) 

WaveLimit_path = '.\';
addpath([WaveLimit_path 'WaveLimit\']) 


options = default_options;
options.SNR_minimum = 0;  %Set to 0 for no SNR minimum
options.include_multiunits = true;  %Save all units even if there is no clear single unit based on ISI violations

ch_to_sort = [];  %use [] to sort all units
options.use_existing_clusters = false;  %True: Sort with centers defined by  input file file, False: use WaveLimit custer identification
options.keep_unsorted_ch_assignments = true; %True: Keep sorting assignments in output file of channels not sorted with WaveLimit call, False: make all waveforms unsorted on channels not sorted by WaveLimit 
options.keep_zero_clusters = true; %Keep units in file even if they have zero spikes because they didn't mean sort criteria, Mainly for use_existing_clusters = true
options.remove_line_noise_units = true;  %Remove units with too much 60Hz in their ISI pattern
options.remove_long_interval_units = true;  
WaveLimit(input_data_file, output_data_file, options,ch_to_sort)



