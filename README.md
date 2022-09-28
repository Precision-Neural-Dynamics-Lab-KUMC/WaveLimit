# WaveLimit
 Spike Sorting Software
WaveLimit is automatic spike-sorting that runs in Matlab
Created by Adam Rouse, v1.4, 9/28/22
It uses an input *.plx, *.pl2, or *.nex file and returns a sorted *.nex file 
An example call script is provided:  WaveLimit_example_call_script.m
options can be specified in the options structure after it is created with default_options.m

To use WaveLimit, you must have the following package:
HowToReadAndWriteNexAndNex5FilesInMatlab (Version Oct. 2019) from NeuroExplorer
Available here:  
Code to Read and Write NeuroExplorer Data Files, Matlab code to read and write .nex and .nex5 files
https://www.neuroexplorer.com/downloadspage/

*NOTE, the writeNexFile function does not create files that read into Plexon Offline Sorter the same as the native export/import .NEX/.NEX5 in Plexon. 
Modified writeNexFile.m and writeNex5File.m are included with WaveLimit that are recommended for use with Plexon Offline Sorter.  
Simply replace writeNexFile.m and writeNex5File.m in the HowToReadAndWriteNexAndNex5FilesInMatlab directory. 

To run, WaveLimit and HowToReadAndWriteNexAndNex5FilesInMatlab must be on your Matlab path.
You can either add it permanantly or include it using the example calls in WaveLimit_example_call_script.m

For reading native Plexon files you must have the OmniPlex and MAP Offline SDK Bundle
https://plexon.com/software-downloads/#software-downloads-SDKs


Current Release 1.4.0
<a href="https://zenodo.org/badge/latestdoi/219063673"><img src="https://zenodo.org/badge/219063673.svg" alt="DOI"></a>
