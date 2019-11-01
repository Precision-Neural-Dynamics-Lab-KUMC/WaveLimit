# WaveLimit
 Spike Sorting Software
WaveLimit is automatic spike-sorting that runs in Matlab
Created by Adam Rouse, v1.1, 11/1/19
It uses an input *.plx, *.pl2, or *.nex file and returns a sorted *.nex file 
An example call script is provided:  WaveLimit_example_call_script.m
options can be specified in the options structure after it is created with default_options.m

To use WaveLimit, you must have the following package:
HowToReadAndWriteNexAndNex5FilesInMatlab from NeuroExplorer
Available here:  
Code to Read and Write NeuroExplorer Data Files, Matlab code to read and write .nex and .nex5 files
https://www.neuroexplorer.com/downloadspage/

*NOTE, the writeNexFile function does not create files that read into Plexon Offline Sorter the same as the native export/import .NEX in Plexon. 
A modified writeNexFile.m is included with WaveLimit that is recommended for use with Plexon Offline Sorter.  
Simply replace writeNexFile.m in the HowToReadAndWriteNexAndNex5FilesInMatlab directory. 

To run, WaveLimit and HowToReadAndWriteNexAndNex5FilesInMatlab must be on your Matlab path.
You can either add it permanantly or include it using the example calls in WaveLimit_example_call_script.m

For reading native Plexon files you must have the OmniPlex and MAP Offline SDK Bundle
https://plexon.com/software-downloads/#software-downloads-SDKs
