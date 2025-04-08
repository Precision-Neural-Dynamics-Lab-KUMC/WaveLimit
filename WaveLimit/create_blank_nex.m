function nexFileData = create_blank_nex()
% create_blank_nex for WaveLimit v1.3
% Author: Adam Rouse, Date: 4/23/20
%
% This function initializes a blank Nex file structure for use with WaveLimit.
% It creates the necessary fields for neurons and waveforms, with default
% values that can be modified later as needed.

% Nex file metadata
nexFileData.version = 104;             % Nex file version
nexFileData.comment = char(zeros(1,0)); % Empty comment field
nexFileData.freq = 40000;              % Sampling frequency in Hz
nexFileData.tbeg = 0;                  % Beginning time of recording
nexFileData.tend = 0;                  % End time of recording

% Initialize empty events and markers
nexFileData.events = {};               % Empty cell array for events
nexFileData.markers = {};              % Empty cell array for markers

% Initialize a default neuron structure
nexFileData.neurons{1}.name = 'sig001U'; % Name of the neuron
nexFileData.neurons{1}.varVersion = 101; % Neuron variable version
nexFileData.neurons{1}.wireNumber = 0;   % Wire number (set to 0 by default)
nexFileData.neurons{1}.unitNumber = 0;   % Unit number (set to 0 by default)
nexFileData.neurons{1}.xPos = 0;         % X position (default: 0)
nexFileData.neurons{1}.yPos = 0;         % Y position (default: 0)
nexFileData.neurons{1}.timestamps = [];  % Empty array for spike timestamps

% Initialize a corresponding waveform structure
nexFileData.waves{1}.name = [nexFileData.neurons{1}.name '_wf']; % Name of the waveform
nexFileData.waves{1}.varVersion = 101;    % Waveform variable version
nexFileData.waves{1}.NPointsWave = 32;    % Number of points per waveform
nexFileData.waves{1}.WFrequency = 40000;  % Waveform sampling frequency in Hz
nexFileData.waves{1}.wireNumber = nexFileData.neurons{1}.wireNumber; % Corresponding wire number
nexFileData.waves{1}.unitNumber = nexFileData.neurons{1}.wireNumber; % Corresponding unit number
nexFileData.waves{1}.ADtoMV = 1e-4;       % Conversion factor from A/D units to millivolts
nexFileData.waves{1}.MVOffset = 0;        % Millivolt offset
nexFileData.waves{1}.timestamps = [];     % Empty array for waveform timestamps
nexFileData.waves{1}.waveforms = zeros(nexFileData.waves{1}.NPointsWave, 0); % Initialize empty waveform matrix

end
