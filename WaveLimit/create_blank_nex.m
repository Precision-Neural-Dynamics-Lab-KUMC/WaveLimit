function nexFileData = create_blank_nex()
%create_blank_nex for WaveLimit v1.3
% Adam Rouse, 4/23/20

nexFileData.version = 104;
nexFileData.comment = char(zeros(1,0));
nexFileData.freq = 40000;
nexFileData.tbeg = 0;
nexFileData.tend = 0;  
nexFileData.events = {};
nexFileData.markers = {};
nexFileData.neurons{1}.name = 'sig001U';
nexFileData.neurons{1}.varVersion = 101;
nexFileData.neurons{1}.wireNumber = 0;
nexFileData.neurons{1}.unitNumber = 0;
nexFileData.neurons{1}.xPos = 0;
nexFileData.neurons{1}.yPos = 0;
nexFileData.neurons{1}.timestamps = [];

nexFileData.waves{1}.name =  [nexFileData.neurons{1}.name '_wf'];
nexFileData.waves{1}.varVersion = 101;
nexFileData.waves{1}.NPointsWave = 32;
nexFileData.waves{1}.WFrequency = 40000;
nexFileData.waves{1}.wireNumber = nexFileData.neurons{1}.wireNumber;
nexFileData.waves{1}.unitNumber = nexFileData.neurons{1}.wireNumber;
nexFileData.waves{1}.ADtoMV = 1e-4;
nexFileData.waves{1}.MVOffset = 0;
nexFileData.waves{1}.timestamps = [];
nexFileData.waves{1}.waveforms = zeros(nexFileData.waves{1}.NPointsWave,0);