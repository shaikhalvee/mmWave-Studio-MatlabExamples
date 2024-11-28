% Doppler_Spectrogram_Drone.m

% Clear workspace and command window
clear;
clc;

% Initialize switch
NORMALIZE_SPECTOGRAM = 0;

% Specify the path to the saved .mat file
pro_path = getenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO');
input_path = strcat(pro_path,'\main\cascade\input\');
testList = strcat(input_path,'testList.txt');
fid = fopen(testList, 'r');
line = fgetl(fid); % Read the first line
fclose(fid);

% Extract the last folder name using fileparts
[~, testName] = fileparts(fileparts(line));
outputFile = ['.\main\cascade\output\newOutput_', testName, '.mat'];
dopplerSpectrogramFile = ['.\main\cascade\output\doppler_spectrogram_', testName, '.mat'];

% Check if the file exists
if ~isfile(outputFile)
    error('The specified output file does not exist: %s', outputFile);
end

if ~isfile(dopplerSpectrogramFile)
    error('%s does not exist', dopplerSpectrogramFile);
end

% Load the required data
loadedData = load(outputFile, 'sig_integrate_all', 'rangeBinSize', 'dopplerBinSize', 'detection_results_all');
doppler_spectrogram = load(dopplerSpectrogramFile, 'doppler_spectrogram');

% Validate data existence
if ~isfield(loadedData, 'sig_integrate_all')
    error('The variable ''sig_integrate_all'' does not exist in the loaded file.');
end

sig_integrate_all = loadedData.sig_integrate_all;

% Retrieve Doppler Bin Size
if isfield(loadedData, 'dopplerBinSize') && loadedData.dopplerBinSize ~= 0
    dopplerBinSize = loadedData.dopplerBinSize; % Velocity per Doppler bin (m/s)
else
    error('dopplerBinSize is not available or is zero in the loaded data.');
end

% Determine the number of frames and Doppler bins
numFrames = length(sig_integrate_all);
if numFrames == 0
    error('No frames found in ''sig_integrate_all''.');
end

% Assuming all frames have the same dimensions
[numRangeBins, numDopplerBins] = size(sig_integrate_all{1});



