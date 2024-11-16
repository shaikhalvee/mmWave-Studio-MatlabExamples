% Doppler_Spectrogram_Updated.m

% Clear workspace and command window
clearvars
close all
% clc;

% Initialize switch
NORMALIZE_SPECTOGRAM = 1;

% Specify the path to the saved .mat file
pro_path = getenv('CASCADE_SIGNAL_PROCESSING_CHAIN_MIMO');
input_path = strcat(pro_path,'\main\cascade\input\');
testList = strcat(input_path,'testList.txt');
fid = fopen(testList, 'r');
line = fgetl(fid); % Read the first line
fclose(fid);
%[inputFolder, ~] = fileparts(line); % Extract the last folder name using fileparts
[~, testName] = fileparts(fileparts(line));
%testName = 'drone_data';
outputFile = ['.\main\cascade\output\newOutput_', testName, '.mat'];

% Check if the file exists
if ~isfile(outputFile)
    error('The specified output file does not exist: %s', outputFile);
end

% Load the required data
loadedData = load(outputFile, 'sig_integrate_all', 'rangeBinSize', 'dopplerBinSize');

% Check if 'sig_integrate_all' exists in the loaded data
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

% Assuming all frames have the same Doppler bin size
[numRangeBins, numDopplerBins] = size(sig_integrate_all{1});

% Initialize a matrix to hold Doppler spectra for all frames
doppler_spectrogram = zeros(numFrames, numDopplerBins);

% Aggregate Doppler spectrum for each frame
for frame = 1:numFrames
    currentSigIntegrate = sig_integrate_all{frame}; % [RangeBins, DopplerBins]

    % Sum across range bins to get Doppler spectrum for the current frame
    doppler_spectrogram(frame, :) = sum(currentSigIntegrate, 1); % [1, DopplerBins]
end

% Normalize the Doppler Spectrogram for better visualization (optional)
% This step enhances contrast by scaling data between 0 and 1
if NORMALIZE_SPECTOGRAM == 1
    doppler_spectrogram = doppler_spectrogram - min(doppler_spectrogram(:));
    doppler_spectrogram = doppler_spectrogram / max(doppler_spectrogram(:));
end

% Create the Doppler Spectrogram plot
figure('Name', 'Doppler Spectrogram (Logarithmic Scale)', 'Position', [400 400 900 600]);


% Use imagesc to create a heatmap
imagesc(doppler_spectrogram);
colormap('parula'); % Choose a colormap (e.g., jet, hot, parula)
colorbar; % Display color scale
axis xy; % Ensure the y-axis starts from the bottom

% Label the axes
xlabel('Doppler Velocity (m/s)');
ylabel('Frame ID');

% Calculate Doppler Velocity Range for X-axis Labels
dopplerIndices = 1:numDopplerBins;
zeroDopplerBin = ceil(numDopplerBins / 2) + 1; % For FFT size 64, zero Doppler is bin 33

% Doppler Velocity Calculation
dopplerVelocities = (dopplerIndices - zeroDopplerBin) * dopplerBinSize;

% Set X-axis Ticks and Labels
numTicks = 8; % Number of ticks on the x-axis
tickInterval = ceil(numDopplerBins / numTicks);

% Generate Tick Positions
%tickPositions = 1:tickInterval:numDopplerBins;
tickPositions = [1:tickInterval:numDopplerBins, numDopplerBins]; % Ensure the last tick is included

% Generate Corresponding Doppler Velocity Labels
%tickLabels = dopplerVelocities(tickPositions);
tickLabels = round(dopplerVelocities(tickPositions), 2); % Round to 2 decimal places for readability

% Apply X-axis Tick Labels
set(gca, 'XTick', tickPositions, 'XTickLabel', tickLabels);

% Customize the y-axis if needed
% Example: Label frames with actual timestamps or frame rates
% Assuming a frame rate 'frameRate' (frames per second)
% frameRate = 30; % Example value
% frameTimes = (1:numFrames) / frameRate;
% set(gca, 'YTick', 1:10:numFrames, 'YTickLabel', round(frameTimes(1:10:end),2));
% ylabel('Time (s)');

% Improve visualization aesthetics
set(gca, 'FontSize', 12);
testName = regexprep(testName, '_', '\\_');
title(sprintf('Doppler Spectrogram of %s', testName), 'FontSize', 14);

% Optional: Add grid lines for better readability
grid on;
