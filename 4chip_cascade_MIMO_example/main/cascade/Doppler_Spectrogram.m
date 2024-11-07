% Doppler_Spectrogram.m

% Clear workspace and command window
clearvars;
close all;
%clc;

% Initialize switch
NORMALIZE_SPECTOGRAM = 1;

% Specify the path to the saved .mat file
testName = 'drone_steady_1';
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
    currentSigIntegrate = sig_integrate_all{frame}; % [256, 64]

    % Sum across range bins to get Doppler spectrum for the current frame
    % You can also use mean instead of sum if preferred
    doppler_spectrogram(frame, :) = sum(currentSigIntegrate, 1); % [1, 64]
end

% Normalize the Doppler Spectrogram (optional)
% This can help in better visualization
if NORMALIZE_SPECTOGRAM == 1
    doppler_spectrogram = doppler_spectrogram - min(doppler_spectrogram(:));
    doppler_spectrogram = doppler_spectrogram / max(doppler_spectrogram(:));
end


% Create the Doppler Spectrogram plot
figure('Name', 'Doppler Spectrogram', 'NumberTitle', 'off');

% Use imagesc to create a heatmap
imagesc(doppler_spectrogram);
colormap('parula'); % Choose a colormap (e.g., jet, hot, parula)
colorbar; % Display color scale
axis xy; % Ensure the y-axis starts from the bottom

% Label the axes
xlabel('Doppler Bins');
ylabel('Frame ID');
title(['Doppler Spectrogram of ' testName]);

% Customize the axes ticks if needed
% Example: Label Doppler bins with actual Doppler velocities
% Assuming you have 'dopplerBinSize' parameter loaded
if isfield(loadedData, 'dopplerBinSize') && loadedData.dopplerBinSize ~= 0
    dopplerBinSize = loadedData.dopplerBinSize;
    dopplerRange = (0:numDopplerBins-1) * dopplerBinSize;
    set(gca, 'XTick', 1:5:numDopplerBins, 'XTickLabel', round(dopplerRange(1:5:end),2));
    xlabel('Doppler Velocity (m/s)');
end

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
sprintf(testName);
