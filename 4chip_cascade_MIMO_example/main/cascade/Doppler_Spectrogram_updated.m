% Doppler_Spectrogram_Updated.m

% Clear workspace and command window
clearvars
close all
% clc;

% Initialize switch
%NORMALIZE_SPECTOGRAM = 0;
INCLUDE_RANGE_BINS = 1;
RANGE_DOPPLER_PRINT = 0;
SHOW_DETECTION_RESULT = 0;

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

% Define the range bins of interest (e.g., where the drone is located)
% Adjust 'rangeStart' and 'rangeEnd' based on where the drone appears in range bins
rangeStart = 1; % Example value (in range bin indices)
rangeEnd = 1;  % Example value (in range bin indices)

% Aggregate Doppler spectrum for each frame
for frame = 1:numFrames
    currentSigIntegrate = sig_integrate_all{frame}; % [RangeBins, DopplerBins]

    if INCLUDE_RANGE_BINS == 1
        % Sum across the range bins where the drone is expected
        % This focuses on the target area and improves SNR
        rangeBinsOfInterest = rangeStart:rangeEnd;
        % doppler_spectrogram(frame, :) = sum(currentSigIntegrate(rangeBinsOfInterest, :), 1); % [1, DopplerBins]
        doppler_spectrogram(frame, :) = currentSigIntegrate(rangeBinsOfInterest, :); % [1, DopplerBins]
    else
        % Sum across range bins to get Doppler spectrum for the current frame
        doppler_spectrogram(frame, :) = sum(currentSigIntegrate, 1); % [1, DopplerBins]
    end
end

% doppler_spectrogram = 10 * log10(doppler_spectrogram + 1);

% Normalize the Doppler Spectrogram for better visualization (optional)
% This step enhances contrast by scaling data between 0 and 1
%if NORMALIZE_SPECTOGRAM == 1
%    doppler_spectrogram = doppler_spectrogram - min(doppler_spectrogram(:));
%    doppler_spectrogram = doppler_spectrogram / max(doppler_spectrogram(:));
%end

% Create the Doppler Spectrogram plot
figure('Name', 'Doppler Spectrogram (Logarithmic Scale)', 'Position', [400 400 900 600], 'NumberTitle', 'off');


% Use imagesc to create a heatmap
imagesc(doppler_spectrogram);
colormap('parula'); % Choose a colormap (e.g., jet, hot, parula)
colorbar; % Display color scale
axis xy; % Ensure the y-axis starts from the bottom

% Label the axes
% xlabel('Doppler Velocity (m/s)');
xlabel('Doppler bins', 'FontWeight', 'bold');
ylabel('Frame ID', 'FontWeight', 'bold');

% Calculate Doppler Velocity Range for X-axis Labels
dopplerIndices = 1:numDopplerBins;
% Adjust if fftshift was applied
zeroDopplerBin = ceil(numDopplerBins / 2) + 1; % For FFT size 64, zero Doppler is bin 33 

% Doppler Velocity Calculation
% dopplerVelocities = (dopplerIndices - zeroDopplerBin) * dopplerBinSize;  % Velocity values in m/s
dopplerBinNumber = (dopplerIndices - zeroDopplerBin);

% Set X-axis Ticks and Labels
xNumTicks = 8; % Number of ticks on the x-axis
tickInterval = ceil(numDopplerBins / xNumTicks);

% Generate Tick Positions
%tickPositions = 1:tickInterval:numDopplerBins;
tickPositions = [1:tickInterval:numDopplerBins, numDopplerBins]; % Ensure the last tick is included

% Generate Corresponding Doppler Velocity Labels
%tickLabels = dopplerVelocities(tickPositions), 2); % Round to 2 decimal places for readability
tickLabels = round(dopplerBinNumber(tickPositions), 2);

% Apply X-axis Tick Labels
set(gca, 'XTick', tickPositions, 'XTickLabel', tickLabels);

% Customize the y-axis
% Example: Label frames with actual timestamps or frame rates
% Assuming a frame rate 'frameRate' (frames per second)
% frameRate = 30; % Example value
% frameTimes = (1:numFrames) / frameRate;
% set(gca, 'YTick', 1:10:numFrames, 'YTickLabel', round(frameTimes(1:10:end),2));
% ylabel('Time (s)');

% Set y-axis labels to represent frame IDs (650 to 1221)
%frameIDs = 650:1221;
%numYTicks = 10; % Number of ticks on the y-axis
%yTickInterval = ceil(length(frameIDs) / numYTicks);
%yTickPositions = 1:yTickInterval:length(frameIDs);
%yTickLabels = frameIDs(yTickPositions);

% Apply Y-axis Tick Labels
%set(gca, 'YTick', yTickPositions, 'YTickLabel', yTickLabels);

% Improve visualization aesthetics
set(gca, 'FontSize', 14, 'LineWidth', 1.5);
%testNewName = regexprep(testName, '_', '\\_');
%title(sprintf('Doppler Spectrogram of %s', testNewName), 'FontSize', 14);
title('Doppler Spectrogram', 'FontSize', 14, 'FontWeight','bold')

% Optional: Add grid lines for better readability
grid on;

% Save Doppler Spectrogram
outputVariablePath = ['.\main\cascade\output\doppler_spectrogram_', testName, '.mat'];
save(outputVariablePath, 'doppler_spectrogram')

% Save image
outputImageFilePath = ['.\main\cascade\output\images\' , testName, '.pdf'];
f = gca;
exportgraphics(gca, outputImageFilePath, "ContentType", 'vector');

%saveas(gcf, ['Doppler_Spectrogram_Drone_', testName, '.png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if RANGE_DOPPLER_PRINT == 1
    % Plot the Range Doppler image for the 100th frame
    figure('Name', 'Range Doppler Image - Frame 100', 'NumberTitle', 'off');

    % Use imagesc to plot the data for the 100th frame
    imagesc(sig_integrate_all{100});

    % Show the colorbar
    c = colorbar;

    % Set colormap
    colormap('parula'); % Optional, choose your preferred colormap
    c.Label.String = 'Relative Power(dB)';

    % Label the axes
    xlabel('Doppler Bins', 'FontSize', 14, 'FontWeight', 'bold'); % Larger font size for x-axis
    ylabel('Range Bins', 'FontSize', 14, 'FontWeight', 'bold'); % Larger font size for y-axis

    % Set ticks to be more visible
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'XTick', tickPositions, 'XTickLabel', tickLabels); % Increase tick font size and line width

    % Saving single range doppler
    title('Range Doppler Image - Frame 100', 'FontSize', 16, 'FontWeight', 'bold'); % Title with larger font size
    title('Range Doppler', 'FontSize', 14, 'FontWeight', 'bold');
    outputImageFilePath = ['.\main\cascade\output\images\' , 'Range Doppler of frame 100', '.pdf'];
    f = gca;
    exportgraphics(gca, outputImageFilePath, "ContentType", 'vector');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if SHOW_DETECTION_RESULT == 1
    if isfield(loadedData, 'detection_results_all')
    hold on;
    for frame = 1:numFrames
        detection_results = loadedData.detection_results_all{frame};
        if ~isempty(detection_results)
            % Extract Doppler indices and convert to velocities
            dopplerInds = arrayfun(@(x) x.dopplerInd, detection_results);
            dopplerVels = (dopplerInds - zeroDopplerBin) * dopplerBinSize;
            % Plot detections
            plot(timeAxis(frame) * ones(size(dopplerVels)), dopplerVels, 'k.', 'MarkerSize', 10);
        end
    end
    hold off;
    end
end
% Optional: Save the figure
% saveas(gcf, ['Doppler_Spectrogram_Drone_', testName, '.png']);
