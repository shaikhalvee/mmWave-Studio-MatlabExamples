% Doppler_Spectrogram_Drone.m

% Clear workspace and command window
clear;
clc;

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
loadVariable = load(dopplerSpectrogramFile, 'doppler_spectrogram');

% Validate data existence
if ~isfield(loadedData, 'sig_integrate_all')
    error('The variable ''sig_integrate_all'' does not exist in the loaded file.');
end

if ~isfield(loadVariable, 'doppler_spectrogram')
    error('the variable ''doppler_spectrogram'' does not exist.')
end

sig_integrate_all = loadedData.sig_integrate_all;
doppler_spectrogram = loadVariable.doppler_spectrogram;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Wrapping Interval

% now static, but implement later to bring from previous dataset
lambda = 3.888536e-03;

% Doppler Velocity Calculation
%dopplerVelocities = (dopplerIndices - zeroDopplerBin) * dopplerBinSize;  % Velocity values in m/s
%frequencies = dopplerVelocities * 2 / lambda;

no_of_blades = 2;

% Rotator frequency interval range
rotator_freq_range =  (2500:7000) / 60; % range from 2500 rpm to 4500 rpm

doppler_bin_range = round((lambda * no_of_blades / (2 * dopplerBinSize)) * rotator_freq_range);
chopping_intervals = unique(doppler_bin_range);

% Create a container to store key-value pairs
% wrap_bin_wrap_value_keypair = containers.Map('KeyType', 'double', 'ValueType', 'double');

wrap_bin_wrap_value = zeros(length(chopping_intervals), 2); % 2 columns: one for wrapping intervals, one for wrapping values

max_wrapped_value = -Inf; % Initialize max folding value for this spectrum
best_wrap_val = 0; % Initialize best wrapping value

% let's consider the 250th frame only
doppler_spectrum = doppler_spectrogram(250, :);

wrapped_matrix_list = {};

index = 1; % Index for filling rows of wrap_bin_wrap_value

for wrap_interval = chopping_intervals
    % Calculate M = floor(length of Doppler spectrum / wrap_val)
    M = floor(length(doppler_spectrum) / wrap_interval);

    % Check if wrapping is feasible
    if M < 2
        continue;
    end

    % Reshape Doppler spectrum into [M, wrap_val] matrix
    wrapped_matrix = reshape(doppler_spectrum(1:M*wrap_interval), [M, wrap_interval]);

    % Compute column-wise average
    column_avg = mean(wrapped_matrix, 1);
        
    % Compute folding value: maximum of column averages
    wrapped_value = max(column_avg);

    % Store wrap_interval and wrap_value in wrap_bin_wrap_value matrix
    wrap_bin_wrap_value(index, :) = [wrap_interval, wrapped_value];

    wrapped_matrix_list{index} = wrapped_matrix;

    index = index + 1;
end

%{
% Preallocate results
wrapping_results = zeros(size(doppler_spectrogram, 1), 1); % Stores P_i for each spectrum
optimal_wrapping = zeros(size(doppler_spectrogram, 1), 1); % Stores the optimal wrapping interval for each spectrum

% Loop through each Doppler spectrum
for frame_idx = 1:size(doppler_spectrogram, 1)
    doppler_spectrum = doppler_spectrogram(frame_idx, :);
    max_wrapped_value = -Inf; % Initialize max folding value for this spectrum
    best_wrap_val = 0; % Initialize best wrapping value

    % Loop through wrapping intervals
    for wrap_val = wrapping_interval
        % Calculate M = floor(length of Doppler spectrum / wrap_val)
        M = floor(length(doppler_spectrum) / wrap_val);
        
        % Check if wrapping is feasible
        if M < 2
            continue;
        end
        
        % Reshape Doppler spectrum into [M, wrap_val] matrix
        wrapped_matrix = reshape(doppler_spectrum(1:M*wrap_val), [M, wrap_val]);
        
        % Compute column-wise average
        column_avg = mean(wrapped_matrix, 1);
        
        % Compute folding value: maximum of column averages
        wrapped_value = max(column_avg);
        
        % Update maximum folding value and wrapping interval
        if wrapped_value > max_wrapped_value
            max_wrapped_value = wrapped_value;
            best_wrap_val = wrap_val;
        end
    end
    
    % Store results
    wrapping_results(frame_idx) = max_wrapped_value; % P_i for this spectrum
    optimal_wrapping(frame_idx) = best_wrap_val; % Optimal wrapping interval
end

% Display results
disp('Folding results (P_i) for each spectrum:');
disp(wrapping_results);

disp('Optimal wrapping intervals for each spectrum:');
disp(optimal_wrapping);
%}

