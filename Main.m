%% Importing the .wav files (3)
clear;

% The folder where the .wav files are located (relative path)
% 56 is the data set allocated to me
dataFolder = 'data/56';

% Join them togheter to make the full folder path
folderPath = fullfile(dataFolder);

% Get a list of .wav files in the folder
fileList = dir(fullfile(folderPath, '*.wav'));

% Initialize array to store audio data
audioDataArray = cell(numel(fileList), 1);

% Loop through each .wav file and import it
for i = 1:numel(fileList)
    filePath = fullfile(folderPath, fileList(i).name);
    
    % Read the .wav file using the audioread function
    [audioData, ~] = audioread(filePath);
    
    % Save the audio data to the cell array
    audioDataArray{i} = audioData;
end

%% Ploting one example of the signals imported (4)

% Amount of times each digit was recorded and the order it was said in
groupSize = 50;
order = 6;

% Number of digits
numGroups = 10;

% Sample rate (in Hz)
sampleRate = 48000;

% Determine the number of rows and columns for the subplot grid
numRows = ceil(sqrt(numGroups));
numCols = ceil(numGroups / numRows);

% Figure where we are ploting the signals
figure(1);

% Iterate through the array of signals
for i = 1:numGroups
    % Calculate the duration of the signal (in seconds)
    duration = size(audioDataArray{i*order}, 1) / sampleRate;
    time = linspace(0, duration, size(audioDataArray{i*order}, 1));
    
    % Create a subplot for each digit
    subplot(numRows, numCols, i);
    
    % Plot the audio data
    plot(time, audioDataArray{i*order});
    xlabel('time (s)');
    ylabel('Amp.');
    title(sprintf('|Digit %d - Order %d|', i-1, order));
    grid on;
end

% Figure title/Resize the figure for better visibility
sgtitle('Audio signals');
set(gcf, 'Position', [100, 100, 1000, 800]);

%% Identifying temporal features to differentiate the digits (4.2)

% Initialize arrays to store the temporal features
averageMaxAmplitude = zeros(numGroups, 1);
averageEnergy = zeros(numGroups, 1);
averageDuration = zeros(numGroups, 1);
% Ploting maximum amplitude
% Ploting settings
numRows = 5;
numCols = 2;
numPlots = numRows*numCols;

% Iterate over the groups
for i = 1:numGroups
    startIndex = (i - 1) * groupSize + 1;
    endIndex = min(startIndex + groupSize - 1, numel(audioDataArray));
    groupData = audioDataArray(startIndex:endIndex);
    
    % Calculate the maximum amplitude, energy and duration for each dataset
    maxValues = zeros(size(groupData));
    energyValues = zeros(size(groupData));
    durationValues = zeros(size(groupData));
    
    for j = 1:numel(groupData)
        maxValues(j) = max(groupData{j});
        energyValues(j) = sum(groupData{j}.^2);        
        % Remove begining and end silence
        index = find(groupData{j} > 0); % Adjust the threshold if needed
        durationValues(j) = index(end) - index(1) + 1;
    end
    % Saving information 
    averageMaxAmplitude(i)= mean(maxValues);
    averageEnergy(i) = mean(energyValues);
    averageDuration(i) = mean(durationValues);
    
    % Calculate the subplot position for the current group
    subplotPosition = mod(i - 1, numPlots) + 1;
    
    figure(2)
    % Create a new subplot for each group
    subplot(numRows, numCols, subplotPosition);
    histogram(maxValues);
    
    % Disable scientific notation on the x-axis
    ax = gca;
    ax.XAxis.Exponent = 0;

    % Add labels and title
    xlabel('Value');
    ylabel('Frequency');
    title(sprintf('|Digit %d|Mean: %.2e|Std: %.2e|Max. value: %.2e|', i-1, median(maxValues), std(maxValues), max(maxValues)));
    grid on;
    
    figure(3)
    % Create a new subplot for each group
    subplot(numRows, numCols, subplotPosition);
    
    histogram(energyValues);
    
    % Disable scientific notation on the x-axis
    ax = gca;
    ax.XAxis.Exponent = 0;
    xlabel('Data index');
    ylabel('Energy');
    title(sprintf('|Digit %d| Mean: %.2e | Std: %.2e | Max. value: %.2e |', i-1, median(energyValues), std(energyValues), max(energyValues)));
    grid on;
    
    figure(4)
    % Create a new subplot for each group
    subplot(numRows, numCols, subplotPosition);
    
    histogram(durationValues);
    
    % Disable scientific notation on the x-axis
    ax = gca;
    ax.XAxis.Exponent = 0;
    xlabel('Data index');
    ylabel('Duration (in sec.)');
    title(sprintf('|Digit %d| Mean: %.3f | Std: %.3f | Max. value: %.3f |', i-1, median(durationValues)/sampleRate, std(durationValues)/sampleRate, max(durationValues)/sampleRate));
    grid on;
end

% Figure title and resize the figure for better visibility
figure(2)
sgtitle('Maximum Amplitude by digit');
set(gcf, 'Position', [100, 100, 1200, 800]);

figure(3);
sgtitle('Energy by digit');
set(gcf, 'Position', [100, 100, 1200, 800]);

figure(4);
sgtitle('Duration by digit');
set(gcf, 'Position', [100, 100, 1200, 800]);

%% Formulating a set of if-then-else rules to differentiate the digits and ploting the decision features (4.3)
% Normalize the features
normalizedAmplitude = (averageMaxAmplitude - min(averageMaxAmplitude(:))) / (max(averageMaxAmplitude(:)) - min(averageMaxAmplitude(:)));
normalizedEnergy = (averageEnergy - min(averageEnergy(:))) / (max(averageEnergy(:)) - min(averageEnergy(:)));

% Plot
digitLabels = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
figure(5);
for i = 1:numel(digitLabels)
    caption = digitLabels{i};
    scatter(normalizedAmplitude(i), normalizedEnergy(i), 'filled');
    x = normalizedAmplitude(i);
    y = normalizedEnergy(i);
    text(x, y, sprintf('%s (%.2f, %.2f)', caption, x, y), 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
    if i==1
        hold on;
    end
end

%{
% This code is to draw lines on the scatter to make it easier to get the
% regions
hold on;
% Define the equations of the lines
line1 = @(x) -2*x + 1;  
line2 = @(x) -2*x + 2;  

% Determine the desired range for the lines
xMin = 0;
xMax = 1;  % Adjust the maximum x-value as needed

% Generate x-values for the lines
xLine = linspace(xMin, xMax, 100);

% Plot the lines
plot(xLine, line1(xLine), 'r-', 'LineWidth', 2);
plot(xLine, line2(xLine), 'g-', 'LineWidth', 2);
%}

% Remove the legend
legend('off');

% Label the axis
xlabel('Normalized Average Maximum Amplitude');
ylabel('Normalized Average Energy');

% Figure title and resize the figure for better visibility
sgtitle('Temporal features for the decision');
set(gcf, 'Position', [100, 100, 1200, 800]);


dividedDigits = zeros(1, numel(numGroups));  % Create an array to store region values

for i = 1:numel(numGroups)
    % Define the point coordinates
    x = normalizedAmplitude(i);
    y = normalizedEnergy(i);

    % Calculate the y-values on the lines
    y1 = -2*x + 1;
    y2 = -2*x + 2;

    % Check if the point is above one line and below the other
    if y < y1
        dividedDigits(i) = 1;  % Assign region value to the array
    elseif y >= y1 && y < y2
        dividedDigits(i) = 2; 
    else
        dividedDigits(i) = 3;  
    end
end

%% Importing data from all folders
clear;
dataFolder = 'data';
folderList = dir(dataFolder);
folderList = folderList([folderList.isdir]);  % Filter out non-directory entries
folderList = folderList(3:end);  % Exclude '.' and '..' entries
numFolders = numel(folderList);

% Initialize cell array to store audio data
audioDataArray = cell(numFolders, 1);

% Loop through each folder and import the audio files
for i = 1:numFolders
    folderPath = fullfile(dataFolder, folderList(i).name);
    fileList = dir(fullfile(folderPath, '*.wav'));
    audioDataCell = cell(numel(fileList), 1);
    
    % Loop through each .wav file in the folder and import it
    for j = 1:numel(fileList)
        filePath = fullfile(folderPath, fileList(j).name);
        [audioData, ~] = audioread(filePath);
        audioDataCell{j} = audioData;
    end
    fprintf('Folder number %d imported\n', i);
    audioDataArray{i} = audioDataCell;
end
disp('Finished importing!');
% %% Save all workspace variables
% save('workspace_data(1).mat');
% disp('Saved!');
% %% Load workspace variables
% clear;
% load('workspace_data(1).mat');
% disp('Loaded!');
 %% (5)
numFilesPerDigit = 50;  % Number of files per digit
numDigits = 10;  % Number of digits (change to compute less digits)
amplitudeSpectrum = cell(numDigits, 1);
%TODO: See if every window works 

for digit = 1:numDigits
    % Create cell
    amplitudeSpectrumCell = cell(numFilesPerDigit*numFolders, 1); 
    count = 1;
    for folder = 1:numFolders
        folderData = audioDataArray(folder);
        for j=1:numFilesPerDigit       
            % Get the current digit
            data = folderData{1,1}{(digit-1)*numFilesPerDigit+j, 1};
            
%             %Hamming window
             %window = hamming(numel(data));
            
            %Rectangular window
           % window = rectwin(numel(data));
%             
             %Blackman window
              %window = blackman(numel(data));
%             
%             %Hanning window
             window = hann(numel(data));

            % Compute amplitude spectrum
            spectrum = abs(fft(data .* window)) / numel(data);
      
            % Only keep positive frequencies
            spectrum = spectrum(1:floor(end/2));
            spectrum(2:end-1) = 2*spectrum(2:end-1);

            % Save the information
            amplitudeSpectrumCell{count} = spectrum;
            count=count+1;
            fprintf('Digit: %d | count: %d\n', digit, count);
        end
    end
    amplitudeSpectrum{digit} = amplitudeSpectrumCell;
end
disp('Finished computing!');
% %% Save all workspace variables
% save('workspace_data(2).mat');
% disp('Saved!');
% %% Load workspace variables
% clear;
% load('workspace_data(2).mat');
% disp('Loaded!');
 %% Getting the Median, Q1 and Q3
normalizedMedianValues = cell(numDigits, 1);
normalizedFirstQ = cell(numDigits, 1);
normalizedThirdQ = cell(numDigits, 1);

for digit = 1:numDigits
    % Define a function to calculate the length of an array inside a cell
    cellData = amplitudeSpectrum{digit, 1};
    getArrayLength = @(cellData) numel(cellData);

    % Apply the function to each cell in the amplitudeSpectrum cell
    lengths = cellfun(getArrayLength, cellData);

    % Get the maximum length
    maxArrayLength = max(lengths);
    
    % Create array of maxArrayLength for mean values
    normalizedMedianValues{digit} = zeros(maxArrayLength, 1);
    
    for itCell = 1:numel(cellData)
        fprintf('Digit: %d | Cell: %d\n', digit, itCell);
        paddedArray = vertcat(cellData{itCell, 1}, zeros(maxArrayLength - length(cellData{itCell, 1}), 1));
        normalizedMedianValues{digit} = normalizedMedianValues{digit} + paddedArray;
        
        currentArray = normalizedFirstQ{digit};
        currentArray = [currentArray, paddedArray]; %TODO: Solve this warning
        normalizedFirstQ{digit} = currentArray;
    end
    %Median
    normalizedMedianValues{digit} = cell2mat(arrayfun(@(row) prctile(normalizedFirstQ{digit}(row, :), 50), 1:size(normalizedFirstQ{digit}, 1), 'UniformOutput', false))';
    %Q3
    normalizedThirdQ{digit} = cell2mat(arrayfun(@(row) prctile(normalizedFirstQ{digit}(row, :), 75), 1:size(normalizedFirstQ{digit}, 1), 'UniformOutput', false))';
    %Q1
    normalizedFirstQ{digit} = cell2mat(arrayfun(@(row) prctile(normalizedFirstQ{digit}(row, :), 25), 1:size(normalizedFirstQ{digit}, 1), 'UniformOutput', false))';
    fprintf('Digit %d processed!\n', digit);
end
% %% Save all workspace variables
% save('workspace_data.mat(3)');
% disp('Saved!');
% %% Load workspace variables
% clear;
% load('workspace_data.mat(3)');
% disp('Loaded!');
%% Plot
fs = 48000;
% Define the number of rows and columns
numRows = 4;
numCols = 3;

% Create the figure
figure(6);

% Loop over each digit
for digit = 1:numDigits
    % Get the normalized mean values for the current digit
    mean = normalizedMedianValues{digit};
    firstQ = normalizedFirstQ{digit};
    thirdQ = normalizedThirdQ{digit};
    
    % Calculate the frequency values for the x-axis
    frequency = linspace(0, fs/2, numel(mean));
    
    % Limit the x-axis range to 0-8000 Hz
    frequencyRange = frequency <= 8000;
    frequencyLimited = frequency(frequencyRange);
    
    % Limit the y-axis range to 0-8000 Hz
    meanLimited = mean(frequencyRange);
    firstQLimited = firstQ(frequencyRange);
    thirdQLimited = thirdQ(frequencyRange);
    
    % Create subplots
    if digit == 1
        subplot(numRows, numCols, 2);
    else
        subplot(numRows, numCols, digit + 2);
    end
    
    % Plot the normalized mean values against the frequency
    plot(frequencyLimited, meanLimited);
    hold on;
    plot(frequencyLimited, firstQLimited, '--r', 'LineWidth', 1.5);
    plot(frequencyLimited, thirdQLimited, ':k', 'LineWidth', 1.5);
    
    % Add a title for each digit
    title(['Digit ', num2str(digit)]);
    
    % Add labels to the plot
    title(sprintf('Digit %d', digit-1));
    xlabel('Frequency(Hz)');
    ylabel('|X_d_f_t(f)|/N');
    legend('Median', 'Q1', 'Q3'); % Add legend
    
    % Add a grid to the plot
    grid on;
end

figure(6);
sgtitle('Meadian Amplitude Spectrum and First and Third Quantile');
set(gcf, 'Position', [100, 100, 1200, 800]);


%% The same as 5 but for STFT (6)

numFilesPerDigit = 50;  % Number of files per digit
folderChosen = 25;     % 1-25
fileChosen = 50;       % 1-50
numDigits = 10;        % Number of digits (change to compute fewer digits)

fs = 48000;            % Sampling rate
windowSize = 0.01*fs;     % Size of the window for each segment
hopSize = windowSize/2;         % Hop size (overlap) between consecutive windows
nfft = 2^9;           % Number of points for the fft


numRows = 4;
numCols = 3;

figure(7)
for digit = 1:numDigits
    % Get the folder
    folderData = audioDataArray(folderChosen);
    % Get the current digit
    data = folderData{1,1}{(digit-1) * numFilesPerDigit + fileChosen, 1};
    if digit==1
        subplot(numRows, numCols, 2);
    else
        subplot(numRows, numCols, digit + 2);
    end
    spectrogram(data, windowSize, hopSize, nfft, fs, 'yaxis');
    title(num2str(digit-1));  % Convert digit to string for title
end

colormap('jet');   % Set the colormap 
colorbar;          % Display the colorbar
sgtitle('Spectrogram of Audio Signal');  % Set the main plot title
set(gcf, 'Position', [100, 100, 1200, 800]);




    