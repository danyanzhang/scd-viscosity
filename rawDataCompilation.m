function rawDataCompilation(filename)


%Data Compilation Code
%Athena Geisness
% Modified by Dan Zhang 09-26-2019
% Modified by Dan Zhang 03-27-2020

% Update as of 03-27-2020: Goal is to rewrite this from the ground up
% Too many uncertain gaps and weird spikes in the data, want to make sure this is all good
% Optimize for new experiments, don't worry about backwards compatibility from old scripts

% Modified by Dan Zhang 12-04-2020
% Removing the NaN lines of oxygen can result in situations where useful velocity data is cut off.
% Rolling back this functionality

% This script generates a warning for some of the table names
warning('off','all')

%% IMPORT DATA
% Import the MATLAB data
% .txt file: Results sfrom the MATLAB script for point detection
data = csvread([filename, '.txt']); % readmatrix is recommended as of R2019b
    % Columns
    % 1-avgvel: Average velocity of the blood flow
    % 2-intensity: Light intensity, may correlate to deoxyhemoglobin
    % 3-index:
    % 4-press: Flow driving air pressure
    % 5-FPS: Sampling rate per second of the image analysis, large deviations indicate skipped frames
    % 6-datetime: DateTime format of Excel form "737820.762224"
    % 7-_
    % 8-18:11 Binned velocities
% Should be 18 columns for 11 bins of velocity

% Import the Oxygen sensor data
% .csv file: Results from the oxygen sensor, separate code
dataOxy = readtable([filename, '.csv']);
% to call the things, use dataOxy.Time, dataOxy.Oxygen etc.
    % Columns
    % 1-Time: DateTime formatted as "1/30/2020  9:14:45 PM"
    % 2-Oxygen: Oxygen as a percentage of atmospheric
    % 3-TauPhaseMethod:
    % 4-SensorTemperature:
    % 5-AirPressure:

% CONVERT OXGYEN DATA INTO MATLAB TIME STAMPS BY INTERPOLATION
[oxygen_mat, cutRows] = ConvertOxyData(filename); % custom function, below

% Turn warnings back on
warning('on','all')


% Cutting out this section 12-02-2020
%{
%% CUT OUT ROWS THAT ARE NAN BECAUSE OF OXYGEN INTERPOLATION
data(cutRows, :) = [];
oxygen_mat(cutRows, :) = [];
% Display how many rows were removed
disp([num2str(length(cutRows)) ' rows of NaN were removed'])
%}

% Date of the experiment
formatOut = 'yyyy-mm-dd';
expDate = datestr(data(1,6),formatOut);

% Convert time of experiment into seconds
time = (data(:,6)-data(1,6))*86400; % time in seconds

% Extracting the time, velocity, pressure, index, intensity, and oxygen
velocity = data(:,1);
intensity = data(:,2);
index = data(:,3); 
pressure = data(:,4);
FPS = data(:,5);
Oxygen1 = oxygen_mat;


%% Extracted Binned velocity data
binVel1 = data(:,8);
binVel2 = data(:,9);
binVel3 = data(:,10);
binVel4 = data(:,11);
binVel5 = data(:,12);
binVel6 = data(:,13);
binVel7 = data(:,14);
binVel8 = data(:,15);
binVel9 = data(:,16);
binVel10 = data(:,17);
binVel11 = data(:,18);

CollectionDate = repmat(string(expDate),[size(data,1),1]);

% Make it into a new table
NewData = table(CollectionDate,time,Oxygen1,pressure,FPS,velocity,binVel1, binVel2, binVel3, binVel4, binVel5, binVel6, binVel7, binVel8, binVel9, binVel10,binVel11,intensity); 

%% 
saveasname = [filename, 'S']; % S for series
if exist('Time Series Data','dir')==7
    writetable(NewData,[pwd, '\Time Series Data\', saveasname, '.csv']);
else
    mkdir('Time Series Data')
    writetable(NewData,[pwd, '\Time Series Data\', saveasname, '.csv']);
end




% OXYGEN TIME-DATE CONVERSION SCRIPT
function [oxygen_mat, cutRows] = ConvertOxyData(filename)


    ploton = 0;
    
    % This script generates a warning for some of the table names
    warning('off','all')
    
    %% IMPORT DATA
    % Import the MATLAB data
    % .txt file: Results sfrom the MATLAB script for point detection
    data = csvread([filename, '.txt']); % readmatrix is recommended as of R2019b
    
    % Import the Oxygen sensor data
    % .csv file: Results from the oxygen sensor, separate code
    dataOxy = readtable([filename, '.csv']);
    
    % Turn warnings back on
    warning('on','all')
    
    %% CONVERT DATA TIME RANGE OF OXY TO MATCH MATLAB
    % General strategy:
    % Convert the date time format in the oxy file to serial datenum
    % 
    time_oxy = datenum(dataOxy.Time); % Time vector in serial datenum format
    oxygen_oxy = dataOxy.Oxygen; % Oxygen vector from sensor
    time_mat = data(:,6); % Excel stuff is already in datenum format
    oxygen_mat = interp1(time_oxy, oxygen_oxy, time_mat, 'linear'); % Linearly interpolated to the time from MATLAB
    % Values that are outside the range are NaN
    
    cutRows = find(isnan(oxygen_mat)); % rows that should be cut from the analysis because nan in O2
    
    if ploton==1
        figure
        plot(time_oxy, oxygen_oxy, 'x-')
        hold on
        plot(time_mat, oxygen_mat, 'x-')
        hold off
        legend('Oxygen Sensor', 'MATLAB Interpolated')
        xlabel('DateTime Format')
        ylabel('Oxygen Tension')
    end
    
    
    end % end of function







end % end of function 