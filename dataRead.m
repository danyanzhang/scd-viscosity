function [t, O2, P, v, vbin, FPS] = dataRead(filename)

%% Read in data
Data = readtable([filename, '.csv']);
t = Data{:,2}; % time
O2 = Data{:,3};
P = Data{:,4};
FPS = Data{:,5};
v = Data{:,6};
vbin = Data{:, 7:17};
    
%v = medfilt1(v,4);
%O2 = medfilt1(O2,4);