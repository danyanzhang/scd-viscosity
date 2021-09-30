function plotTimeSeries(filename)
% Reads in CSV file from Sickle Cell Disease experiment and plots

[t, O2, P, v, vbin, FPS] = dataRead(filename);

% No FPS data in these files so only applying median filtering
% Median filter the rest of velocity and light intensity

% FPS filter

vmean = mean(vbin,2);
vmax = max(vbin,[],2);
if size(vbin,2)== 11
    skewness = mean(vbin(:,1:5),2)./mean(vbin(:,6:11),2);
elseif size(vbin,2)== 5
    skewness = mean(vbin(:,1:2),2)./mean(vbin(:,4:5),2);
else
    error('Unexpected number of velocity bins.')
end

FPSfilter = true;

if FPSfilter == true
    %figure
    %plot(FPS)
    % exclude anything outside 2 standard devs
    mu = mean(FPS);
    stdev = std(FPS);
    hiFPSidx = find(FPS>mu+1*stdev | FPS<mu-1*stdev);
    t(hiFPSidx) = [];
    O2(hiFPSidx) = [];
    P(hiFPSidx) = [];
    v(hiFPSidx) = [];
    vbin(hiFPSidx,:) = [];
    FPS(hiFPSidx) = [];
end



figure('Position', [100, 100, 1280, 720])
    subplot(3,1,1) % AVERAGE VELOCITY
plot(t,v,'k.')
ylabel('Velocity (um/s)')
axis([0, max(t), 0, nanmedian(v)*2])
title(filename(1:end-1),'Interpreter','None')

subplot(3,1,2) % OXYGEN TENSION
plot(t,O2,'k.')
ylabel('Oxygen (%)')
axis([0 max(t) 0 22])

subplot(3,1,3) % PRESSURE
plot(t,P,'k.')
ylabel('Pressure (PSI)')
axis([0 max(t) 0 max(P)+0.5])
xlabel('Time (s)')

end