function createJSONtcoord(filename)
close all
currDir = pwd;
saveDir = [pwd, filesep, 'Time Series Data'];

[t, O2, P, v, vbin, FPS] = dataRead(filename);
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

% Custom figure for doing the drawing
figure('Position', [100, 100, 1280, 720])
subplot(2,1,1)
yyaxis right
plot(t, O2, 'Color', [0.5, 0.5, 0.5]) % plot in grey
axis([0 max(t) -2 23])
ylabel('Oxygen (%)')

yyaxis left
plot(t, v, '.k')
axis([0, max(t), 0, nanmedian(v)*2])
ylabel('Velocity (um/s)')

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

subplot(2,1,2)
plot(t, P, '.k')
ylabel('Pressure (PSI')
xlabel('Time (s)')
axis([0 max(t) 0 max(P)+0.5])
subplot(2,1,1)


lines = drawSelection;

% Uses lines from drawSelection
tcoord = zeros(length(lines), 2); % start and end time
for i = 1:length(lines)
    tcoord(i,:) = lines{i}.Position(:,1);
end

a = struct;
a.timeSeries = filename;
a.notes = [];
a.tcoord = tcoord;

b = jsonencode(a);

% Save the file
fid = fopen([saveDir, filesep, filename, '.json'], 'w');
fprintf(fid, '%s', b);
fclose(fid);

pause(1)

doubleCheck = true;
if doubleCheck ==  true
    close all
    lineIndexing(filename)
end