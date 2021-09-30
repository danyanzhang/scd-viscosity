function lineIndexing(filename)

a = jsondecode(fileread([filename '.json'])); % JSON must exist
[t, O2, P, v, vbin, FPS] = dataRead(filename); % CSV must exist

plotTimeSeries(filename)

tcoord = a.tcoord;

% Needs t vector from original data, finds closest match
icoord = tcoord;
for i = 1:numel(tcoord)
    [~, idx] = min(abs(t-tcoord(i))); % get index of closest match
    icoord(i) = idx;
end

% Plots oxygen points on chart as a test
for i = 1:size(icoord,1)
    istart = icoord(i,1);
    iend = icoord(i,2);
    
    subplot(3,1,1)
    hold on
    plot(t(istart:iend), v(istart:iend), '.g')
    
    subplot(3,1,2)
    hold on
    plot(t(istart:iend), O2(istart:iend), '.g')

    subplot(3,1,3)
    hold on
    plot(t(istart:iend), P(istart:iend), '.g')
end


% Averages O2 values
% Use median because there will probably be some overlap
O2_SS = zeros(size(icoord,1),1);
for i = 1:size(icoord,1)
    istart = icoord(i,1);
    iend = icoord(i,2);
    O2_SS(i) = nanmedian(O2(istart:istart+10));
end
O2_SS = round(O2_SS);


% Get arrays for velocity
v_SS_all = cell(size(icoord,1),1);
O2_SS_all = v_SS_all;
v_SS = zeros(size(icoord,1),1);
for i = 1:size(icoord,1)
    istart = icoord(i,1);
    iend = icoord(i,2);
    v_SS_all{i} = v(istart:iend);
    v_SS(i) = nanmedian(v(istart:iend));
    O2_SS_all{i} = repmat(O2_SS(i), size(v_SS_all{i},1), size(v_SS_all{i},2));
end


figure
plot(O2_SS, v_SS, 'o')
hold on
for i = 1:size(icoord,1)
    plot(O2_SS_all{i}, v_SS_all{i},'k.')
end

figure
plot(O2_SS, v_SS, 'o')

end