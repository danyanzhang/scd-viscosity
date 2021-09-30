clear;clc;close all

filename = 'CHC028GS';
lineIndexing(filename)

timeSeriesMetadata = jsondecode(fileread([filename '.json'])); % JSON must exist
tcoord = timeSeriesMetadata.tcoord;

disp(timeSeriesMetadata.notes)

[t, O2, P, v, vbin, FPS] = dataRead(filename); % CSV must exist
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

v = medfilt1(v); % filter

% Needs t vector from original data, finds closest match
icoord = tcoord;
for i = 1:numel(tcoord)
    [~, idx] = min(abs(t-tcoord(i))); % get index of closest match
    icoord(i) = idx;
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
t_SS_all = v_SS_all;
for i = 1:size(icoord,1)
    istart = icoord(i,1);
    iend = icoord(i,2);
    v_SS_all{i} = v(istart:iend);
    t_SS_all{i} = t(istart:iend);
    v_SS(i) = nanmedian(v(istart:iend));
    O2_SS_all{i} = repmat(O2_SS(i), size(v_SS_all{i},1), size(v_SS_all{i},2));
end



% Check if second file exists
if isfile([pwd, filesep, 'Time Series Data', filesep, filename(1:end-1), '2', 'S','.csv'])
    disp('Second file exists')
else
end



% Attempt to interpolate based on the full data
figure
hold on
for i = 1:size(icoord,1)
    plot(t_SS_all{i}, v_SS_all{i},'.')
end
hold off

% Change anything 20 and above to be 21
O2_SS(O2_SS>=19) = 21;

idx21 = find(O2_SS==21); % gets indices of points that are 21% oxygen
v21_raw = vertcat(v_SS_all{idx21});
t21_raw = vertcat(t_SS_all{idx21});

v21_extrap = interp1(t21_raw, v21_raw, t, 'linear');

hold on
plot(t, v21_extrap, '-k')





figure
hold on
for i = 1:size(icoord,1)
    plot(t_SS_all{i}, v_SS_all{i},'.')
end
hold off


% Try doing the next point pair only
if O2_SS(1)~=21
    disp('First point is not 21')
end

% If the last ramp is 21 then ignore it
if O2_SS(end)==21
    disp('Last point was 21')
    O2_SS = O2_SS(1:end-1);
end

idx21 = find(O2_SS==21);
idxSS = find(O2_SS~=21);

v21 = v_SS(idx21);
vnorm = v_SS(idxSS)./v21;
O2norm = O2_SS(idxSS);

% add additional points
vnorm = [vnorm; 1];
O2norm = [O2norm; 21];

vnorm2 = [vnorm; 1; 1; 1];
O2norm2 = [O2norm; 19; 20; 21];

% if the last point is 21, then remove it?

% Also need to double check the hill fit between adding additional points in 19, 20


figure
plot(O2norm, vnorm, '-o')

coeffs = hillFit(O2norm, vnorm);
coeffs2 = hillFit(O2norm2, vnorm2);
[LBP, RBP] = bendPoints(coeffs);

hold on
plotHill(coeffs)
plotHill(coeffs2)

BPY = hillEval([LBP, RBP], coeffs);
plot([LBP, RBP], BPY, 'x-')
v90 = (1-coeffs(1));
plot(hillEvalY([0.9, 0.95], coeffs),[0.90, 0.95],'o')





% Plot the slope at the EC50 as well?
v = hillDEval(coeffs(3), coeffs);
yEC50 = (1-coeffs(1))/2 + coeffs(1);
b = hillEval(coeffs(3), coeffs) - v.*coeffs(3); % b = y-mx
xplot = linspace(0, coeffs(3)*2, 10); % 2X the EC50;
yplot = v.*xplot + b; 
plot(xplot, yplot, '-')

% Plot the slope at the inflection
inflection = hillInflect(coeffs);
inflecttest = hillInflect([100 1.405 0.01455^(-1/1.405) 0]);
inflectionY = hillEval(inflection, coeffs);
plot(inflection, inflectionY,'x');
inflectv = hillDEval(inflection, coeffs);
inflectb = inflectionY - inflectv.*inflection; % b = y-mx
yplotinflection = inflectv.*xplot + inflectb;
plot(xplot, yplotinflection, '--')

hold off

xlabel('Oxygen (%)')
ylabel('Velocity Ratio')
legend('Raw Points', 'Hill Fit', 'Hill Fit with 19,20','Point of Sickling','Location','Southeast')


pAUC = hillAUC(coeffs)
jsonencode(O2norm)
jsonencode(vnorm)

% Fit a Hill function with lsqnonlin
% HillFitFast(x,y)
% outputs coeffs = [a b c d]
% change, hill coefficient, starting asymptote, half-max
% Up to 4X faster than HillFunctionFit code using "fit"

function fitcoeffs = hillFit(x,y)
    % Use the 4PL formulation
    % a is left asymptote
    % b is hill coefficient
    % c is EC50
    % d is right asymptote
    fun = @(a) (a(1)-a(4)) ./ (1 + (x ./ a(3)).^a(2)) + a(4) - y; % where a are the coefficients
    x0 = [0.5, 2, 3, 1];
    lb = [0, 0, 0, 0];
    ub = [1.3, 20, 20, 1.3];
    options = optimoptions(@lsqnonlin, 'Display', 'None');
    fitcoeffs = lsqnonlin(fun, x0, lb, ub, options);
    %fitcoeffs = lsqnonlin(fun, x0, lb, ub);
end



function plotHill(coeffs)
    x = linspace(0, 21, 100);
    a = coeffs;
    y = (a(1)-a(4)) ./ (1 + (x ./ a(3)).^a(2)) + a(4);
    plot(x,y,'-')
end


function [BPU, BPL] = bendPoints(coeff)
    % Assuming Hill coefficients as follows
    % a = change
    % b = hill coefficient
    % c = starting point
    % d = half maximal point
    k = 4.6805; % for bend point equation
    % From Sebaugh 2003
    % Bend points are c*(k)^(1/b) and c*(1/k)^(1/b)
    BPU = coeff(3)*k^(1/coeff(2));
    BPL = coeff(3)*(1/k)^(1/coeff(2));
end


function y = hillEval(x, coeff)
    a = coeff;
    y = (a(1)-a(4)) ./ (1 + (x ./ a(3)).^a(2)) + a(4);
end

function x = hillEvalY(y, coeff)
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);
    x = ((a-d)./(y-d) - 1).^(1./b) .* c;
end


% Integration section
function pAUC = hillAUC(coeff)
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);
    fun = @(x) (a-d) ./ (1 + (x./c).^b) + d;
    pAUC = integral(fun, 0, 12) ./ 12;
end

function v = hillDEval(x, coeff)
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);
    v = -(b.*(a - d).*(x./c).^(b - 1)) ./ (c.*((x./c).^b + 1).^2);
end

function x = hillInflect(coeff)
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);
    x = ((b-1) ./ ((c.^-b).*(b + 1))).^(1./b);
end