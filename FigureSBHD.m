clear;clc;close all
%%
% For SBHD figure
filename = '202009CHC023';
a = fileread([filename, '.json']);
sample = jsondecode(a);
numTreatments = length(sample.experiments);



% Get the hill coefficients
coeffNative = hillFit(sample.experiments(1).O2, sample.experiments(1).vnorm);
coeffGBT = hillFit(sample.experiments(2).O2, sample.experiments(2).vnorm);

if length(sample.experiments) == 3
    coeffDMSO = hillFit(sample.experiments(3).O2, sample.experiments(3).vnorm);
end



%% Plot colors
plotcolor.aqua = [58	130	161]./255.*1.10;
plotcolor.gold = [239	201	112]./255;
plotcolor.purple = [106	80	121]./255;
plotcolor.lime = [170	213	116]./255;
plotcolor.red = [187	80	57]./255;
plotcolor.orange = [242	163	100]./255;
plotcolor.seafoam = [159	219	208]./255;
plotcolor.grey = [0.5, 0.5, 0.5];

%% Fill area of GBT effect
hold on
xplot = linspace(0, 21, 100);

fillx = [xplot, fliplr(xplot)];
filly = [hillEval(xplot, coeffNative), fliplr(hillEval(xplot, coeffGBT))];
fill(fillx, filly, [1, 240/255, 240/255])


% Plot the actual points
plot(sample.experiments(1).O2, sample.experiments(1).vnorm, 'b.','Color', plotcolor.aqua, 'MarkerSize', 16) % Native
plot(sample.experiments(2).O2, sample.experiments(2).vnorm, 'r.','Color', plotcolor.red, 'MarkerSize', 16) % GBT

if length(sample.experiments) == 3
    %plot(sample.experiments(3).O2, sample.experiments(3).vnorm, 'm.','MarkerSize', 16) % DMSO
end

%%
% Plot the hill fit lines
plot(xplot, hillEval(xplot, coeffNative), 'b-', 'Color', plotcolor.aqua, 'LineWidth', 2.5)
plot(xplot, hillEval(xplot, coeffGBT), 'r-', 'Color', plotcolor.red, 'LineWidth', 2.5)


%%
if length(sample.experiments) == 3
    %plot(xplot, hillEval(xplot, coeffDMSO), 'm-', 'LineWidth', 2)
end

hold off
%% Legend

axis([0, 22, 0, 1.04])

%{
if length(sample.experiments) == 3
    legend(sample.experiments(1).treatment, sample.experiments(2).treatment, sample.experiments(3).treatment, 'location', 'Southeast')
else
    legend(sample.experiments(1).treatment, sample.experiments(2).treatment, 'location', 'Southeast')
end
%{b.hospital}' % cell array
% vertcat(b.hosptal)
%}
%% Title
%title(filename)
xlabel('Oxygen (%)', 'FontSize', 10)
ylabel('Velocity Ratio', 'FontSize', 10)
saveas(gcf,[filename '.png'])

clinicalData = readtable('CBC Parameters.xlsx');
clinicalSample = clinicalData(sample.clinicalDataIdx,:);

if length(sample.experiments) == 3
    coeffs = vertcat(coeffNative, coeffGBT, coeffDMSO);
else
    coeffs = vertcat(coeffNative, coeffGBT);
end

coeffs = array2table(coeffs);
coeffs.Properties.VariableNames{1} = 'E0';
coeffs.Properties.VariableNames{2} = 'hill';
coeffs.Properties.VariableNames{3} = 'EC50';
coeffs.Properties.VariableNames{4} = 'Emax';

if numTreatments == 3
    treatment = {'Native'; 'GBT'; 'DMSO'};
else
    treatment = {'Native'; 'GBT'};
end
treatment = table(treatment);
results = [treatment, coeffs];

% Other table variables
V90 = zeros(numTreatments, 1);
bendPoint = zeros(numTreatments, 1);
V95 = zeros(numTreatments, 1);
VR90 = zeros(numTreatments, 1);
VR95 = zeros(numTreatments, 1);
BPL = zeros(numTreatments, 1);
BPU = zeros(numTreatments, 1);
inflection = zeros(numTreatments, 1);
slopeEC50 = zeros(numTreatments, 1);
slopeBP = zeros(numTreatments, 1);
slopeInflection = zeros(numTreatments,1);
AUC12 = zeros(numTreatments, 1);
AUC21 = zeros(numTreatments, 1);
etaOxy = zeros(numTreatments, 1);
etaDeoxy = zeros(numTreatments, 1);


for i = 1:numTreatments
    coeffTemp = results{i, 2:5}; % get the hill coefficients
    V90(i) = hillEvalY(0.9, coeffTemp);
    bendPoint(i) = ECempir(coeffTemp);
    V95(i) = hillEvalY(0.95, coeffTemp);
    VR90(i) = hillEvalY((1-coeffTemp(1)).*0.90 + coeffTemp(1), coeffTemp);
    VR95(i) = hillEvalY((1-coeffTemp(1)).*0.95 + coeffTemp(1), coeffTemp);
    [BPU(i), BPL(i)] = bendPoints(coeffTemp);
    inflection(i) = hillInflect(coeffTemp);
    slopeEC50(i) = hillDEval(coeffTemp(3), coeffTemp);
    slopeBP(i) = hillDEval(BPU(i), coeffTemp);
    slopeInflection(i) = hillDEval(inflection(i), coeffTemp);
    AUC12(i) = hillAUC(coeffTemp, 12);
    AUC21(i) = hillAUC(coeffTemp, 21);
    [etaOxy(i) etaDeoxy(i)] = getEta(filename, i);
end

clinicalDataIdx = table(repmat(sample.clinicalDataIdx,numTreatments,1));
clinicalDataIdx.Properties.VariableNames{1} = 'clinicalDataIdx';
clinicalSample = repmat(clinicalSample,numTreatments,1);

sampleID = table(repmat(filename,numTreatments,1));
sampleID.Properties.VariableNames{1} = 'sampleID';


results = [sampleID, clinicalDataIdx, results, table(V90), table(V95), ...
            table(VR90), table(VR95), table(BPL), table(BPU), ...
            table(inflection), table(slopeEC50), table(slopeBP), ...
            table(slopeInflection), table(AUC12), table(AUC21), table(etaOxy), table(etaDeoxy), clinicalSample, table(bendPoint)];

%%
%{
% Plot the bend points
hold on
plot(bendPoint(1), hillEval(bendPoint(1), coeffNative), 'ok')
plot(bendPoint(2), hillEval(bendPoint(2), coeffGBT), 'ok')

if length(sample.experiments) == 3
    plot(bendPoint(3), hillEval(bendPoint(3), coeffDMSO), 'ok')
end

hold off
%}

%% Set size of the plot
x0=50;
y0=50;
width=300;
height=280;
set(gcf,'position',[x0,y0,width,height])


%%
% Axis labels
xticks([0 2 4 6 8 10 12 14 16 18 20]);
xtick_mmHg = xticks./100.*760;
xtick_label = strtrim(cellstr(num2str(xtick_mmHg','%.0f'))');
xticklabels(xtick_label)
xlabel('Oxygen Tension (mmHg)')

ylabel('Normalized Velocity')
axis([0 13 0 1.1])
     
%%
% CUSTOM FUNCTIONS
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
    options = optimoptions('lsqnonlin', 'Display', 'None');
    fitcoeffs = lsqnonlin(fun, x0, lb, ub, options);
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
    BPU = coeff(3).*k.^(1./coeff(2));
    BPL = coeff(3).*(1./k).^(1./coeff(2));
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
function pAUC = hillAUC(coeff, B)
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);
    fun = @(x) (a-d) ./ (1 + (x./c).^b) + d;
    pAUC = integral(fun, 0, B) ./ B;
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

function ECany = ECF(coeffs, F)
    a = coeffs(1);
    b = coeffs(2); % hill
    c = coeffs(3); % EC50
    d = coeffs(4);
    ECany = (F./(100-F)).^(1./b).*c;
end


function ECval = ECempir(coeffs)
    empirECF = log10(coeffs(2)).*39.41 + 51.72
    ECval = ECF(coeffs, empirECF);
end
