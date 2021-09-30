function results = coeff2results(coeffs)
% Set this up to be able to take a matrix of coeffs as input? Not today

% coeffs is given in the old style
% Cdiff, Hill, C0, CP50
Cdiff = coeffs(1);
Hill = coeffs(2);
C0 = coeffs(3);
CP50 = coeffs(4);

% New order is
% E0, hill, EC50, Emax
E0 = C0/(Cdiff + C0);
hill = Hill;
EC50 = CP50;
Emax = 1;
coeffNew = [E0, hill, EC50, Emax]

numTreatments = 1;

figure
hold on
% Plot the hill fit lines
xplot = linspace(0, 21, 100);
plot(xplot, hillEval(xplot, coeffNew), 'b-', 'LineWidth', 2)
hold off

axis([0, 22, 0, 1.04])

xlabel('Oxygen (%)')
ylabel('Velocity Ratio')

coeffs = coeffNew;

coeffs = array2table(coeffs);
coeffs.Properties.VariableNames{1} = 'E0';
coeffs.Properties.VariableNames{2} = 'hill';
coeffs.Properties.VariableNames{3} = 'EC50';
coeffs.Properties.VariableNames{4} = 'Emax';

treatment = {'Native'};

treatment = table(treatment);
results = [treatment, coeffs];

% Other table variables
V90 = zeros(numTreatments, 1);
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


for i = 1:numTreatments
    coeffTemp = results{i, 2:5}; % get the hill coefficients
    V90(i) = hillEvalY(0.9, coeffTemp);
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
end


results = [results, table(V90), table(V95), ...
            table(VR90), table(VR95), table(BPL), table(BPU), ...
            table(inflection), table(slopeEC50), table(slopeBP), ...
            table(slopeInflection), table(AUC12), table(AUC21)];



            
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

end