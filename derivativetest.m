close all


%coeffs = [0.57715, 0.47941, 3.0279, 1.169]; % GBT sample
%coeffs = [0.5292, 4.7857, 5.9612, 0.99995]; % Native sample
%coeffs = [0.5, 6, 3, 1];
a = coeffs(1);
b = coeffs(2); % hill
c = coeffs(3); % EC50
d = coeffs(4);

x = linspace(0, 21, 200);
y = hillEval(x, coeffs);
y1 = derivative1(x, coeffs);
y2 = derivative2(x, coeffs);

figure
plot(x,y)

figure
plot(x, y1)

figure
plot(x, y2)

curvature = abs(y2)./((1 + y1.^2).^(3/2));
roc = 1./curvature;
figure
yyaxis left
plot(x, curvature)
yyaxis right
plot(x, roc)




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

function z = dwrtb(coeffs, x)
    a = coeffs(1);
    b = coeffs(2); % hill
    c = coeffs(3); % EC50
    d = coeffs(4);
    z = (d-a).*(x./c).^b.*log(x./c)./(1+(x./c).^b).^2;
end

function y1 = derivative1(x, coeffs)
    a = coeffs(1);
    b = coeffs(2); % hill
    c = coeffs(3); % EC50
    d = coeffs(4);
    num = -b.*(a-d).*(x./c).^b;
    dem = x.*((x./c).^b + 1).^2;
    y1 = num./dem;
end

function y2 = derivative2(x, coeffs)
    a = coeffs(1);
    b = coeffs(2); % hill
    c = coeffs(3); % EC50
    d = coeffs(4);
    num = b.*(a-d).*(x./c).^b.*(b.*(x./c).^b + (x./c).^b - b + 1);
    dem = x.^2 .* ((x./c).^b + 1).^3;
    y2 = num./dem;
end