coeffs = [0.57715, 0.47941, 3.0279, 1.169]; % GBT sample
%coeffs = [0.5292, 4.7857, 5.9612, 0.99995]; % Native sample
%coeffs = [0.89191, 0.27277, 17.168, 1.0901]; % weird GBT from 202003UMN022
% E0, hill, EC50, Emax
%coeffs = [0.8, 0.2, 3, 1];



x = linspace(0, 21, 500);
y = hillEval(x, coeffs);


EC10 = ECF(coeffs, 10)
EC20 = ECF(coeffs, 20)
EC30 = ECF(coeffs, 30)
EC40 = ECF(coeffs, 40)
EC50 = ECF(coeffs, 50)
EC60 = ECF(coeffs, 60)
EC70 = ECF(coeffs, 70)
EC80 = ECF(coeffs, 80)
EC90 = ECF(coeffs, 90)




figure
plot(x, y)
axis([0 21 0 1])

hold on

plot(EC10, hillEval(EC10, coeffs), 'o')
plot(EC20, hillEval(EC20, coeffs), 'o')
plot(EC30, hillEval(EC30, coeffs), 'o')
plot(EC40, hillEval(EC40, coeffs), 'o')
plot(EC50, hillEval(EC50, coeffs), 'x')
plot(EC60, hillEval(EC60, coeffs), 'o')
plot(EC70, hillEval(EC70, coeffs), 'o')
plot(EC80, hillEval(EC80, coeffs), 'o')
plot(EC90, hillEval(EC90, coeffs), 'o')
[BPU, BPL] = bendPoints(coeffs);
plot(BPU, hillEval(BPU, coeffs), 'x')
plot(BPL, hillEval(BPL, coeffs), 'x')



empirECF = log10(coeffs(2)).*39.41 + 51.72
ECval = ECF(coeffs, empirECF);

plot(ECval, hillEval(ECval, coeffs), '.', 'MarkerSize', 30)


hold off





% https://www.graphpad.com/support/faqid/169/


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