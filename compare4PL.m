% Goal is to make sure the equations I'm using are absolutely correct
close all
clear
clc

upperAsy = 0;
lowerAsy = 100;
hill = 2;
EC50 = 4;

coeff = [upperAsy, hill, EC50, lowerAsy];
coeff2 = [lowerAsy, hill, EC50, upperAsy];

x = linspace(0, 21, 100);
y1 = hillDan(x, coeff);
y2 = hillSebaugh(x, coeff);

[BPU, BPL] = bendPoints(coeff)
BPUY = hillSebaugh(BPU, coeff);
BPLY = hillSebaugh(BPL, coeff);


figure
hold on
plot(x,y2)
plot(BPU, BPUY, 'x')
plot(BPL, BPLY, 'x')

%set(gca,'XScale','log')


hold off
legend('Dan', 'Sebaugh', 'Location', 'Southeast')

function y = hillDan(x, coeff)
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);
    y = (a-d).*x.^b ./ (c.^b + (x.^b)) + d;
end

function y = hillSebaugh(x, coeff)
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);
    y = (a-d) ./ (1 + (x./c).^b) + d;
end
% in sebaugh's definition
% a is the lower plateau
% d is the upper plateau
% b is the slope pfactor
% c is the midway


function [BPU, BPL] = bendPoints(coeff)
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);
    % Assuming Hill coefficients as follows
    % a = change
    % b = hill coefficient
    % c = starting point
    % d = half maximal point
    k = 4.6805; % for bend point equation
    % From Sebaugh 2003
    % Bend points are c*(k)^(1/b) and c*(1/k)^(1/b)
    BPU = c*k^(1/b);
    BPL = c*(1/k)^(1/b);
end

function y = hillEval(x, coeff)
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);
    y = (a-d) ./ (1 + (x./c).^b) + d;
end