clear;clc; close all
a = 0.9;
b = 2;
Vmax = 2;
R = 1;
x = linspace(-1,1);

v = Vmax.*(1-(1-a).*abs(x./R).^b);

plot(x,v)