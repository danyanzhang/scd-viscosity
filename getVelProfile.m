clear;clc

filename = '202103CHC028';

a = fileread([filename, '.json']);
sample = jsondecode(a);
numTreatments = length(sample.experiments);

% get Native file first
O2_points = sample.experiments(1).O2;
v_points = sample.experiments(1).vnorm;
b = sample.experiments(1).rawDataFile;
b = b{1}; % string with name of file

% all the vectors that I care about
[t, O2, P, v, vbin, FPS] = dataRead(b);

% need to get the timestamps to go between
c = fileread([b, '.json']);
timestamps = jsondecode(c);
tcoord = timestamps.tcoord;

indstart = findclosest(t,tcoord(3,1));
indend = findclosest(t,tcoord(3,2));

% Goal: find deoxy bluntness and oxy bluntness

%{
 % idealized Poisueille flow profile
    % vz = A/(4*eta) * (R^2 - r^2)
    R = 7.5;
    r = linspace(-R,R,100);
    vz = 2.*(R.^2 - r.^2);
    rnew = (r+7.5)./(15).*10 + 1;
    vznew = vz + 320;
    vznew2 = vz.*3.3 + 550;
    plot(rnew,vznew)
    plot(rnew,vznew2)
%}

% Show two samples that prove that eta1 is related to eta2

%{
eta1 = zeros(indend-indstart + 1,1);
eta2 = eta1;
for i = 1:length(eta1)
    [eta1(i), eta2(i)] = bluntness(vbin(i+indstart-1,:));
end
%}

% look at just vmean/vmax
eta_all = mean(vbin,2)./max(vbin,[],2);
eta = zeros(size(tcoord,1),1);
v_actual = eta;
for i = 1:size(tcoord,1)
    indstart = findclosest(t,tcoord(i,1));
    indend = findclosest(t,tcoord(i,2));
    eta(i) = median(eta_all(indstart:indend));
    v_actual(i) = median(v(indstart:indend));
end

% To insert spaces in between a series of things (very useful)
%x  = 1:4;
%x1 = cat( 1, x, zeros(size(x)) );
%x1 = transpose(x1(:));
O2_actual = O2_points';
O2_actual = cat(1, ones(size(O2_actual)).*21, O2_actual);
O2_actual = transpose(O2_actual(:));
O2_actual = O2_actual(1:end-2);
figure
yyaxis right
plot(O2_actual, v_actual,'o')
yyaxis left
plot(O2_actual, 1./eta,'o')

% Does voxelotor change viscosity response
% Possibility that this is purely related to velocity / shear rate

function [fitresult, gof] = fitVelProfile(x,y)
    [xData, yData] = prepareCurveData( x, y );
    
    % Set up fittype and options.
    ft = fittype( 'vmax*(1-(1-alpha)*abs((x-x0)/Whalf)^beta)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    % Whalf, alpha, beta, vmax, x0
    opts.Lower = [2 -Inf -Inf -Inf 2];
    opts.StartPoint = [5, 0.1869, 0.4898, max(y), 0.64];
    opts.Upper = [10 Inf Inf max(y)*2 12];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
end

function V = evalVelProfile(x,coeff)
    Whalf = coeff(1);
    alpha = coeff(2);
    beta = coeff(3);
    vmax = coeff(4);
    x0 = coeff(5);
    V = vmax.*(1-(1-alpha).*abs((x-x0)./Whalf).^beta);
end


function [eta1 eta2] = bluntness(testvprof)
    Vmax = max(testvprof);
    Vmax = 2;
    Vslip = 0.3;
    r = linspace(-1,1);
    R = 1;
    B = 2.159;

    V = Vmax.*( 1 - (1-Vslip).*(abs(r./R).^B));
    %hold on
    %plot(r,V);
    eta1 = max(testvprof)./mean(testvprof);

    x = linspace(1.5, 13.5, 11);
    [fitresult, gof] = fitVelProfile(x,testvprof);
    %plot(x,testvprof,'o')
    %hold on

    %xfit = linspace(0, 15);
    coeff=coeffvalues(fitresult);
    %vfit = evalVelProfile(xfit,coeff);

    %plot(xfit,vfit)
    %hold off

    % compare eta = alphabeta or vmean/vmax
    eta2 = (coeff(3)+1)/(coeff(3)+coeff(2));
end


function ind = findclosest(t,tcoord)
    [minValue,closestIndex] = min(abs(t-tcoord));
    ind = closestIndex;
end