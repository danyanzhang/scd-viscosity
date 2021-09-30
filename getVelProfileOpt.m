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

% for all time points
x = linspace(1.5, 13.5, 11);
coeff = zeros(1000,5);
coeff2 = coeff;
tic

% 9.5 seconds
for i = 1:1000  % length(vbin)
    coeff(i,:) = fitVelProfile2(x,vbin(i,:));
end
toc

%{
% 68.8 seconds
tic
for i = 1:1000
    fitresult = fitVelProfile(x,vbin(i,:));
    coeff2(i,:)=coeffvalues(fitresult);
end
toc
%}

%{
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
%}

function drugFit1 = fit1DHill(cells)
    % Needs a cellProps object, outputs a drugSingle object
    y = cells.viability;
    x = cells.C1;
    fun = @(a) MuSyC1(a,x) - y; % where a are the coefficients
    x0 = [1, 1, 1, 1]; % initial guess
    lb = [0, 0, 0, min(x)]; % lower bound
    ub = [100, 100, 100, max(x)]; % upper bound
    %options = optimoptions('lsqnonlin','Display','None'); % disables text display
    options = optimoptions('lsqnonlin');
    fitcoeffs = lsqnonlin(fun,x0,lb,ub,options); % solves
    
    a = fitcoeffs;
    
    Emax = a(1);
    hill = a(2);
    E0 = a(3);
    EC50 = a(4);
    
    drugFit1 = drugSingle(Emax,hill,E0,EC50);
end

function coeff = fitVelProfile2(x,y)
    %fun = @(a) a(4).*(1-(1-a(2)).*abs((x-a(5))./a(1)).^a(3)) - y;
    fun = @(a) a(4).*(1-(1-a(2)).*abs((x)./a(1)+a(5)).^a(3)) - y;
    x0 = [5, 0.19, 0.51, max(y), 0.23]; % initial guess
    lb = [2 -Inf -Inf -Inf -8]; % lower bound
    ub = [10 Inf Inf max(y)*2 8]; % upper bound
    options = optimoptions('lsqnonlin','Display','None'); % disables text display
    %options = optimoptions('lsqnonlin');
    fitcoeffs = lsqnonlin(fun,x0,lb,ub,options); % solves

    coeff = fitcoeffs;
    
    Whalf = coeff(1);
    alpha = coeff(2);
    beta = coeff(3);
    vmax = coeff(4);
    x0 = coeff(5);
end


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