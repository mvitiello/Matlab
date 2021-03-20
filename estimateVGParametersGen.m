function paramsTableStructure = estimateVGParametersGen(prices, path)

% function that generates the parameters of the process
% Variance Gamma (VG) from the financial time series
% we use the method in the paper
% ""Loregian, A., Mercuri, L., and Rroji, E., (2011),
% Approximation of the variance gamma model with a
% finite mixture of normals"".
% input
% prices = historical prices
% path = folder to save info
% output
% paramsTableStructure = Structure Table with n tables with
% the estimated Parameters

% get the number of tickers
m = size(prices, 2);

for i = 1 : m
    % prices and returns
    ret = price2ret(prices(:,m));
    % initial values: delta, mu, sigma, a
    mu = mean(ret);
    var = std(ret);
    skew = skewness(ret);
    kurt = kurtosis(ret,0)- 3;
    % Initial values by MoM
    f = @(x)calculateVGMoM(x,mu,var,skew,kurt);
    x0 = [0.002 -0.001 0.005 1.39];
    s = fsolve(f,x0);
    k1 = s(1);
    k2 = s(2);
    k3 = s(3);
    k4 = s(4);
    % initial values: delta, mu, sigma, a
    p0 = [k1 k2 k4 k3];
    % a rule of thumb is needed for the kurtosis
    if p0(4) > 10
        p0 = [k1 k2 k4 10];
    end
    % Laguerre polynomials
    % number of roots
    N = 75;
    % iterations
    steps = 500;
    % roots
    u = roots(LaguerrePoly(N));
    % evaluate the polynomial
    L = polyval(LaguerrePoly(N + 1), u);
    % weights
    w = u ./ ((N + 1)^2 .* L .^ 2);
    % obtain the parameters
    disp('Estimating Variance Gamma parameters');
    [a, mu, sigma, d] = parmemvg(ret, p0(1), p0(2), 1, p0(3),...
        p0(4), u, w, steps);
    % we save the parameters in a xls file
    prm = {'a'; 'mu'; 'sigma'; 'delta'};
    params = [a; mu; sigma; d];
    T = table(params, 'RowNames', prm);
    paramsTableStructure(i).Params = T;
    writetable(T, [char(path), '\', 'Params','_', ...
        datestr(today()), '.xls'],'Sheet', 1,...
        'WriteRowNames',1)
    % delete old portfolio params file
    if exist(char(strcat(path, '\', 'Params','_', ...
        datestr(today()-1), '.xls')), 'file')
        disp(['Deleting old portfolio Params file at ',...
            char(datetime('now','format','HH:mm'))])
        delete(char(strcat(path, '\', 'Params','_', ...
        datestr(today()-1), '.xls')))
    end
    
    % open Activex server
    e = actxserver('Excel.Application');
    % open file (enter full path!)
    ewb = e.Workbooks.Open([char(path), '\', 'Params','_', ...
        datestr(today()), '.xls']);
    % rename sheets
    ewb.Worksheets.Item(1).Name = 'VG_Params' ;
    % save to the same file
    ewb.Save
    ewb.Close(true)
    e.Quit
end
end
% obtain parameters by MoM
function F = calculateVGMoM(x,mu,var,skew,kurt)
F(1) = x(1) + x(2)*x(3) - mu;
F(2) = x(3)*x(4)^2 - var;
F(3) = (3*x(4)^2*x(2))/(sqrt(x(3))*(x(4)^3)) - skew;
F(4) = 3*(1 + 1/x(3)) - kurt;
end
% Approximation of VG pdf
function [a, m, sig, m0] = parmemvg(x, m00, mu0, beta, ...
    sig0, a0, u, w, passi)
% initial input
a = [a0; zeros(passi, 1)];
m = [mu0; zeros(passi, 1)];
m0 = [m00; zeros(passi, 1)];
sig = [sig0; zeros(passi, 1)];
phi = [psi(a0); zeros(passi, 1)];
T = length(x);
% calculate the parameters E y M step
for h = 2 : passi
    postDistributionValue = postdistrib(x, m0(h-1), ...
        m(h-1), beta, sig(h-1), a(h-1), u, w);
    phi(h) = sum(log(u)' * postDistributionValue) / T;
    % a parameter
    a(h) = real(invpsi(phi(h)));
    % drift of the vg process
    m0(h) = (((u.^-1)' * postDistributionValue * x) - ...
        T * sum(x)/sum(u' * postDistributionValue))/...
        (sum((u.^-1)' * postDistributionValue)- ...
        T^2/sum(u' * postDistributionValue));
    % mean parameter
    m(h) = (sum(x) - T * m0(h)) * beta/(sum(u' * postDistributionValue));
    % sigma parameter
    sig(h) = sqrt(1/T * sum(sum((((x * ones(1, length(u)))' - ...
        (u * ones(1, length(x)))* m(h)/beta - m0(h)).^2)./u .* ...
        postDistributionValue)));
end

% check the stimated parameters, if the estimated parameter is nan, we use
% the parameter obtained by MoM method

m2 = m;
m = m(passi);
if isnan(m) || isinf(m)
    m = m2(1);
end

sig2 = sig;
sig = sig2(passi);

if isnan(sig)|| isinf(sig)
    sig = sig2(1);
end

a2 = a;
a = a(passi);
if isnan(a)|| isinf(a)
    a = a2(1);
end

m02 = m0;
m0 = m0(passi);
if isnan(m0) || isinf(m0)
    m0 = m02(1);
end
end

function Y = invpsi(X)
% Y = INVPSI(X)
% Inverse digamma (psi) function.  The digamma function is the
% derivative of the log gamma function.  This calculates the value
% Y > 0 for a value X such that digamma(Y) = X.
% This algorithm is from Paul Fackler:
% http://www4.ncsu.edu/~pfackler/
L = 1;
Y = exp(X);
while L > 10e-8
    Y = Y + L*sign(X - psi(real(Y)));
    L = L / 2;
end
end

function y = postdistrib(x, d, m, beta, sig, a, u, w)
% aproximation of the convolution distribution (posterior distribution)
pesi = approxgammadens(a, u, w);
num = zeros(length(u), length(x));
y = zeros(length(u), length(x));
for i=1:length(x)
    num(:,i) = normpdf(x(i), d + m*u/beta, ...
        sig*sqrt(u/beta)) .* pesi;
    y(:,i) = num(:,i)/sum(num(:,i));
end
end

function y = approxgammadens(a, u, w)
% aproximation of the gamma density
y_num = u.^(a-1) .* w;
den = sum(y_num);
y = y_num/den;
end
