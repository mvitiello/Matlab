function tableBacktesting = bcktsting2pnl(ret, vr, p, model, ticker)

% this function obtains a table with the results 
% from a VAR's backtesting:
% we calculate the Kupiec test, Cristofersen text and Haas test 
% the inputs are: 
% ret, the returns, a column n x 1 array,
% vr, VaR value, a  1 x m matrix    
% p,  percentile value, a single 1 x 1 value, e.g., 0.05 (5%)

% starting index
n  = size(ret, 1);
m = size(vr, 2);
q = size(vr,1);
I = zeros(n, m);
I1 = zeros(n, m);
I2 = zeros(n, m);
I3 = zeros(n, m);
I4 = zeros(n, m);

expecExceptions = p * n * ones(m , 1) ;
ticker = repmat(ticker, m, 1);
% obtaining exceptions

% this is in case the Var is a single value
if m == 1
    vr = vr .* ones(n,1);
end

if q == 1
    vr = vr .* ones(n,m);
end

% obtaining exceptions
for i = 1:size(I, 1)
    for j = 1:m
        if ret(i, :) < vr(i, j)
            I(i, j)= 1;
        end
    end
end

% number of exceptions 
numberExceptions  = sum(I);

% kupiec test

% failure ratio
failRatio = numberExceptions / n;
% kupiec numerator
kn = p.^ numberExceptions .* (1 - p) .^ (n - numberExceptions);
% kupiec denominator
kd = failRatio .^ numberExceptions .*...
    (1 - failRatio) .^ (n - numberExceptions);
% Kupiec Ratio
K = - 2 * log(kn ./ kd);

% Cristofersen test

%calculate n's and probabilities

for i = 2:size(I,1)
    for j = 1:m
        if (I(i-1,j) == 1 && I(i,j) == 1)
            I1(i,j)= 1;
        end
        if (I(i-1,j) == 1 && I(i,j) == 0)
            I2(i,j)= 1;
        end
        if (I(i-1,j) == 0 && I(i,j) == 1)
            I3(i,j)= 1;
        end
        if (I(i-1,j) == 0 && I(i,j) == 0)
            I4(i,j)= 1;
        end
    end
end

% states and probabilities
n11 = sum(I1);
n10 = sum(I2);
n01 = sum(I3);
n00 = sum(I4);

n0  = n00 + n01;
n1  = n10 + n11;
np  = n00 + n10;
npi = n01 + n11;

p0 = n01 / n0;
p1 = n11 / n1;
pp = npi / n;

% Cristofersen numerator
cn = (1 - pp).^ np .* pp .^ npi ;
% Cristofersen denomiantor
cd = (1 - p0).^ n00 .* p0 .^ n01 .* (1-p1) .^ n10 .* p1 .^ n11;
% Cristofersen ratio
C = -2 * log(cn ./ cd);

% Haas test 

Ht = zeros(m,1);

for jj = 1:m
    % Exceptions Positions
    k(jj).k = find(I(:,jj));
    % time Between Exceptions
    tbe(jj).tbe = diff(k(jj).k);
    LR(jj).LR = zeros(numberExceptions(jj),1);
    for j = 1:size(LR(jj).LR,1) - 1
        % Haas numerator
        hn(jj).hn(j, 1) = p .* (1 - p) .^ (tbe(jj).tbe(j, 1)- 1);
        % Haas denominator
        hd(jj).hd(j, 1) = tbe(jj).tbe(j, 1).^ -1 * (1 - 1 ./ ...
            tbe(jj).tbe(j, 1)) .^ (tbe(jj).tbe(j, 1) - 1);
        % Likelihood Ratio
        LR(jj).LR(j, 1) = - 2 * log(hn(jj).hn(j, 1) ./...
            hd(jj).hd(j, 1));
        H(jj).H = sum(LR(jj).LR);
        Ht(jj,1) = H(jj).H(1);
    end
end

% critical values and percentiles
cv1 = chi2inv(1 - p, 1) .* ones(m, 1);
cv2 = chi2inv(1 - p, numberExceptions').* ones(m, 1);
p = p .* ones(m, 1);


% create the table

Names = {'TickerID', 'ExpectedExceptions', 'Exceptions', 'VaRLevel', 'FailureRatio',...
    'KupiecRatio', 'CristofersenValue', 'CriticalValueCtest', 'Haas', ...
    'CriticalValueHtest', 'Model'};
tableBacktesting = table(ticker, expecExceptions, numberExceptions',...
    p, failRatio', K', C', cv1, Ht, cv2, model', 'VariableNames', Names);
end 