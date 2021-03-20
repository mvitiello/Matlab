function S = simulaVG(S0, k, sigma, mu, d, N, n, plotg)

% function that generates the VG process with desired parameters:
% S0 = Spot price
% sigma = volatility
% mu = drift
% d = the drift in the Brownian motion with drift
% k = the variance rate of the gamma time change
% N = number of observations
% n = number of simulations
% plotg, to plot the graphs

% initial parameters
a = 1/k;
b = 1/k;
T = N;
h = T/N;
X = zeros(N + 1, n);
I = zeros(N + 1, n);
% generate the process
for j = 1 : n
    for i = 1 : N
        I(i, j) = gamma1a(a * h, b);
        X(i + 1, j) = X(i, j) + d + mu*I(i, j) + sigma * sqrt(I(i, j))*randn;
    end
end
% X = X(1:n,1:N);
% X = X';
S = S0 .* exp(X);
if plotg
    % plot the process
    xdata = linspace(1, N, N);
    plot(xdata, S), title('VG process'),
    xlabel('time'), ylabel('S(t)'), xlim([0 N])
end
end

% function that generates the random number of a gamma process
function y = gamma1a(a, b)
if a == 0
    answ = 0;
elseif a <= 1
    answ = gamma1(a);
else
    answ = gamma2(a);
end
y = answ/b;
end

% Ahrens-Dieter’s Gamma generator (if a<=1)

function y = gamma1(a)
e = exp(1);
c = (a+e)/e;
flag = 0;
while flag == 0
    U1 = rand;
    U2 = rand;
    Y = c*U1;
    if Y <= 1
        Z = Y^(1/a);
        if U2 < exp(-Z)
            flag = 1;
        end
    else
        Z = -log((c-Y)/a);
        if U2 <= Z^(a-1)
            flag = 1;
        end
    end
end
y = Z;
end

%Fishman (Cheng and Feast)’s Gamma generator (if a>1)

function y = gamma2(a)
a2 = a-1;
c = (a-(1/(6*a)))/a2;
m = 2/a2;
d = m + 2;
flag = 0;
while flag == 0
    U1 = rand;
    U2 = rand;
    V = c*U2/U1;
    if m*U1 - d + V + (1/V) <= 0
        flag = 1;
    elseif m*log(U1) - log(V)+ V-1 <= 0
        flag = 1;
    end
end
y = a2*V;
end