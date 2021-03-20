function[portfMetrics, portfMetricsinMoney, rhoT, rhoG] =...
    CalculateVaRandESdinamicwithCopula(returns, Pnl, wtick, pnlperc)

% this function calculate the portfolio diversified VAR and ES with a copula
% input:
% returns = table with the returns of each ticker in the portfolio
% Pnl = money value of the portfolio
% wtick = expusore of each ticker
% output
% portfMetrics = Var, Es with t and gaussian copulas, in percentage
% portfMetricsinMoney = VaR and ES in monetary units
% rhoT = correlation matrix using t - copula
% rhoG = correlation matrix using gaussian copula

% size
nAssets = width(returns);
tdays = height(returns);
% NaN variables
mresiduals = NaN(tdays, nAssets);
mvariances = NaN(tdays, nAssets);
mfit = cell(nAssets, 1);
% model and options
model     = arima('AR', NaN, 'Distribution', 't', 'Variance', gjr(1,1));
opt   = optimoptions(@fmincon, 'Display', 'off', ...
    'Diagnostics', 'off', 'Algorithm', 'interior-point', 'TolCon', 1e-7);
disp(['Asymmetric GARCH model at ', char(datetime('now','format','HH:mm'))])
% stimate the GARCH model and fet the residuals
for i = 1:nAssets
    mfit{i} = estimate(model, returns{:,i}, 'Display', 'off', ...
        'Options', opt);
    [mresiduals(:,i), mvariances(:,i)] = infer(mfit{i}, ...
        returns{:,i});
end
% standarize residuals
% residuals are modeled as a standardized Student's t distribution
standardmresiduals = mresiduals ./ sqrt(mvariances);
disp(['Copula analysis at ', char(datetime('now','format','HH:mm'))]);
% decimal fraction allocated to each tail
tailFraction = 0.05;
marginal = cell(nAssets,1);
% modelling the distribution with pareto the tail and a kernel the body
for k = 1:nAssets
    marginal{k} = paretotails(standardmresiduals(:,k), ...
        tailFraction, 1 - tailFraction, 'kernel');
end
for k = 1:nAssets
    % transform each margin to uniform
    U(:,k) = marginal{k}.cdf(standardmresiduals(:,k));
end
% obtain the df
param = fitdist(price2ret(Pnl(2:end)), 'tlocationscale');
% find unity values
U(U >= 1) = 0.99;
% covariance
%for k = 15:15:size(returns,1)
c = cov(returns{1:end,1:end});
try
    % fit copulae
    [rhoT] = copulafit('t', U(:,:), 'Method', 'ApproximateML');
catch ME
    if (strcmp(ME.identifier,'stats:copulafit:RhoRankDeficient'))
        warning('The estimate of Rho has become rank-deficient. ',...
            'You may have too few data, or strong dependencies among variables.');
        % using Gaussian copula if the t copula rank is deficient
        disp(['So, using Gaussian Copula at ',...
            char(datetime('now','format','HH:mm'))]);
        [rhoT] = copulafit('Gaussian', U(:,:));
    end
end
[rhoG] = copulafit('Gaussian', U(:,:));
% VAR / ES with diversification
% VAR with t student without t-copula
portfvar = abs(sqrt((param.nu - 2)/param.nu))* tinv(tailFraction,...
    param.nu) .* sqrt(wtick' * c * wtick);
% ES with t student without t-copula
portfes =-((param.nu + tinv(tailFraction, param.nu)^2) / (param.nu - 1))...
    * abs(sqrt((param.nu - 2)/param.nu)) * ...
    (tpdf(tinv(tailFraction, param.nu), param.nu)/(tailFraction)) .*...
    sqrt(wtick' * c * wtick);
% VAR with t student with t-copula
portfvarwr = abs(sqrt((param.nu - 2)/param.nu))* tinv(tailFraction,...
    param.nu).* wtick' * real(sqrt(c + rhoT)) * wtick;
% ES with t student with t-copula
portfeswr =-((param.nu + tinv(tailFraction, param.nu)^2) / (param.nu - 1))...
    * abs(sqrt((param.nu - 2)/param.nu)) * ...
    (tpdf(tinv(tailFraction, param.nu), param.nu)/(tailFraction)) .*...
    sqrt(wtick' * (c + rhoT) * wtick);
% VaR with normal distribution and gaussian copula
portfvarwg = norminv(tailFraction) .* wtick' *...
    real(sqrt(c + rhoG)) * wtick;
% ES with normal distribution and gaussian copula
portfeswg = -(1/(tailFraction*sqrt(2*pi)))* exp(-0.5 *...
    norminv(tailFraction)^2).* wtick' * real(sqrt(c + rhoG))* wtick;
%end
% VaR/ES in Money
vares = repmat([portfvar', portfvarwr', portfvarwg', portfes', ...
    portfeswr', portfeswg'], length(Pnl), 1);

% standarize Var in money
pnlmu = repmat(mean(Pnl),length(Pnl), 1);
pnlstd = movstd(Pnl,15);
portfvarval = vares(:,1) .* (pnlmu + pnlstd .* pnlperc);
portfvarwrval = vares(:,2) .* (pnlmu + pnlstd .* pnlperc);
portfvarwgval = vares(:,3) .* (pnlmu + pnlstd .* pnlperc);
portfesval = vares(:,4) .* (pnlmu + pnlstd .* pnlperc);
portfeswrval = vares(:,5) .* (pnlmu + pnlstd .* pnlperc);
portfeswgval = vares(:,6) .* (pnlmu + pnlstd .* pnlperc);

% output
portfMetrics = vares;
portfMetricsinMoney = [portfvarval, portfvarwrval, portfvarwgval, ...
    portfesval, portfeswrval, portfeswgval];
end