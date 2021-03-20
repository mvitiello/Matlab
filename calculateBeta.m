function beta = calculateBeta(prices, benchmark, areprices)

% this function calculates the beta of the asset
% input = prices, benchmark and areprices, a 0 or 1 indicator variable
% output = beta of the asset

% if the input are prices we require a simple change in
% the lenght of the data

if areprices
    ret = [0;price2ret(prices)];
    retbench = [0;price2ret(benchmark)];
else
    ret = prices;
    retbench = benchmark;
end

% obtain standard deviation
volPrices = movstd(ret, length(prices)-1);
volBench = movstd(retbench, length(prices)-1);
% obtain correlation
correlation = movcorr(ret, retbench, length(prices)-1);
% obtain beta using two methods
betawd = correlation .* volPrices./volBench;
betawr = regress(retbench,ret);
% chose the minimum between two methods
beta = min(betawd, betawr);
end