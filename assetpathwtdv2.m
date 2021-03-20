function [tpath, nu]  = assetpathwtdv2(prices, frecuency, n, m)

% function that generates m paths of t distributed prices
% input:
% m = numbers of paths
% n = number of days of each serie
% frecuency = differential time 
% prices = historical prices
% output:
% tpath = m*n matrix with the simulation paths 
% nu = dof of the t distrubution obtaining from modelling the returns 
% as a t distribution

% get returns 
ret = price2ret(prices);
% drift part
drift = (mean(ret) - 0.5 .* std(ret).^ 2).* frecuency .*...
    ones(length(ret) + n, m);
% stochastic part
param = fitdist(ret, 'tlocationscale');
nu = param.nu;
% dof 
if param.nu < 1
    nu = 2 + param.nu;
end
if param.nu < 2 && param.nu > 1
    nu = 1 + param.nu;
end 
% random numbers 
rr = trnd(nu, length(ret) + n, m);
% difussion 
difussion = (std(ret).*sqrt(frecuency)).* rr;
% both part are added and exponentialized
browniantProcess  = exp(drift + difussion);
spotPrice   = prices(1) .* ones(1, m);
simulation  = [spotPrice; browniantProcess];
tpath = cumprod(simulation);
end