function [bestSimul, ksresult,maxpvalue, hresult] = ...
    getBestSimulation(simulation,realserie)

% Function that select the best simulation using the 
% Kolgomorov - Smirnov (KS) test
% input:
% simulation = the simulation paths
% realserie = real prices
% output:
% bestSimul = the selected path with maximum p value
% ksresult = KS ratio of the selected path
% maxpvalue = p value of the selected variable
% hresult = desicion of the KS test 

% number of paths
nPaths = size(simulation, 2);

% analysing each serie 
for k = 1:nPaths
    simulation2check = simulation(:,k);
    % KS test
    [h,p,ks] = kstest2(real(price2ret(simulation2check)),...
        real(price2ret(realserie)),'Alpha',0.05);
    % desicion test
    hr(k) = h;
    % p value
    pr(k) = p;
    % ks ratio
    ksr(k) = ks;
end
% obtaning the maximum p value and other results
[maxpvalue, impvalue] = max(pr);
ksresult = ksr(impvalue);
hresult = hr(impvalue);
% obtaining the best simulation 
bestSimul = simulation(:,impvalue);
end