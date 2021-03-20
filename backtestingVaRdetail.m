function [VaRsumUptable, TraffLightTable, ProportionOfFailuresTable,...
    ConditionalCoverageTable, TimeBFailuresTable]...
    = backtestingVaRdetail(ID, retseries, VaRseries, model, level)

% Function to obtain a backtesting table
% This function obtain VaR backtesting 
% The following test are done: 
% 1) Traffic light test,  Basel Compliance test
% 2) POF test, proportion of failures test, "classic" Kupice test
% 3) CCI test, conditional coverage independent test, Cristoffersen test  
% 4)TBFI test, time between failures independent test, Haas Test
% The input:
% ID = ticker identifier, name identifier
% retseries = returns to backtest the VaR
% VaRseries = value at risk estimated  
% model = model used to calculate the VaR: normal, t, cornish-fisher,
% variance gamma....
% level = confidence level
% Output, 4 tables, one of each test and a 5th table with a sum up 

% Initial tables
VaRsumUptable = table();
TraffLightTable = table();
ProportionOfFailuresTable = table();
ConditionalCoverageTable = table();
TimeBFailuresTable = table();
% number of tickers / VaR 
numticker = numel(ID);
% analysis for each ticker / VAR calculated 
for k = 1: numticker
    % get each return serie
    ret = retseries(:,k);
    retv = ret(ret~=0);
    % get each var seria 
    VaR = VaRseries(:,k);
    VaRv =  VaR(ret~=0); 
    % create object variable 
    vbt = varbacktest(retv,-VaRv,...
        'PortfolioID',ID(k),...
        'VaRID',model(k),...
        'VaRLevel',level);
    % runtest 
    vartest = runtests(vbt);
    % sum up table
    VaRsumUptable = [VaRsumUptable;vartest];
    trltable = tl(vbt);
    % 1) Traffic Light test
    TraffLightTable = [TraffLightTable;trltable];
    % 2) POF test
    propOfFailstable = pof(vbt,'TestLevel',level);
    ProportionOfFailuresTable = [ProportionOfFailuresTable; propOfFailstable];
    % 3) CCI test
    condCovtable = cci(vbt);
    ConditionalCoverageTable = [ConditionalCoverageTable;condCovtable];
    % 4) TBFI test 
    tbftable = tbfi(vbt);
    TimeBFailuresTable = [TimeBFailuresTable; tbftable];
end
end