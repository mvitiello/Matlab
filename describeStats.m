function tableStats = describeStats(financialserie)

% tableStats = describeStats(financialserie)
% obtains the mean, standard desviation, skewness and kurotsis of the
% financial serie, prices and returns. 
% financialserie (as input) must be a column vector of the prices 
%(then the returns are obtained) and tableStats (as output) 
% is a table with the measures described above.

% obtain returns
ret = price2ret(financialserie);
% describe stats from the financial serie
meanP = mean (financialserie,'omitnan');
standardDeviationP = std(financialserie,'omitnan' ); 
skewP = skewness(financialserie);
% flag = 0 to fix  the bias
kurtosisP = kurtosis(financialserie, 0) - 3; 
% describe stats from the returns
meanR = mean (ret,'omitnan');
standardDeviationR = std(ret,'omitnan' ); 
skewR = skewness(ret);
kurtosisR = kurtosis(ret, 0) - 3; 
% create the table
Stats = {'Mean'; 'Standard_Deviation'; 'Skewness';...
    'Exc. Kurtosis'};
Prices = [meanP; standardDeviationP; skewP; kurtosisP];
Returns = [meanR; standardDeviationR; skewR; kurtosisR];
tableStats   = table(Prices, Returns, 'RowNames', Stats);
end