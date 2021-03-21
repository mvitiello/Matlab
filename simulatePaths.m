function  [nPath, simulpath] = simulatePaths(date, histPrices, hisRet, n,...
    path2save, ticker)

% Function that simulate paths using montecarlo method. 
% A) first we predict n days using different models:
% Linear extrapolation and makima model extrapolation and 
% Geometric brownian motion and Variance Gamma model 
% for both models 100 paths are generated and a mean is obtained
% B) second, a entire path is simulated using different models:
% Geometric Brownian motion, a t distribution and a Variance Gamma model 
% For the models, 100 paths are generated and only one is choosed usign the
% Kolgomorov-Smirnov test with the highest p value
% A stats table for these paths is also calculated 
% And finally a plot is done with the paths and the distributions of these
% paths. 
% Inputs:
% date = time array
% historical prices for each ticker  
% historical returns for each ticker
% n = numbers of days to predict
% path2save = folder to save the info
% ticker = name, identifier for each instrument
% Output:
% nPath = table with predicted prices
% simulpath = table with entire path simulation 

% A) predict n days
% Simulate 100 paths
N = 100;
% grid and dates for extrapolation
x = linspace(1, length(histPrices), length(histPrices));
xq = linspace(length(histPrices)+1, length(histPrices) + n, n);
endDate = date(end);
xqDate = linspace(endDate+1, endDate + n, n);
% extrapolation
histpl = interp1(x, histPrices, xq, 'linear', 'extrap');
histpm = interp1(x, histPrices, xq, 'makima', 'extrap');
% brownian motion
stockpathBM = mean(reshape(assetpath(histPrices(end), mean(hisRet(hisRet~=0)),...
    std(hisRet(hisRet~=0)), 1, n, N, 'direct', '', 0),[], N), 2);
% Variance Gamma process
if exist([char(path2save), '\', 'Params','_', ...
        datestr(today()), '.xls'], 'file')
    paramsVG = readtable([char(path2save), '\', 'Params','_', ...
        datestr(today()), '.xls']);
    disp(['paramsVG at execution date already exist in ',...
        char(path2save),', we are reading that file at ',...
        char(datetime('now','format','HH:mm'))])
    stockpathVG = mean(simulaVG(histPrices(end), paramsVG.params(1),...
        paramsVG.params(3), paramsVG.params(2), ...
        paramsVG.params(4), n + 1, N, 0), 2);
else
    paramsVG = estimateVGParametersGen(histPrices, path2save);
    stockpathVG = mean(simulaVG(histPrices(end), paramsVG.Params.params(1),...
        paramsVG.Params.params(3), paramsVG.Params.params(2), ...
        paramsVG.Params.params(4), n + 1 , N, 0), 2);
end
% delete the yesterday file
if exist([char(path2save), '\', 'Params','_', ...
        datestr(today()-1), '.xls'], 'file')
    disp(['Deleting  yesterday paramsVG file in ', char(path2save),...
        ' at ', char(datetime('now','format','HH:mm'))])
    delete([char(path2save), '\', 'Params','_', ...
        datestr(today()-1), '.xls'])
end

% output table
nPath = timetable(xqDate', histpl', histpm',...
    stockpathBM(end-1:end), real(stockpathVG(end-1:end)));
nPath.Properties.VariableNames = {'PrExtrapLinear', 'PrExtrapMakimaMethod',...
    'PricesBrownianMotion', 'PricesVarianceGamma'};

% B) simulate entire path and predict price
% grid and dates
xt = linspace(1, length(histPrices) + n, length(histPrices) + n);
xtDate = linspace(date(1), endDate + n,...
    length(histPrices)+ n);
% brownian motion
stockpathBMt = assetpath(histPrices(1), mean(hisRet(hisRet~=0)),...
    std(hisRet(hisRet~=0)), 1, xt(end)-1, N, 'direct', '', 0);

% select the best brownian motion simulation
stockpathBMt = getBestSimulation(stockpathBMt, histPrices);
% t student 
[assetpathwt, nu] = assetpathwtdv2(histPrices, 1, 2, N);
% select the best t simulation
assetpathwt = getBestSimulation(assetpathwt, histPrices);

% Variance gamma
try
    stockpathVGt = simulaVG(histPrices(1), paramsVG.Params.params(1),...
        paramsVG.Params.params(3), paramsVG.Params.params(2), ...
        paramsVG.Params.params(4), xt(end), N, 0);
catch Me
    stockpathVGt = simulaVG(histPrices(1), paramsVG.params(1),...
        paramsVG.params(3), paramsVG.params(2), ...
        paramsVG.params(4), xt(end), N, 0);
end
% select the best VG simulation
stockpathVGt = getBestSimulation(stockpathVGt, histPrices);

% output table
tPath = timetable(xtDate',[histPrices;histpl'],...
    stockpathBMt, real(stockpathVGt(1:length(xtDate))), assetpathwt);
simulpath = timetable(xtDate', real(stockpathVGt(1:length(xtDate))),...
    assetpathwt);

% write output table 
tPath.Properties.VariableNames = {'RealPrices', 'PricesBrownianMotion',...
    'PricesVarianceGamma', 'PriceswithTstudent'};
writetimetable(tPath, [char(path2save), '\', 'SimulatedPrices', '.xls'], ...
    'filetype', 'spreadsheet')
% stats table 
tablerealp = describeStats(histPrices);
writetable(tablerealp, [char(path2save), '\', 'SimulatedPrices', '.xls'], ...
    'filetype', 'spreadsheet', 'sheet', 2, 'WriteRowNames', 1)
tableBmprices = describeStats(stockpathBMt);
writetable(tableBmprices, [char(path2save), '\', 'SimulatedPrices', '.xls'], ...
    'filetype', 'spreadsheet', 'sheet', 3, 'WriteRowNames', 1)
tableVGprices = describeStats(real(stockpathVGt));
writetable(tableVGprices, [char(path2save), '\', 'SimulatedPrices', '.xls'], ...
    'filetype', 'spreadsheet', 'sheet', 4, 'WriteRowNames', 1)
tableTprices = describeStats(assetpathwt);
writetable(tableTprices, [char(path2save), '\', 'SimulatedPrices', '.xls'], ...
    'filetype', 'spreadsheet', 'sheet', 5, 'WriteRowNames', 1)
% open Activex server
e = actxserver('Excel.Application');
% open file (enter full path!)
ewb = e.Workbooks.Open([char(path2save), '\', 'SimulatedPrices', '.xls']);
% rename sheets
ewb.Worksheets.Item(1).Name = 'Simulated_Prices' ;
ewb.Worksheets.Item(2).Name = 'Real_Prices_Stats' ;
ewb.Worksheets.Item(3).Name = 'BM_Prices_Stats' ;
ewb.Worksheets.Item(4).Name = 'VG_Prices_Stats' ;
ewb.Worksheets.Item(5).Name = ['T-',num2str(nu), '_df_Prices_Stats'];
% save to the same file
ewb.Save
ewb.Close(true)
e.Quit
% plot
g = figure('visible', 'off', 'units', 'normalized',...
    'outerposition', [0 0 1 1]);
% plot, legend, title...
% BM prices
subplot(2,4,1)
plot(xtDate(2:end), price2ret(stockpathBMt))
xlim([xtDate(1) xtDate(end)])
datetick('x', 'ddmmyyyy', 'keeplimits');
title(['Returns BM-Simulated from ', char(ticker), ' in ' ,...
    datestr(xtDate(1), 'dd-mm-yyyy'), ' to ',...
    datestr(xtDate(end), 'dd-mm-yyyy')])
% VG prices
subplot(2,4,2)
plot(xtDate(2:end), price2ret(real(stockpathVGt(1:length(xtDate)))))
xlim([xtDate(1) xtDate(end)])
datetick('x', 'mmyyyy', 'keeplimits');
title('VG-Returns ')
% T prices
subplot(2,4,3)
plot(xtDate(2:end), price2ret(assetpathwt))
xlim([xtDate(1) xtDate(end)])
datetick('x', 'ddmmyyyy', 'keeplimits');
title(['Returns t-Simulated with ', num2str(nu), ' df'])
% real prices
subplot(2,4,4)
plot(xtDate(2:end), price2ret([histPrices;histpl']))
xlim([xtDate(1) xtDate(end)])
datetick('x', 'ddmmyyyy', 'keeplimits');
title(['Real returns from ', char(ticker)])
% BM returns
subplot(2,4,5)
histogram(price2ret(stockpathBMt), 40);
title(['Histogram - BM-Returns from ', char(ticker), ' in ',...
    datestr(xtDate(1), 'dd-mm-yyyy'), ' to ',...
    datestr(xtDate(end), 'dd-mm-yyyy')]);
xlabel('value'), ylabel('frecuency')
% VG returns
subplot(2,4,6)
histogram(price2ret(real(stockpathVGt)), 40);
title('VG-Returns');
xlabel('value'), ylabel('frecuency')
% T returns 
subplot(2,4,7)
histogram(price2ret(assetpathwt), 40);
title(['t-returns with ', num2str(nu), ' df']);
xlabel('value'), ylabel('frecuency')
% real returns 
subplot(2,4,8)
histogram(price2ret([histPrices;histpl']), 40);
title(['Real Returns from ', char(ticker)]);
xlabel('value'), ylabel('frecuency')
% save plot 
saveas(g, char(strcat(path2save, '\', ...
    'SimulationPlot')), 'jpg')
end