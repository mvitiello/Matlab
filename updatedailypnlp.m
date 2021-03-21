function updatedailypnlp(pathfilewm, varargin)
tic
%% Individual Analysis
% Open, High, Low, Close Prices xls
xlsfileh = 'OHLData.xlsx';
xlsfilehpath = fullfile(pathfilewm, xlsfileh);
[~,nsheets] = xlsfinfo(xlsfilehpath);
% empty time table to fill it
AllHistoricalTable = timetable();
AllHistoricalRetTable = timetable();
% without varargin the prices are updated
if isempty(varargin)
    disp(['Reading and updating historical prices at ',...
        char(datetime('now','format','HH:mm'))])
    % update prices for each ticker
    for k = 1:length(nsheets)
        % read historical data
        historicalData = readtable(xlsfilehpath, 'filetype', ...
            'spreadsheet', 'sheet', k);
        % we 'add' weekend prices
        newtableprices = addweekendprices(historicalData);
        % read xls with last price and historical purchases
        xlsfileacc = 'Account.xlsx';
        xlsfileaccpath = fullfile(pathfilewm, xlsfileacc);
        % last price
        dailyinfo = readtable(xlsfileaccpath, 'filetype',...
            'spreadsheet', 'sheet', 2);
        % get tickers to update
        namestoupdate = matlab.lang.makeValidName(dailyinfo.Symbol);
        [~, updateIdx] = ismember(matlab.lang.makeValidName(newtableprices.Properties.VariableNames(1)),...
            namestoupdate);
        % ticker and name
        ticker = char(nsheets(k));
        [~, idName] = ismember(ticker, dailyinfo.Symbol);
        if idName > 0
            nameTicker = char(dailyinfo.Name(idName));
        elseif updateIdx > 0
            nameTicker = char(dailyinfo.Name(updateIdx));
        end
        % update price
        if iscell(dailyinfo.Open(updateIdx)) || ...
                iscell(dailyinfo.Low(updateIdx))
            updatedprice = timetable(newtableprices.Date(end)+1, ...
                dailyinfo.Last(updateIdx), ...
                str2double(dailyinfo.Open(updateIdx)),...
                dailyinfo.High(updateIdx),...
                str2double(dailyinfo.Low(updateIdx)));
        else
            updatedprice = timetable(newtableprices.Date(end)+1, ...
                dailyinfo.Last(updateIdx), dailyinfo.Open(updateIdx),...
                dailyinfo.High(updateIdx),dailyinfo.Low(updateIdx));
        end
        updatedprice.Properties.VariableNames = ...
            newtableprices.Properties.VariableNames;
        % warning message if the prices are updated
        if updatedprice.Time > datestr(today())
            warning(['Price from ',ticker, ' // ', nameTicker,...
                ' is already updated. Check. Maybe is an error']);
            newtablepricesupdated = newtableprices;
            tablehReturns = timetable(newtableprices.Date,...
                [0;price2ret(newtablepricesupdated{:,1})]);
        else
            disp(['Updating price from ',ticker, ' // ', nameTicker, ...
                ' at ', char(datetime('now','format','HH:mm'))]);
            newtablepricesupdated = [newtableprices; updatedprice];
            tablehReturns = timetable(newtableprices.Date,...
                price2ret(newtablepricesupdated{:,1}));
        end
        % returns table
        tablehReturns.Properties.VariableNames = ...
            updatedprice.Properties.VariableNames(1);
        % write table
        writetimetable(newtablepricesupdated, xlsfilehpath, ...
            'filetype', 'spreadsheet', 'sheet', k)
        AllHistoricalTable = [AllHistoricalTable, ...
            newtablepricesupdated(:, 1)];
        AllHistoricalRetTable = [AllHistoricalRetTable,  tablehReturns];
    end
    % with 'No update' the price are not updated
elseif strcmp(varargin(1), 'No Update')
    disp(['Reading historical prices at ', ...
        char(datetime('now','format','HH:mm'))])
    for k = 1:length(nsheets)
        % read historical data
        historicalData = readtable(xlsfilehpath, 'filetype', ...
            'spreadsheet', 'sheet', k);
        % we 'add' weekend prices
        newtableprices = addweekendprices(historicalData);
        % read xls with last price and with historical purchases
        xlsfileacc = 'Account.xlsx';
        xlsfileaccpath = fullfile(pathfilewm, xlsfileacc);
        % last price
        dailyinfo = readtable(xlsfileaccpath, 'filetype', ...
            'spreadsheet', 'sheet', 2);
        % get tickers to update
        names = matlab.lang.makeValidName(dailyinfo.Symbol);
        [~, Idx] = ismember(matlab.lang.makeValidName(newtableprices.Properties.VariableNames(1)),...
            names);
        tick = matlab.lang.makeValidName(newtableprices.Properties.VariableNames(1));
        if Idx
            lastPriceTick = newtableprices{:,tick};
        end
        lastPriceTickttable = timetable(newtableprices.Date,lastPriceTick);
        lastPriceTickttable.Properties.VariableNames(1) = tick;
        tablehReturns = timetable(newtableprices.Date(2:end-1),...
            price2ret(lastPriceTick(1:end-1)));
        tablehReturns.Properties.VariableNames(1) = tick;
        AllHistoricalTable = [AllHistoricalTable, lastPriceTickttable];
        AllHistoricalTable.Properties.DimensionNames(1) = {'Date'};
        AllHistoricalRetTable = [AllHistoricalRetTable,  tablehReturns];
        AllHistoricalRetTable.Properties.DimensionNames(1) = {'Date'};
    end
end
% read daily data and actual portfolio purchases
dailyinfoportf =  readtable(xlsfileaccpath, 'filetype', ...
    'spreadsheet', 'sheet', 1);
namestocharge =  intersect(unique(matlab.lang.makeValidName(dailyinfo.Symbol)),...
    unique(matlab.lang.makeValidName(dailyinfoportf.symbol)));
% time tables to fill it
AllRiskHistoricaltable = timetable();
AllRiskHistoricalRetTable = timetable();
tomorrowPriceticker = timetable();
VGPriceTable = timetable();
tDistributedTable = timetable();
% individual analysis & risk analysis
disp(['Individual analysis and risk analysis at ', ...
    char(datetime('now','format','HH:mm'))])
for h = 1: length(namestocharge)
    % read ticker
    ticker = namestocharge(h);
    % path to save
    path2save = fullfile(pathfilewm, ticker);
    if ~exist(char(path2save), 'dir')
        mkdir(char(path2save))
    end
    if ~exist(char(strcat(path2save, '\',...
            ticker,'_',datestr(today()),'.xlsx')), 'file')
        % historical data
        historicalData = AllHistoricalTable{:, ticker};
        tickerColumn = repmat(ticker, size(historicalData,1), 1);
        if AllHistoricalRetTable{1, ticker} == 0
            historicalRetDatat = [AllHistoricalRetTable{:, ticker}];
        else
            historicalRetDatat = [0; AllHistoricalRetTable{:, ticker}];
        end
        % pnl
        [~, pos] = ismember(ticker, ...
            matlab.lang.makeValidName(dailyinfoportf.symbol));
        datebuy = dailyinfoportf.Date(pos);
        [~, posbuy] = ismember(datebuy, AllHistoricalTable.Date);
        % change currency if it is necessary
        currency = dailyinfoportf.currency(pos);
        if ~strcmp(currency, 'EUR')
            historicalData = historicalData .*  ...
                AllHistoricalTable.USD_EUR_USDollarEuro;
        end
        % actual shares
        posbuyh = ismember(matlab.lang.makeValidName(dailyinfoportf.symbol),...
            ticker);
        datebuyh = dailyinfoportf.Date(posbuyh);
        datebuyhnum = datenum(datebuyh);
        buyhvec = ismember(AllHistoricalTable.Date, datebuyh);
        actualpos = dailyinfoportf.shares(posbuyh);
        buyhvecpos = zeros(length(buyhvec),1);
        % fees
        feespos = dailyinfoportf.fee(posbuyh);
        feesvec = zeros(length(buyhvec),1);
        % fees and shares history
        for j = posbuy:length(buyhvec)
            idx = cumsum(buyhvec);
            buyhvecpos(j) = actualpos(idx(j)).* buyhvec(j);
            feesvec(j) = feespos(idx(j)).* buyhvec(j);
        end
        % fees and shares acummulated
        shares = cumsum(buyhvecpos);
        feescum = cumsum(feesvec);
        % acumm pnl
        pnlData =  shares .* historicalData;
        vpnl = [0;diff(pnlData)];
        histData2ret = zeros(length(historicalData) ,1);
        histData2ret(posbuy:end) = historicalData(posbuy:end);
        % beta of the asset
        % identify index / benchmark
        bench = matlab.lang.makeValidName(dailyinfoportf.BenchSymbol(pos));
        historicalDataBench = AllHistoricalTable{:, bench};
        betaAsset = calculateBeta(historicalData,historicalDataBench, 1);
        % alpha, TE and Inforatio
        rethData = [0; price2ret(historicalData)];
        retbData = [0; price2ret(historicalDataBench)];
        alphavxs = rethData - retbData;
        % initial values
        inforatiov = zeros(length(rethData), 1);
        TrackingErrorv = zeros(length(rethData),1);
        % period to calculate the moving variables
        period = 14;
        for k = 1: length(rethData) - period
            [inforatiov(k:k+period),TrackingErrorv(k:k+period)] = ...
                inforatio(rethData(k:k+period), retbData(k:k+period));
        end
        % returns
        ret = [0; diff(histData2ret)./histData2ret(1:end-1)];
        ret = standardizeMissing(ret, inf);
        ret = fillmissing(ret, 'constant', 0);
        retacc = cumsum(ret);
        % pnl(0) - money invested (with fees) and return
        pnl0nf = cumsum(diff([0;shares]) .* buyhvec .* historicalData);
        fees0 = cumsum(diff([0;feescum]) .* buyhvec);
        pnl0 = pnl0nf + fees0;
        % middle price
        middlePrice = pnl0 ./ shares ;
        middlePrice = standardizeMissing(middlePrice, inf);
        middlePrice = fillmissing( middlePrice, 'constant', 0);
        retbyMiddlePrice = log(historicalData ./ middlePrice);
        retbyMiddlePrice = standardizeMissing(retbyMiddlePrice, inf);
        retbyMiddlePrice = fillmissing(retbyMiddlePrice, 'constant', 0);
        pnl0t = (retbyMiddlePrice .* pnl0);
        pnl0t = standardizeMissing(pnl0t, inf);
        pnl0t = standardizeMissing(pnl0t, -inf);
        pnl0t = fillmissing(pnl0t, 'constant', 0); %+ lossv;
        vpnlperc0t = retbyMiddlePrice ;
        vpnlperc0t = fillmissing(vpnlperc0t, 'constant', 0);
        % risk analysis
        accvol = movstd(retacc,5);
        % value at risk and expected shortfall with their backtesting
        % value at risk
        valueAtRiskt = calcParamVaR(historicalRetDatat,...
            't', 0.05, 5);
        % backtesting var
        [tickerSumUpt, tickerTLtablet, tickerPofTablet,...
            tickerCCITablet, tickerTBFITablet] =...
            backtestingVaRdetail(ticker, historicalRetDatat, valueAtRiskt,...
            {'t'}, 0.95);
        tableBacktestingVAR = table();
        tableBacktestingVAR.TickerID = tickerSumUpt.PortfolioID;
        tableBacktestingVAR.Model = tickerSumUpt.VaRID;
        tableBacktestingVAR.VaRLevel = tickerSumUpt.VaRLevel;
        tableBacktestingVAR.Failures =  tickerTLtablet.Failures;
        tableBacktestingVAR.Observations = tickerPofTablet.Observations;
        tableBacktestingVAR.TrafficLight =  tickerTLtablet.TL;
        tableBacktestingVAR.POF = tickerPofTablet.POF;
        tableBacktestingVAR.LRatioPOF = tickerPofTablet.LRatioPOF;
        tableBacktestingVAR.PValuePOF = tickerPofTablet.PValuePOF;
        tableBacktestingVAR.CCI = tickerCCITablet.CCI;
        tableBacktestingVAR.LRatioCCI = tickerCCITablet.LRatioCCI;
        tableBacktestingVAR.PValueCCI = tickerCCITablet.PValueCCI;
        tableBacktestingVAR.TBFI = tickerTBFITablet.TBFI;
        tableBacktestingVAR.LRatioTBFI = tickerTBFITablet.LRatioTBFI;
        tableBacktestingVAR.PValueTBFI = tickerTBFITablet.PValueTBFI;
        % write table
        writetable(tableBacktestingVAR, char(strcat(path2save, '\',...
            'Backtesting_', ticker, '.xlsx')), 'filetype', 'spreadsheet', ...
            'sheet', 1)
        % Expected shortfall
        expectedShortfallt = calcParamES(historicalRetDatat,...
            't', 0.05, path2save, 5);
        % backtesting ES
        tableBacktestingES = bcktstingexepshrtfall2pnl(historicalRetDatat,...
            historicalData, valueAtRiskt, expectedShortfallt,{'t'}, 0.05,...
            path2save);
        writetable(tableBacktestingES, char(strcat(path2save, '\',...
            'Backtesting_', ticker, '.xlsx')), 'filetype', 'spreadsheet', ...
            'sheet', 3)
        % VaR  with Cornish-Fisher Expansion using normal distribution
        cornishFisherVaRn = calcParamVaR(historicalRetDatat,...
            'CFn', 0.05, 5);
        % bcktesting VaR
        [tickerSumUpcf, tickerTLtablecf, tickerPofTablecf,...
            tickerCCITablecf, tickerTBFITablecf] =...
            backtestingVaRdetail(ticker, historicalRetDatat, cornishFisherVaRn,...
            {'cf'}, 0.95);
        tableBacktestingCFVAR = table();
        tableBacktestingCFVAR.TickerID = tickerSumUpcf.PortfolioID;
        tableBacktestingCFVAR.Model = tickerSumUpcf.VaRID;
        tableBacktestingCFVAR.VaRLevel = tickerSumUpcf.VaRLevel;
        tableBacktestingCFVAR.Failures =  tickerTLtablecf.Failures;
        tableBacktestingCFVAR.Observations = tickerPofTablecf.Observations;
        tableBacktestingCFVAR.TrafficLight =  tickerTLtablecf.TL;
        tableBacktestingCFVAR.POF = tickerPofTablecf.POF;
        tableBacktestingCFVAR.LRatioPOF = tickerPofTablecf.LRatioPOF;
        tableBacktestingCFVAR.PValuePOF = tickerPofTablecf.PValuePOF;
        tableBacktestingCFVAR.CCI = tickerCCITablecf.CCI;
        tableBacktestingCFVAR.LRatioCCI = tickerCCITablecf.LRatioCCI;
        tableBacktestingCFVAR.PValueCCI = tickerCCITablecf.PValueCCI;
        tableBacktestingCFVAR.TBFI = tickerTBFITablecf.TBFI;
        tableBacktestingCFVAR.LRatioTBFI = tickerTBFITablecf.LRatioTBFI;
        tableBacktestingCFVAR.PValueTBFI = tickerTBFITablecf.PValueTBFI;
        % write table
        writetable(tableBacktestingCFVAR, char(strcat(path2save, '\',...
            'Backtesting_', ticker, '.xlsx')), 'filetype', 'spreadsheet', ...
            'sheet', 2)
        % Cornish Fisher ES
        cornishFisherESn = real(calcParamES(historicalRetDatat,...
            'CFn', 0.05, path2save, 5));
        % backtesting ES
        tableBacktestingEScf = bcktstingexepshrtfall2pnl(historicalRetDatat,...
            historicalData, cornishFisherVaRn, cornishFisherESn, {'cf'}, ...
            0.05, path2save);
        writetable(tableBacktestingEScf, char(strcat(path2save, '\',...
            'Backtesting_', ticker, '.xlsx')), 'filetype', 'spreadsheet', ...
            'sheet', 4)
        % return adjusted by risk (Vol and ES)
        RetAdjVol = retacc ./ accvol;
        RetAdjVol = fillmissing(RetAdjVol, 'constant', 0);
        RetAdjESt = historicalRetDatat ./ expectedShortfallt;
        RetAdjESt = fillmissing(RetAdjESt, 'constant', 0);
        % tables
        RiskHistoricaltable = timetable(AllHistoricalTable.Date, tickerColumn, ...
            historicalData, ret, retacc, middlePrice, retbyMiddlePrice,...
            pnlData, vpnl, shares, pnl0, pnl0t, vpnlperc0t, ...
            betaAsset, accvol, alphavxs, inforatiov,TrackingErrorv, ...
            valueAtRiskt, expectedShortfallt, cornishFisherVaRn, ...
            cornishFisherESn, RetAdjVol, RetAdjESt);
        columnames = {'ID_', 'Hist_', 'Ret_', 'RetAcc_', 'MiddlePrice_', ...
            'RetbyMiddlePrice_', 'PNL_', 'Vpnl_', 'Shares_', 'PNL0_', 'PNL0T_', ...
            'VPNL0T_', 'Beta_' , 'VolAcc_', 'Alpha_',...
            'InfoRatio_','TE_', 'VAR_','ES_', 'VaRCF_', 'ESCF_',...
            'RetAdVol_', 'RetAdES_'};
        RiskHistoricaltable.Properties.VariableNames = strcat(columnames,...
            ticker);
        historicalRetData = timetable(AllHistoricalTable.Date,...
            historicalRetDatat);
        historicalRetData.Properties.VariableNames = strcat('hRet_',ticker);
        % write table
        disp(['Writing historical risk table from ', char(ticker), ' at ',...
            char(datetime('now','format','HH:mm'))]);
        writetimetable(RiskHistoricaltable, char(strcat(path2save, '\',...
            ticker,'_',datestr(today()),'.xlsx')), 'filetype', ...
            'spreadsheet', 'sheet', 1)
        % open Activex server
        e = actxserver('Excel.Application');
        % open file (enter full path!)
        ewb = e.Workbooks.Open(char(strcat(path2save, '\',...
            ticker,'_',datestr(today()),'.xlsx')));
        % rename sheets
        ewb.Worksheets.Item(1).Name = 'Main_Values' ;
        ewb.Worksheets.Item(2).Name = 'WorstDay' ;
        ewb.Worksheets.Item(3).Name = 'MetricAverages' ;
        % save to the same file
        ewb.Save
        ewb.Close(true)
        e.Quit
        % find the worst day
        [~, minRet] = min(ret);
        worstDayTable = RiskHistoricaltable(minRet,:);
        writetimetable(worstDayTable, char(strcat(path2save, '\',...
            ticker,'_',datestr(today()),'.xlsx')), 'filetype', ...
            'spreadsheet', 'sheet', 2)
        % calculate mean values for most important metrics
        averageMetricsTable = table();
        averageMetricsTable.Ticker = ticker;
        averageMetricsTable.averagewrtMiddlePrice = mean(retbyMiddlePrice, 'omitnan');
        averageMetricsTable.averageVAR = mean(valueAtRiskt, 'omitnan');
        averageMetricsTable.averageES = mean(expectedShortfallt, 'omitnan');
        averageMetricsTable.averageVARCF = mean(cornishFisherVaRn, 'omitnan');
        averageMetricsTable.averageESCF = mean(cornishFisherESn, 'omitnan');
        averageMetricsTable.averageVOL = mean(accvol, 'omitnan');
        writetable(averageMetricsTable, char(strcat(path2save, '\',...
            ticker,'_',datestr(today()),'.xlsx')), 'filetype', ...
            'spreadsheet', 'sheet', 3)
        % delete old ticker risk analysis
        if exist(char(strcat(path2save, '\',...
                ticker,'_',datestr(today()-1),'.xlsx')), 'file')
            disp('Deleting old portfolio file')
            delete(char(strcat(path2save, '\',...
                ticker,'_',datestr(today()-1),'.xlsx')))
        end
        
        % simulating paths
        [nFuturePrices, simulprices] = simulatePaths(AllHistoricalTable.Date,...
            historicalData, historicalRetDatat, 2, path2save, ticker);
        tticker = timetable(nFuturePrices.Time(1), ticker);
        futurepricetableticker = [tticker, nFuturePrices(1,:)];
        tomorrowPriceticker = [tomorrowPriceticker; futurepricetableticker];
        % store simulated prices
        % Variance - Gamma Prices
        VGprices = simulprices(:,1);
        tdistributedprices = simulprices(:,2);
        VGPriceTable2store = VGprices(1:height(AllHistoricalTable),:);
        columnamesVGtable = {'VGPrices_'};
        VGPriceTable2store.Properties.VariableNames = strcat(columnamesVGtable,...
            ticker);
        VGPriceTable = [VGPriceTable, VGPriceTable2store];
        % Prices t - distributed
        tPriceTable2store = tdistributedprices(1:height(AllHistoricalTable),:);
        columnamesTtable = {'tPrices_'};
        tPriceTable2store.Properties.VariableNames = strcat(columnamesTtable,...
            ticker);
        tDistributedTable = [tDistributedTable, tPriceTable2store];
        
        % plot variables for each ticker
        % last date
        endDate = datenum(RiskHistoricaltable.Time(end));
        % figure
        g = figure('visible', 'off', 'units', 'normalized',...
            'outerposition', [0 0 1 1]);
        % plot, legend, title...
        subplot(1,2,1)
        plotTimeSeries(RiskHistoricaltable.Time, historicalData,...
            ['Historical Prices from ', char(ticker), ' in ' ,...
            datestr(RiskHistoricaltable.Time(1), 'dd-mm-yyyy'), ' to ',...
            datestr(RiskHistoricaltable.Time(end), 'dd-mm-yyyy')], ...
            char(ticker), 'Date', 'Prices')
        % plot the data with minimum of 10 days of history
        if (endDate - datebuyhnum(1) > 10)
            subplot(1,2,2)
            idxzero = (retacc~=0);
            plotTimeSeries(RiskHistoricaltable.Time(idxzero), ...
                [rethData(idxzero), accvol(idxzero),...
                expectedShortfallt(idxzero), cornishFisherESn(idxzero)],...
                ['Returns and Risk from ', char(ticker), ' in ' ,...
                datestr(datebuyhnum(1), 'dd-mm-yyyy'), ' to ',...
                datestr(RiskHistoricaltable.Time(end), 'dd-mm-yyyy')],...
                {'Return','Volatiltiy', 'ESt', 'VaRnCF'}, 'Date', '%')
        end
        saveas(g, char(strcat(path2save, '\', ...
            'IndividualAnalysis')), 'jpg')
        % plot histogram pf returns
        gg = figure('visible', 'off', 'units', 'normalized',...
            'outerposition', [0 0 1 1]);
        subplot(1,2,1)
        plotTimeSeries(RiskHistoricaltable.Time(rethData~=0), ...
            rethData(rethData~=0),['Historical returns from ', ...
            char(ticker), ' in ' , datestr(RiskHistoricaltable.Time(1), ...
            'dd-mm-yyyy'), ' to ', datestr(RiskHistoricaltable.Time(end),...
            'dd-mm-yyyy')], char(ticker), 'Date', 'Returns')
        subplot(1,2,2)
        % distribution - fit returns
        h1 = histfit(rethData(rethData~=0), 40, 'normal');
        h1(2).LineWidth = 6;
        h1(2).LineStyle = '--';
        hold on
        h2 = histfit(rethData(rethData~=0), 40, 'tlocationscale');
        h2(2).Color = [.2 .2 .2];
        h2(2).LineStyle = ':';
        legend('returns','normal fit', 'returns', 't student fit')
        title(['returns and distribution fit from ', char(ticker), ' in ' ,...
            datestr(RiskHistoricaltable.Time(1), 'dd-mm-yyyy'), ' to ',...
            datestr(RiskHistoricaltable.Time(end), 'dd-mm-yyyy')]);
        xlabel('value'), ylabel('frecuency')
        saveas(gg, char(strcat(path2save, '\', ...
            'DistributionPlot_', char(ticker))), 'jpg')
        % plot returns with VaR/ ES
        gh = figure('visible', 'off', 'units', 'normalized',...
            'outerposition', [0 0 1 1]);
        subplot(1,2,1)
        plotTimeSeries(RiskHistoricaltable.Time(rethData~=0),...
            [rethData(rethData~=0),valueAtRiskt(rethData~=0), ...
            expectedShortfallt(rethData~=0)], ['Historical returns from ',...
            char(ticker), ' in ' , datestr(RiskHistoricaltable.Time(1), ...
            'dd-mm-yyyy'), ' to ', datestr(RiskHistoricaltable.Time(end), ...
            'dd-mm-yyyy'),', with t-VaR and t-ES'],{'returns','VARt', ...
            'ESt'},'Time','%')
        subplot(1,2,2)
        plotTimeSeries(RiskHistoricaltable.Time(rethData~=0),...
            [rethData(rethData~=0), cornishFisherVaRn(rethData~=0),...
            cornishFisherESn(rethData~=0)], ['Historical returns from ', ...
            char(ticker), ' in ' , datestr(RiskHistoricaltable.Time(1), ...
            'dd-mm-yyyy'), ' to ', datestr(RiskHistoricaltable.Time(end),...
            'dd-mm-yyyy'),', with CF Normal VaR and ES'], ...
            {'returns','VARnCF', 'ESnCF'}, 'Time', 'Returns')
        saveas(gh, char(strcat(path2save, '\', ...
            'ReturnswtVndESndVndESCF_', char(ticker))), 'jpg')
        % plot returns with VaR/ ES and compare Methods
        gj = figure('visible', 'off', 'units', 'normalized',...
            'outerposition', [0 0 1 1]);
        subplot(1,2,1)
        plotTimeSeries(RiskHistoricaltable.Time(rethData~=0),...
            [rethData(rethData~=0), valueAtRiskt(rethData~=0),...
            cornishFisherVaRn(rethData~=0)],['Historical returns from ',...
            char(ticker), ' in ' , datestr(RiskHistoricaltable.Time(1),...
            'dd-mm-yyyy'), ' to ', datestr(RiskHistoricaltable.Time(end),...
            'dd-mm-yyyy'),', with t-VaR and CF normal VaR'],{'returns',...
            'VARt', 'VARnCF'}, 'Time', '%')
        subplot(1,2,2)
        plotTimeSeries(RiskHistoricaltable.Time(rethData~=0),...
            [rethData(rethData~=0), expectedShortfallt(rethData~=0),...
            cornishFisherESn(rethData~=0)],['Historical returns from ', ...
            char(ticker), ' in ' , datestr(RiskHistoricaltable.Time(1), ...
            'dd-mm-yyyy'), ' to ', datestr(RiskHistoricaltable.Time(end),...
            'dd-mm-yyyy'),', with t-ES and CF normal ES'], {'returns',...
            'ESt', 'ESnCF'}, 'Time', '%')
        saveas(gj, char(strcat(path2save, '\', ...
            'CompareReturnswtVndESndVndESCF_', char(ticker))), 'jpg')
        
        % join all historical data in a single table
        AllRiskHistoricaltable = [AllRiskHistoricaltable,RiskHistoricaltable];
        AllRiskHistoricalRetTable = [AllRiskHistoricalRetTable, ...
            historicalRetData];
    else %
        disp([char(strcat(path2save, '\',...
            ticker,'_',datestr(today()),'.xlsx')), ' already exist'])
        RiskHistoricaltable = readtimetable(char(strcat(path2save, '\',...
            ticker,'_',datestr(today()),'.xlsx')),'filetype', ...
            'spreadsheet', 'sheet', 1);
        historicalData = RiskHistoricaltable{:,strcat('Hist_', char(ticker))};
        historicalRetDatat = RiskHistoricaltable{:,strcat('Ret_', char(ticker))};
        % simulating paths
        [nFuturePrices, simulprices] = simulatePaths(RiskHistoricaltable.Time,...
            historicalData, historicalRetDatat, 2, path2save, ticker);
        tticker = timetable(nFuturePrices.Time(1), ticker);
        futurepricetableticker = [tticker, nFuturePrices(1,:)];
        tomorrowPriceticker = [tomorrowPriceticker; futurepricetableticker];
        % store simulated prices
        % Variance - Gamma Prices
        VGprices = simulprices(:,1);
        tdistributedprices = simulprices(:,2);
        VGPriceTable2store = VGprices(1:height(AllHistoricalTable),:);
        columnamesVGtable = {'VGPrices_'};
        VGPriceTable2store.Properties.VariableNames = strcat(columnamesVGtable,...
            ticker);
        VGPriceTable = [VGPriceTable, VGPriceTable2store];
        % Prices t - distributed
        tPriceTable2store = tdistributedprices(1:height(AllHistoricalTable),:);
        columnamesTtable = {'tPrices_'};
        tPriceTable2store.Properties.VariableNames = strcat(columnamesTtable,...
            ticker);
        tDistributedTable = [tDistributedTable, tPriceTable2store];
        AllRiskHistoricaltable = [AllRiskHistoricaltable,RiskHistoricaltable];
    end
end

%% technically analysis
disp(['Technicall analysis at ', char(datetime('now','format','HH:mm'))])
for k = 1:length(nsheets)
    ticker = char(nsheets(k));
    % read historical prices and OHL prices
    historicalData = readtable(xlsfilehpath, 'filetype', ...
        'spreadsheet', 'sheet', k);
    historicalData.Properties.VariableNames = ...
        matlab.lang.makeValidName(historicalData.Properties.VariableNames);
    % period assumed 14 days
    period = 14;
    % technicall analysis table
    path2saveTA = fullfile(pathfilewm, 'TA', ticker);
    if ~exist(char(path2saveTA), 'dir')
        mkdir(char(path2saveTA))
    end
    if ~exist(char(strcat(path2saveTA, '\',...
            ticker,'_TA', '_', datestr(today()), '.xlsx')), 'file')
        [tablePPD, tableTechAnalysis] = Tecnicallyanalysis(historicalData,...
            period, path2saveTA, ticker);
        % write table
        disp(['Writing technical table from ', char(ticker), ' at ',...
            char(datetime('now','format','HH:mm'))]);
        rowb = 2;
        columnLetters = char(xlscol(rowb));
        cellReference = sprintf('%s2', columnLetters);
        writetable(tablePPD, char(strcat(path2saveTA, '\',...
            ticker,'_TA', '_', datestr(today()), '.xlsx')), ...
            'filetype', 'spreadsheet', 'WriteRowNames', 1, 'sheet', 1,...
            'Range', cellReference)
        row2b = rowb + 2 + height(tablePPD);
        cellReference2b = strcat(columnLetters, num2str(row2b));
        writetable(tableTechAnalysis, char(strcat(path2saveTA, '\',...
            ticker,'_TA', '_', datestr(today()), '.xlsx')), ...
            'filetype', 'spreadsheet', 'sheet', 1, 'Range', cellReference2b)
    else
        disp(['Technical table from ', char(ticker),', at ',...
            datestr(today()), ' already exist. At ', ...
            char(datetime('now','format','HH:mm'))])
    end
    % delete old file
    if exist(char(strcat(path2saveTA, '\',...
            ticker,'_TA', '_', datestr(today()-1), '.xlsx')), 'file')
        disp(['Deleting  old technical file in ', char(path2saveTA),...
            ' at ', char(datetime('now','format','HH:mm'))])
        delete(char(strcat(path2saveTA, '\',...
            ticker,'_TA', '_', datestr(today()-1), '.xlsx')))
    end
end

% the pnl file is created only if this file does not exist
if ~exist(char(strcat(pathfilewm, '\', 'Portfolio_PNL','_', ...
        datestr(today()) ,'.xlsx')), 'file')
    %% Individual Analysis Sum Up
    disp(['Writing Sum up table at ', char(datetime('now','format','HH:mm'))]);
    % A) Main Table
    % get buy prices
    tSymbolndPrice = table();
    tSymbolndPrice.symbol = dailyinfoportf.symbol;
    tSymbolndPrice.price =  dailyinfoportf.price .* dailyinfoportf.fx;
    [UniqueSymbol, idp] = unique(tSymbolndPrice.symbol);
    buyPrices = tSymbolndPrice.price(idp);
    table0sym = table(UniqueSymbol, buyPrices);
    otherUpdateIdx = ismember(matlab.lang.makeValidName(table0sym.UniqueSymbol),...
        namestocharge);
    table0sym = table0sym(otherUpdateIdx > 0, :);
    % delete pos USD/EUR
    [~, idusd] = ismember(table0sym.UniqueSymbol,...
        'USD/EUR - US Dollar Euro');
    table0sym(idusd > 0, :) = [];
    % get actual prices
    uniqueSymbol = unique(dailyinfoportf.symbol);
    symbolInDailyinfo = ismember(dailyinfo.Symbol, uniqueSymbol);
    lastPrice = dailyinfo.Last(symbolInDailyinfo);
    nsymbol = dailyinfo.Symbol(symbolInDailyinfo);
    tabletsym = table(nsymbol, lastPrice);
    tabletsym. Properties.VariableNames(1) = table0sym.Properties.VariableNames(1);
    tabletsym = sortrows(tabletsym,'UniqueSymbol');
    % delete pos USD/EUR
    [~, idusd] = ismember(tabletsym.UniqueSymbol,...
        'USD/EUR - US Dollar Euro');
    tabletsym(idusd > 0, :) = [];
    % table
    tableprices0t = join(table0sym, tabletsym);
    % change currency
    idxUsd = strcmp(dailyinfoportf.currency, 'USD');
    tickerUsd = unique(dailyinfoportf(idxUsd, {'Name' ,'symbol'}),'rows');
    idu = ismember(tableprices0t.UniqueSymbol, tickerUsd.symbol);
    USDEURv = idu * AllHistoricalTable.USD_EUR_USDollarEuro(end);
    USDEURv(USDEURv == 0) = 1;
    tableprices0t.lastPrice = tableprices0t.lastPrice .*  USDEURv;
    % get additional information
    [GroupsShares, Grnames] = findgroups(dailyinfoportf.symbol);
    tbl = table();
    % table with all information (including sold positions)
    tbl.UniqueSymbol = Grnames ;
    tbl.sharesBySec = splitapply(@sum, dailyinfoportf.shares,...
        GroupsShares);
    tbl.Accvalue = splitapply(@sum, dailyinfoportf.loc_quantity,...
        GroupsShares);
    % sold positions are not needed
    tableprices0t = innerjoin(tableprices0t, tbl);
    tableprices0t.Accvalue = tableprices0t.Accvalue;
    tableprices0t.Accvaluet = tableprices0t.sharesBySec .* ...
        tableprices0t.lastPrice;
    totalValuePortfT = sum(tableprices0t.Accvaluet);
    totalValuePortf = sum(tableprices0t.Accvalue);
    % percentage by each ticker
    tableprices0t.PercByTickt = (tableprices0t.Accvalue) ./ totalValuePortf;
    tableprices0t.PercByTickT = tableprices0t.Accvaluet ./ ...
        totalValuePortfT;
    % return using price
    tableprices0t.PercUsingPrice = tableprices0t.lastPrice ./...
        tableprices0t.buyPrices - 1;
    tableprices0t.PercValUsingPrice = tableprices0t.PercUsingPrice ...
        .* tableprices0t.Accvalue;
    % get beta
    names2filteridxBeta = startsWith(AllRiskHistoricaltable.Properties.VariableNames,...
        'Beta_');
    names2filterBeta = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxBeta);
    names2filterBeta(strcmp(names2filterBeta, 'Beta_USD_EUR_USDollarEuro')) = [];
    tableprices0t.Beta = AllRiskHistoricaltable{end, names2filterBeta}';
    % exit price wrt buy price
    tableprices0t.ExitPriceUp = 1.8 .* tableprices0t.buyPrices;
    tableprices0t.ExitPriceDown = 0.6 .* tableprices0t.buyPrices;
    % middle price table
    tablepricesMiddlePrice = table();
    tablepricesMiddlePrice.Name = tableprices0t.UniqueSymbol;
    % middle price
    tablepricesMiddlePrice.MidPrice = tableprices0t.Accvalue ./ ...
        tableprices0t.sharesBySec;
    % return by middle price
    tablepricesMiddlePrice.PercByMidPr = tableprices0t.lastPrice ./ ...
        tablepricesMiddlePrice.MidPrice - 1;
    % Money Weigthed Return - MWR, Time Weigthed Return - TWR
    [TWRt, MWRt] = calculateWR(tablepricesMiddlePrice.Name,...
        AllRiskHistoricaltable, dailyinfoportf, tabletsym);
    tablepricesMiddlePrice.PercByTick = tableprices0t.PercByTickt;
    tablepricesMiddlePrice.TWR = TWRt.twr;
    tablepricesMiddlePrice.MWR = MWRt.mwr;
    tablepricesMiddlePrice.PercValMidPr = tablepricesMiddlePrice.PercByMidPr .*...
        tableprices0t.Accvalue;
    % exit price wrt middle price
    tablepricesMiddlePrice.ExitMdPriceUp = 1.8 .* tablepricesMiddlePrice.MidPrice;
    tablepricesMiddlePrice.ExitMdPriceDown = 0.6 .* tablepricesMiddlePrice.MidPrice;
    % B) Main values with dividends and fees
    tablepricesDivndFee = table();
    tablepricesDivndFee. Name = tablepricesMiddlePrice.Name;
    tablepricesDivndFee. Accvalue = tableprices0t.Accvalue;
    % fee
    fee = splitapply(@sum, dailyinfoportf.fee,...
        GroupsShares);
    tblf = table();
    tblf.Name = Grnames;
    tblf.fee = fee;
    % we do not need fee from sold positions
    tablepricesDivndFee = innerjoin(tablepricesDivndFee,tblf);
    tablepricesDivndFee. AccvalueNofee = tablepricesDivndFee. Accvalue...
        - tablepricesDivndFee.fee ;
    tablepricesDivndFee. Percfee = 100*(tablepricesDivndFee. fee  ./ ...
        tablepricesDivndFee. Accvalue);
    % dividends
    dividends =  readtable(xlsfileaccpath, 'filetype', 'spreadsheet',...
        'sheet', 4);
    % search exchange price
    for otherUpdateIdx = 1: height(dividends)
        if strcmp(dividends.currency(otherUpdateIdx),'USD')
            dividends.ExhType(otherUpdateIdx) = ...
                AllHistoricalTable.USD_EUR_USDollarEuro(dividends.Date(otherUpdateIdx));
            dividends.Valuelm(otherUpdateIdx) = ...
                dividends.value(otherUpdateIdx) .* ...
                dividends.ExhType(otherUpdateIdx);
        else
            dividends.ExhType(otherUpdateIdx) = 1;
            dividends.Valuelm(otherUpdateIdx) = ...
                dividends.value(otherUpdateIdx).* ...
                dividends.ExhType(otherUpdateIdx);
        end
    end
    [GroupsDividens, Gnames] = findgroups(dividends.symbol);
    dividensTable = table();
    dividensTable.Name = Gnames;
    % accumulated dividends by ticker
    dividensTable.AccValDiv = splitapply(@sum, dividends.Valuelm,...
        GroupsDividens);
    % insert dividend
    for j = 1:height(tablepricesDivndFee)
        ticker = tablepricesDivndFee.Name(j);
        [h,iddiv] = ismember(ticker, dividensTable.Name);
        if h > 0
            tablepricesDivndFee. Div(j) = dividensTable.AccValDiv(iddiv);
        else
            tablepricesDivndFee. Div(j) = 0;
        end
    end
    tablepricesDivndFee. PercDiv = tablepricesDivndFee. Div ./ ...
        tablepricesDivndFee.Accvalue;
    % pnl + dividend
    tablepricesDivndFee. Accvaluetwdiv = tablepricesDivndFee. Div +...
        tableprices0t.Accvaluet ;
    % return with dividend and fee
    tablepricesDivndFee.PercTickUsingPnl = tablepricesDivndFee. Accvaluetwdiv...
        ./ tablepricesDivndFee.Accvalue - 1;
    % to calculate TWR with dividends and fee
    [TWRtwd, MWRtwd] = calculateWRwDivnfee(tablepricesMiddlePrice.Name,...
        AllRiskHistoricaltable,dividends, dailyinfoportf, tabletsym);
    tablepricesDivndFee.TWRtwd = TWRtwd .twr;
    tablepricesDivndFee.MWRtwd = MWRtwd.mwr;
    
    % B.1) table with info divided asset type (Crypto, ETF,...)
    % find groups in the new table
    [GroupsAssetType, GrAssetType] = findgroups(dailyinfoportf.Category1);
    tblAssetType = table();
    % table with all information (including sold positions)
    tblAssetType.AssetType = GrAssetType ;
    ValueAccAssetType = splitapply(@sum, dailyinfoportf.loc_quantity,...
        GroupsAssetType);
    tblAssetType.ValueAccAssetType = ValueAccAssetType;
    % get type for each ticker
    dailyinfoportf.AssetTypeID = GroupsAssetType;
    [~,tickeridxpos] = ismember(tablepricesDivndFee.Name, dailyinfoportf.symbol);
    valuetbyassetType = table();
    valuetbyassetType.Name = tablepricesDivndFee.Name;
    valuetbyassetType.Accvaluetwdiv = tablepricesDivndFee. Accvaluetwdiv;
    valuetbyassetType.AssetType = dailyinfoportf.Category1(tickeridxpos);
    valuetbyassetType.AssettypeID = dailyinfoportf.AssetTypeID(tickeridxpos);
    % table for each asset type
    tblAssetType.PercbyAssetType = 100*(tblAssetType.ValueAccAssetType /...
        sum(tblAssetType.ValueAccAssetType));
    tblAssetType.ValuetAccAssetType = splitapply(@sum, valuetbyassetType.Accvaluetwdiv,...
        valuetbyassetType.AssettypeID);
    tblAssetType.RetbyAssetType = (tblAssetType.ValuetAccAssetType ./ ...
        tblAssetType.ValueAccAssetType) - 1;
    % B.2) table with info divided exchange (Crypto, StockExchange)
    [GroupsExchange, GrExch] = findgroups(dailyinfoportf.Category2);
    tblExch = table();
    % table with all information (including sold positions)
    tblExch.Exchange = GrExch ;
    ValueAccExchange = splitapply(@sum, dailyinfoportf.loc_quantity,...
        GroupsExchange);
    tblExch.ValueAccExch = ValueAccExchange;
    % get type for each ticker
    dailyinfoportf.ExchangeTypeID = GroupsExchange;
    [~,tickeridxposexch] = ismember(tablepricesDivndFee.Name, dailyinfoportf.symbol);
    valuetbyExch = table();
    valuetbyExch.Name = tablepricesDivndFee.Name;
    valuetbyExch.Accvaluetwdiv = tablepricesDivndFee. Accvaluetwdiv;
    valuetbyExch.Exchange = dailyinfoportf.Category2(tickeridxposexch);
    valuetbyExch.ExchtypeID = dailyinfoportf.ExchangeTypeID(tickeridxposexch);
    % table for each asset type
    tblExch.PercbyExchange = 100*(tblExch.ValueAccExch /...
        sum(tblExch.ValueAccExch));
    tblExch.ValuetbyExch = splitapply(@sum, valuetbyExch.Accvaluetwdiv,...
        valuetbyExch.ExchtypeID);
    tblExch.RetbyExhange = (tblExch.ValuetbyExch ./ ...
        tblExch.ValueAccExch) - 1;
    
       % C) Sector, Country
    tblSectCountr = unique(dailyinfoportf(:, {'symbol', 'Sector', ...
        'Country'}), 'rows');
    % add sector and country
    [~,tickSectContr] = ismember(tablepricesDivndFee.Name, tblSectCountr.symbol);
    % table with country and sector
    valuetbySector = table();
    valuetbySector.Name = tablepricesDivndFee.Name;
    valuetbySector.Accvaluefee = tablepricesDivndFee.Accvalue;
    valuetbySector.Accvaluetwdiv = tablepricesDivndFee.Accvaluetwdiv;
    valuetbySector.Sector = tblSectCountr.Sector(tickSectContr);
    valuetbySector.Country = tblSectCountr.Country(tickSectContr);
    % group by sector and country
    [GroupSector, GrSector] = findgroups(valuetbySector.Sector);
    [GroupCountry, GrCountry] = findgroups(valuetbySector.Country);
    valuetbySector.SectorID  = GroupSector;
    valuetbySector.CountryID = GroupCountry;
    % table by sector
    tblSector = table();
    tblSector.Sector = GrSector;
    tblSector.ValuebySector = splitapply(@sum, valuetbySector.Accvaluefee,...
        valuetbySector.SectorID);
    tblSector.PercbySector = 100*(tblSector.ValuebySector ./ ...
        sum(tblSector.ValuebySector));
    tblSector.ValuetbySector = splitapply(@sum, valuetbySector.Accvaluetwdiv,...
        valuetbySector.SectorID);
    tblSector.RetbySector = (tblSector.ValuetbySector ./...
        tblSector.ValuebySector) - 1;
    % table by country
    tblCountry = table();
    tblCountry.Country = GrCountry;
    tblCountry.ValuebyCountry = splitapply(@sum, valuetbySector.Accvaluefee,...
        valuetbySector.CountryID);
    tblCountry.PercbyCountry = 100*(tblCountry.ValuebyCountry ./ ...
        sum(tblCountry.ValuebyCountry));
    tblCountry.ValuetbyCountry = splitapply(@sum, valuetbySector.Accvaluetwdiv,...
        valuetbySector.CountryID);
    tblCountry.RetbyCountry = (tblCountry.ValuetbyCountry ./...
        tblCountry.ValuebyCountry) - 1;
    
    % D) Risk Table
    riskTablet = table();
    % get Risk Metrics (VAR, ES) for each ticker
    names2filteridxES = startsWith(AllRiskHistoricaltable.Properties.VariableNames,...
        'ES_');
    names2filteridxVaR = startsWith(AllRiskHistoricaltable.Properties.VariableNames,...
        'VAR_');
    names2filteridxVARCF = startsWith(AllRiskHistoricaltable.Properties.VariableNames,...
        'VaRCF_');
    names2filteridxESCF = startsWith(AllRiskHistoricaltable.Properties.VariableNames,...
        'ESCF_');
    % Get names
    names2filterES = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxES);
    names2filterES(strcmp( names2filterES, 'ES_USD_EUR_USDollarEuro')) = [];
    names2filterVAR = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxVaR);
    names2filterVAR(strcmp(names2filterVAR, 'VAR_USD_EUR_USDollarEuro')) = [];
    names2filterVARCF = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxVARCF);
    names2filterVARCF(strcmp(names2filterVARCF, 'VaRCF_USD_EUR_USDollarEuro')) = [];
    names2filterESCF = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxESCF);
    names2filterESCF(strcmp(names2filterESCF, 'ESCF_USD_EUR_USDollarEuro')) = [];
    % table
    riskTablet.Name = tableprices0t.UniqueSymbol;
    riskTablet.ExpShortFall = AllRiskHistoricaltable{end, names2filterES}';
    riskTablet.VaR = AllRiskHistoricaltable{end, names2filterVAR}';
    riskTablet.VaRCF = AllRiskHistoricaltable{end, names2filterVARCF}';
    riskTablet.ExpShortFallCF = AllRiskHistoricaltable{end, names2filterESCF}';
    % A VAR-ES sum up table
    sumUpVaRTable = table();
    sumUpESTable = table();
    % Read Backtesting files
    for h = 1: length(namestocharge)
        for j = 1:2
            ticker = namestocharge(h);
            path2save = fullfile(pathfilewm, ticker);
            tableVar = readtable(char(strcat(path2save, '\',...
                'Backtesting_', ticker, '.xlsx')), 'filetype', 'spreadsheet', ...
                'sheet', j);
            sumUpVaRTable = [sumUpVaRTable; tableVar];
        end
        for k = 3:4
            tableES = readtable(char(strcat(path2save, '\',...
                'Backtesting_', ticker, '.xlsx')), 'filetype', 'spreadsheet', ...
                'sheet', k);
            sumUpESTable = [sumUpESTable; tableES];
        end
    end
    
    % A worst day and average values sum up table
    sumUpWorstDayTable = table();
    sumUpAverageValuesTable = table();
    for h = 1: length(namestocharge)
        ticker = namestocharge(h);
        path2save = fullfile(pathfilewm, ticker);
        % worst day table
        tableWorstDay = readtable(char(strcat(path2save, '\',...
            ticker, '_',datestr(today()),'.xlsx')), 'filetype', 'spreadsheet', ...
            'sheet', 2);
        variableNames = {'ID_','Hist_', 'Ret_', 'VAR_', 'ES_', 'VaRCF_', ...
            'ESCF_'};
        variableNamestTick = ['Time', strcat(variableNames, char(ticker))];
        tableWorstDay = tableWorstDay(:,  variableNamestTick);
        tableWorstDay.Properties.VariableNames = {'Worst_day','ID',...
            'HistPrice', 'Ret', 'VAR', 'ES', 'VaRCF','ESCF'};
        sumUpWorstDayTable = [sumUpWorstDayTable; tableWorstDay];
        % average values table
        tableAverageValues = readtable(char(strcat(path2save, '\',...
            ticker,'_',datestr(today()),'.xlsx')), 'filetype', 'spreadsheet', ...
            'sheet', 3);
        sumUpAverageValuesTable = [sumUpAverageValuesTable; tableAverageValues];
    end
    
    %% History - Portfolio Analysis
    disp(['Portfolio Analysis at ', char(datetime('now','format','HH:mm'))]);
    % All portfolio analysis pnl(t) - all portfolio money earned
    names2filteridxpnl0t = startsWith(AllRiskHistoricaltable.Properties.VariableNames,...
        'PNL0T_');
    names2filterpnl0t = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxpnl0t);
    names2filterpnl0t(strcmp(names2filterpnl0t, 'PNL0T_USD_EUR_USDollarEuro')) = [];
    Pnl0tallData = AllRiskHistoricaltable{:, names2filterpnl0t};
    sumPnl0t = sum(Pnl0tallData,2, 'omitnan');
    % pnl(0) - all portfolio money invested
    names2filteridxpnl0 = ...
        startsWith(AllRiskHistoricaltable.Properties.VariableNames,'PNL0_');
    names2filterpnl0 = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxpnl0);
    names2filterpnl0(strcmp(names2filterpnl0, 'PNL0_USD_EUR_USDollarEuro')) = [];
    Pnl0allData = AllRiskHistoricaltable{:, names2filterpnl0};
    sumPnl0 = sum(Pnl0allData, 2);
    % pnl - all portfolio money
    names2filteridxpnl = ...
        startsWith(AllRiskHistoricaltable.Properties.VariableNames, 'PNL_');
    names2filterpnl = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxpnl);
    names2filterpnl(strcmp(names2filterpnl, 'PNL_USD_EUR_USDollarEuro')) = [];
    PnlallData = AllRiskHistoricaltable{:, names2filterpnl};
    sumPnl = sum(PnlallData, 2);
    % historical percentage earned by positions
    histweigths = PnlallData(:,:) ./ sumPnl ;
    names2filteridxretacc = startsWith(AllRiskHistoricaltable.Properties.VariableNames,...
        'RetAcc_');
    names2filterRetAcc = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxretacc);
    RetAccallData = AllRiskHistoricaltable(:, names2filterRetAcc);
    RetAccallData(:,'RetAcc_USD_EUR_USDollarEuro') = [];
    pnlpercbypos = sum(histweigths .* table2array(RetAccallData),2);
    % VAR / ES with historical weigths no diversification
    hExpShortFall = AllRiskHistoricaltable{:, names2filterES};
    hExpShortFallbypos = -sqrt(sum((histweigths).^2 .* (hExpShortFall).^2, 2));
    hVaR = AllRiskHistoricaltable{:, names2filterVAR};
    hVaRbypos = -sqrt(sum((histweigths).^2 .* (hVaR).^2, 2));
    hVaRCF = AllRiskHistoricaltable{:, names2filterVARCF};
    hVaRCFbypos = -sqrt(sum((histweigths).^2 .* (hVaRCF).^2, 2));
    hExpShortFallCF = AllRiskHistoricaltable{:, names2filterESCF};
    hExpShortFallCFbypos = -sqrt(sum((histweigths).^2 .* (hExpShortFallCF).^2, 2));
    % backtesting
    % VAR
    [hVARaggsumUp, hVARaggTLtable, hVARaggPofTable,...
        hVARaggCCITable, hVARaggTBFITable] =...
        backtestingVaRdetail({'hVaRbypos' ; 'hVaRCFbypos'}, [pnlpercbypos,pnlpercbypos] , ...
        [hVaRbypos, hVaRCFbypos], {'t', 'cf'}, 0.95);
    % ES
    % to backtest Portfolio ES modelling the portfolio is needed
    pricesPortf = sumPnl./sumPnl0;
    path2saveportf = fullfile(pathfilewm,'Portfolio');
    ESportf = calcParamES(price2ret(pricesPortf),...
        'CFn', 0.05, path2saveportf, 5);
    tableBacktestinghESaggt = ...
        bcktstingexepshrtfall2pnl(pnlpercbypos, pricesPortf,...
        [hVaRbypos, hVaRCFbypos],[hExpShortFallbypos, hExpShortFallCFbypos],...
        {'t', 'cf'}, 0.05, path2saveportf);
    ColumNamesIdentifiershES = {'hExpShortFallbypos';'hExpShortFallbypos';...
        'hExpShortFallCFbypos'};
    tableBacktestinghESaggWithId = table();
    tableBacktestinghESaggWithId.Identifier = ColumNamesIdentifiershES;
    tableBacktestinghESaggWithId = [tableBacktestinghESaggWithId, ...
        tableBacktestinghESaggt];
    tableBacktestinghESaggWithId = rmmissing(tableBacktestinghESaggWithId);
    % percentage earned
    pnlperc0 = pricesPortf - 1; %incorrect: sumPnl0t ./ sumPnl0;
    pnlvald = diff(sumPnl0t);
    pnlperc0 = fillmissing(pnlperc0, 'constant', 0);
    % syntetic benchmark and performance analysis
    idx = ismember(tbl.UniqueSymbol, tableprices0t.UniqueSymbol);
    ticksold = tbl.UniqueSymbol(~idx);
    idxt = ismember(dailyinfoportf.symbol, ticksold);
    hPositions = dailyinfoportf(~idxt, :);
    tablewithBench = tableprices0t(:,{'UniqueSymbol','PercByTickT'});
    benchUnique = unique(hPositions(:,{'symbol','BenchSymbol'}));
    tablewithBench.Benchmark = benchUnique.BenchSymbol; 
    [GroupsBench, Bench] = findgroups(tablewithBench.Benchmark);
    investedbyExch = splitapply(@sum, tablewithBench.PercByTickT, GroupsBench);
    histBench = AllHistoricalTable(:, matlab.lang.makeValidName(Bench));
    newBench = [0; cumsum(price2ret(histBench{:, :} * investedbyExch))];
    % daily money earned by category
    % category defined in GrAssetType
    VpnlByCategory = timetable(AllHistoricalTable.Date);
    pnlPercByCategory = timetable(AllHistoricalTable.Date);
    vpnlDataInvest = timetable(AllHistoricalTable.Date);
    vpnlDataMV = timetable(AllHistoricalTable.Date);
    for k = 1:numel(GrAssetType)
        % filter asset type
        idxAt = ismember(hPositions.Category1, GrAssetType(k));
        % get unique tickers
        tickers = matlab.lang.makeValidName(unique(hPositions{idxAt, 'symbol'}));
        % add label 
        tickersFilterHT = strcat('Vpnl_', tickers);
        tickersFilterIM = strcat('PNL0_', tickers);
        tickersFilterMV = strcat('PNL_', tickers);
        % get vnpl data
        vpnlData = AllRiskHistoricaltable(:, tickersFilterHT);
        % get invested money and MV
        vpnlDataInvestIn = AllRiskHistoricaltable(:, tickersFilterIM);
        vpnlDataMVIn = AllRiskHistoricaltable(:, tickersFilterMV);
        % sum vpnl in each category
        vpnlBycat = sum(vpnlData{:,:},2);
        vpnlInvestBycat = sum(vpnlDataInvestIn{:,:},2);
        tvpnlInvestBycat = table(vpnlInvestBycat);
        vpnlBycat = table(vpnlBycat - [0;diff(vpnlInvestBycat)]);
        vpnlMVBycat = sum(vpnlDataMVIn{:,:},2);
        tvpnlMVBycat = table(vpnlMVBycat);
        percBycat = table(vpnlMVBycat ./ vpnlInvestBycat -1);
        % names
        vpnlBycat.Properties.VariableNames = {strcat('vpnl_', char(GrAssetType(k)))}; 
        percBycat.Properties.VariableNames = {strcat('Percpnl_', char(GrAssetType(k)))};
        tvpnlMVBycat.Properties.VariableNames = {strcat('MV_', char(GrAssetType(k)))}; 
        tvpnlInvestBycat.Properties.VariableNames = {strcat('IM_', char(GrAssetType(k)))};
        % concat tables
        VpnlByCategory = [VpnlByCategory, vpnlBycat];
        pnlPercByCategory  = [pnlPercByCategory, percBycat];
        vpnlDataInvest = [vpnlDataInvest, tvpnlInvestBycat];
        vpnlDataMV = [vpnlDataMV, tvpnlMVBycat];
    end 
    % fillmissing
    VpnlByCategory = fillmissing(VpnlByCategory, 'constant', 0);
    pnlPercByCategory = fillmissing(pnlPercByCategory, 'constant', 0);
    % Percentage
    VpnlByCategoryPerc = array2table(VpnlByCategory{:,:} ./ sum(VpnlByCategory{:,:},2));
    VpnlByCategoryPerc = fillmissing(VpnlByCategoryPerc, 'constant', 0);
    VpnlByCategoryPerc.Properties.VariableNames = strcat('perc_', VpnlByCategory.Properties.VariableNames);
    % Returns in no-crypto
    VpnlPercInNoCrypto = sum(vpnlDataMV{:,{'MV_Stock', 'MV_ETF', 'MV_Fund'}}, 2) ./ ...
        sum(vpnlDataInvest{:,{'IM_Stock', 'IM_ETF', 'IM_Fund'}}, 2) - 1;
    pnlPercByCategory.PercInNoCrypto = VpnlPercInNoCrypto; 
    % all in one 
    gnPnlByCategory = [VpnlByCategory, VpnlByCategoryPerc, pnlPercByCategory];
    
    % performance analysis
    % Jensen' alpha and Information Ratio
    alphavec = pnlperc0 - newBench;
    period = 14;
    for k = 1:length(pnlperc0) - period
        [inforatioc(k:k+period), TrackingError(k:k+period)] = ....
            inforatio(pnlperc0(k:k+period), newBench(k:k+period));
    end
    % VaR /ES portfolio
    % weigth for each ticker
    [GroupsTickers,GN] = findgroups(dailyinfoportf.symbol); %
    tblb = table();
    tblb.Name = GN;
    tblb.investedbytick = splitapply(@sum, dailyinfoportf.quantity,...
        GroupsTickers);
    tblb.holdP = splitapply(@sum, dailyinfoportf.shares,...
        GroupsTickers);
    tblb(tblb.holdP <= 0.01,:) = [];
    tblb = innerjoin(tblb, tablepricesMiddlePrice);
    wtick = tblb.investedbytick ./ sum(hPositions.quantity);
    retallData = AllHistoricalRetTable(:, matlab.lang.makeValidName(tblb.Name)');
    % we need returns and the evolution of the portfolio
    %     retallData = timetable2table(retallData);
    if length(sumPnl) > height(retallData)
        patch = retallData(1, :);
        patch.Time = patch.Time -1;
        retallData = [patch; retallData];
    end
    % GARCH and copulae
    [portfVaRndESinPerc, portfVaRndESinMoney, rhoT, rhoG] = ...
        CalculateVaRandESdinamicwithCopula(retallData, sumPnl,...
        wtick, [0;pnlvald]./sumPnl);
    % VaR/ES in Money
    hExpShortFallbyposval = hExpShortFallbypos .* sumPnl;
    hVaRbyposval =  hVaRbypos .* sumPnl;
    hVaRCFbyposval = hVaRCFbypos .* sumPnl;
    hExpShortFallCFbyposval = hExpShortFallCFbypos .* sumPnl;
    portfvarval =  portfVaRndESinMoney(:, 1);
    portfvarwrval = portfVaRndESinMoney(:, 2);
    portfvarwgval = portfVaRndESinMoney(:, 3);
    portfesval = portfVaRndESinMoney(:, 4);
    portfeswrval = portfVaRndESinMoney(:, 5);
    portfeswgval = portfVaRndESinMoney(:, 6);
    % backtesting VAR
    pnlperc2backtest = [0;pnlvald] ./ sumPnl;
    [VARaggsumUp, VARaggTLtable, VARaggPofTable,...
        VARaggCCITable, VARaggTBFITable] =...
        backtestingVaRdetail({'PortfTVaR'; 'PortfT_copVaR'; 'PortfG_copVaR'}, ...
        [pnlperc2backtest, pnlperc2backtest, pnlperc2backtest] , ...
        portfVaRndESinPerc(:,1:3), ...
        {'t', 't', 'normal'}, 0.95);
    % backtesting ES
    tableBacktestingESagg = bcktstingexepshrtfall2pnl(pnlperc2backtest,pricesPortf,...
        portfVaRndESinPerc(:,1:3), portfVaRndESinPerc(:,4:6),...
        {'t', 't', 'normal'}, 0.05, path2saveportf);
    ColumNamesIdentifiersES = {'PortfT_ES'; 'PortfT_ES'; ...
        'PortfT_copES'; 'PortfT_copES'; 'PortfG_copES'; 'PortfG_copES'};
    tableBacktestingESaggWithId = table();
    tableBacktestingESaggWithId.Identifier = ColumNamesIdentifiersES;
    tableBacktestingESaggWithId = [tableBacktestingESaggWithId, ...
        tableBacktestingESagg];
    tableBacktestingESaggWithId = rmmissing(tableBacktestingESaggWithId);
    % sharpe modified (with VAR)
    sharpemod = alphavec ./ abs(hVaRbypos) ;
    % treynor ratio
    portfolioBeta = calculateBeta(pnlperc0, newBench, 0);
    treynorRatio = alphavec ./ portfolioBeta;
    % table
    RiskPortfolioAnalysis = timetable(AllHistoricalTable.Date, ...
        sumPnl-sumPnl0, [0;diff(sumPnl-sumPnl0)],...
        [0;diff(pnlperc0)],... %(sumPnl-sumPnl0)./sumPnl , ... %incorrect
        movstd([0;pnlvald]./sumPnl, 4), mov2std([0;pnlvald]./sumPnl), ...
        sumPnl0, sumPnl, pnlperc0, newBench, pnlpercbypos,...
        inforatioc', hVaRbyposval, hExpShortFallbyposval, portfvarwrval, ...
        portfeswrval, treynorRatio, sharpemod);
    columnamesRiskAnalysis = {'vPNL0','ValPnl', 'ValPnlperc' ,'vol_4d',...
        'volAcc_4d', 'PNL0', 'PNL', 'PNLperc', 'nBENC', 'PNLpercbyPos', ...
        'IR', 'hparmtVaR', 'hparmtES', 'CopVaRtVal','copEStVal', ...
        'SharpeRatioMod', 'TreynorRatio'};
    RiskPortfolioAnalysis.Properties.VariableNames = columnamesRiskAnalysis;
    % VAR and ES table
    VaRandESPortfolioAnalysis = timetable(AllHistoricalTable.Date,...
        hExpShortFallbyposval, hVaRbyposval, hVaRCFbyposval,...
        hExpShortFallCFbyposval, portfvarval, portfvarwrval, portfvarwgval,...
        portfesval, portfeswrval, portfeswgval);
    columnamesVaRnESAnalysis = {'hESbyPos', 'hVaRbyPos', 'hVaRCFbyPos', ...
        'hESCFbyPos','PortfTVaR', 'PortfT_copVaR', 'PortfG_copVaR',...
        'PortfT_ES', 'PortfT_copES', 'PortfG_copES' };
    VaRandESPortfolioAnalysis.Properties.VariableNames = columnamesVaRnESAnalysis;
    % calculate VAR / ES with copulae with simulated returns
    % Variance Gamma
    VGpricetableToCopula = table();
    VGpricetableToCopula.Time = VGPriceTable.Time(2:end);
    VGPriceTablef = VGPriceTable{:,strcat('VGPrices_', ...
        matlab.lang.makeValidName(tblb.Name))};
    VGPriceTablef = array2table(price2ret(VGPriceTablef));
    VGPriceTablef.Properties.VariableNames = strcat('VGPrices_', ...
        matlab.lang.makeValidName(tblb.Name));
    %     VGPriceTablef.VGPrices_USD_EUR_USDollarEuro = [];
    VGpricetableToCopula =  [VGpricetableToCopula , VGPriceTablef];
    [portfSimulVGVaRndESinPerc, portfSimulVGVaRndESinMoney, rhoTsVG, rhoGsVG] = ...
        CalculateVaRandESdinamicwithCopula(table2timetable(VGpricetableToCopula),...
        sumPnl,wtick, [0;pnlvald]./sumPnl);
    VGtableVARndES = timetable(AllHistoricalTable.Date,...
        portfSimulVGVaRndESinMoney(:,1), portfSimulVGVaRndESinMoney(:,2), ...
        portfSimulVGVaRndESinMoney(:,3), portfSimulVGVaRndESinMoney(:,4),...
        portfSimulVGVaRndESinMoney(:,5), portfSimulVGVaRndESinMoney(:,6)) ;
    VGtableVARndES.Properties.VariableNames = {'PortfTVaRwVG', ...
        'PortfT_copVaRwVG', 'PortfG_copVaRwVG',...
        'PortfT_ESwVG', 'PortfT_copESwVG', 'PortfG_copESwVG'};
    %     backtesting
    %     VAR
    [VARaggwVGsumUp, VARaggwVGTLtable, VARaggwVGPofTable,...
        VARaggwVGCCITable, VARaggwVGTBFITable] =...
        backtestingVaRdetail({'PortfTVaRwVG'; 'PortfT_copVaRwVG'; ...
        'PortfG_copVaRwVG'}, ...
        [pnlperc2backtest, pnlperc2backtest, pnlperc2backtest] , ...
        portfSimulVGVaRndESinPerc(:,1:3), ...
        {'vg', 'vg', 'vg'}, 0.95);
    % ES
    tableBacktestingESsimulVG = bcktstingexepshrtfall2pnl(pnlperc2backtest,...
        pricesPortf, portfSimulVGVaRndESinPerc(:,1:3), ...
        portfSimulVGVaRndESinPerc(:,4:6), {'vg', 'vg', 'vg'}, 0.05,...
        path2saveportf);
    ColumNamesIdentifiersESwVG = {'PortfT_ESwVG'; 'PortfT_copESwVG'; ...
        'PortfG_copESwVG'};
    tableBacktestingESwVGaggWithId = table();
    tableBacktestingESwVGaggWithId.Identifier =  ColumNamesIdentifiersESwVG ;
    tableBacktestingESwVGaggWithId = [tableBacktestingESwVGaggWithId, ...
        tableBacktestingESsimulVG];
    tableBacktestingESwVGaggWithId = rmmissing(tableBacktestingESwVGaggWithId);
    %t distributed prices
    tpricetableToCopula = table();
    tpricetableToCopula.Time = tDistributedTable.Time(2:end);
    trTablef = tDistributedTable{:,strcat('tPrices_', ...
        matlab.lang.makeValidName(tblb.Name))};
    trTablef = array2table(price2ret(trTablef));
    trTablef.Properties.VariableNames = strcat('tPrices_',...
        matlab.lang.makeValidName(tblb.Name));
    %     trTablef.tPrices_USD_EUR_USDollarEuro = [];
    tpricetableToCopula =  [tpricetableToCopula , trTablef];
    [portfSimultVaRndESinPerc, portfSimulTVaRndESinMoney, rhoTt, rhoGt] = ...
        CalculateVaRandESdinamicwithCopula(table2timetable(tpricetableToCopula),...
        sumPnl, wtick, [0;pnlvald]./sumPnl);
    TtableVARndES = timetable(AllHistoricalTable.Date,...
        portfSimulTVaRndESinMoney(:,1), portfSimulTVaRndESinMoney(:,2), ...
        portfSimulTVaRndESinMoney(:,3), portfSimulTVaRndESinMoney(:,4),...
        portfSimulTVaRndESinMoney(:,5), portfSimulTVaRndESinMoney(:,6)) ;
    TtableVARndES.Properties.VariableNames = {'PortfTVaRwT', ...
        'PortfT_copVaRwT', 'PortfG_copVaRwT',...
        'PortfT_ESwT', 'PortfT_copESwT', 'PortfG_copESwT'};
    simulatedMetricsTable = [VGtableVARndES,TtableVARndES];
    % backtesting
    % VAR
    [VARaggwTsumUp, VARaggwTTLtable, VARaggwTPofTable,...
        VARaggwTCCITable, VARaggwTTBFITable] =...
        backtestingVaRdetail({'PortfTVaRwT'; ...
        'PortfT_copVaRwT'; 'PortfG_copVaRwT'}, ...
        [pnlperc2backtest, pnlperc2backtest, pnlperc2backtest] , ...
        portfSimultVaRndESinPerc(:,1:3), ...
        {'t', 't', 't'}, 0.95);
    % ES
    tableBacktestingESsimulT = bcktstingexepshrtfall2pnl(pnlperc2backtest,...
        pricesPortf, portfSimultVaRndESinPerc(:,1:3),...
        portfSimultVaRndESinPerc(:,4:6), {'t', 't', 't'}, 0.05, path2saveportf);
    ColumNamesIdentifiersESwT = {'PortfT_ESwT'; 'PortfT_ESwT'; ...
        'PortfT_copESwT'; 'PortfT_copESwT'; 'PortfG_copESwT';'PortfG_copESwT'};
    tableBacktestingESwTaggWithId = table();
    tableBacktestingESwTaggWithId.Identifier =  ColumNamesIdentifiersESwT;
    tableBacktestingESwTaggWithId = [tableBacktestingESwTaggWithId, ...
        tableBacktestingESsimulT];
    tableBacktestingESwTaggWithId = rmmissing(tableBacktestingESwTaggWithId);
    
    %% monthly - portfolio // aggr. portfolio metrics
    disp(['Monthly - portfolio analysis at ', char(datetime('now','format','HH:mm'))])
    % obtain a historical pnl
    names2filteridxpnl = startsWith(AllRiskHistoricaltable.Properties.VariableNames,...
        'PNL_');
    names2filterpnl = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxpnl);
    names2filterpnl(strcmp(names2filterpnl0, 'PNL_USD_EUR_USDollarEuro')) = [];
    PnlallData = AllRiskHistoricaltable{:, names2filterpnl};
    sumPnl = timetable(AllRiskHistoricaltable.Time,sum(PnlallData,2));
    % quantity invested every month
    tickhold = tbl.UniqueSymbol(idx);
    idxh = ismember(dailyinfoportf.symbol, tickhold);
    dailyinfoportfh = dailyinfoportf(idxh, :);
    monthlyQuantity = groupsummary(table2timetable(dailyinfoportfh),...
        'Date','month', 'sum', 'loc_quantity');
    monthlyQuantity.AccumQuantity = cumsum(monthlyQuantity.sum_loc_quantity);
    % calculate pnl end of month
    m = months(datestr(sumPnl.Time(1)), datestr(sumPnl.Time(end)));
    ti = dateshift(sumPnl.Time(1) + calmonths(0:m),'end','month');
    if datestr(today()) < dateshift(datetime(datestr(today())), 'end', 'month')
        tdate = [ti,datetime(datestr(today()))];
    else
        tdate = ti;
    end
    monthlyPnl = timetable2table(sumPnl(tdate,:));
    % aggregate and calculate monthly return
    monthlyQuantity.eom_PNL = monthlyPnl{:,2};
    % pnlVal
    monthlyQuantity.eom_PNLval =  monthlyQuantity.eom_PNL - ...
        monthlyQuantity.AccumQuantity ;
    % return
    monthlyQuantity.eom_return = monthlyQuantity.eom_PNL ./...
        monthlyQuantity.AccumQuantity - 1;
    newBenchtt = timetable(AllRiskHistoricaltable.Time, newBench);
    monthlyBench = timetable2table(newBenchtt(tdate,:));
    monthlyQuantity.eom_BenchReturn = monthlyBench{:, 2};
    % MWR  as IRR
    cashflow2irr = [-table2array(monthlyQuantity(:, 3));...
        monthlyQuantity.eom_PNL(end)];
    % MWR with fees and dividends
    monthlyDividend = groupsummary(table2timetable(dividends), 'Date',...
        'month', 'sum', 'value');
    % fees
    monthlyQuantitywfee = groupsummary(table2timetable(dailyinfoportfh),...
        'Date','month', 'sum', 'fee');
    info2othfees = readtable(xlsfileaccpath, 'filetype', 'spreadsheet',...
        'sheet', 5);
    info2othfees.Date = datetime(info2othfees.Date, 'InputFormat',...
        'dd-MM-yyyy');
    monthlyQuantitywotherfee = groupsummary(info2othfees,...
        'Date','month', 'sum', 'value');
    monthlyQuantitywfee.Quantity = monthlyQuantity.sum_loc_quantity;
    monthlyQuantitywfee.Quantitywfee = monthlyQuantitywfee.Quantity +...
        monthlyQuantitywfee.sum_fee - monthlyQuantitywotherfee.sum_value;
    monthlyQuantitywfee.AccumQuantitywfee = cumsum(monthlyQuantitywfee.Quantitywfee);
    % pnl and dividends
    monthlyQuantitywfee.eom_PNL = monthlyQuantity.eom_PNL;
    monthlyQuantitywfee.totDiv = monthlyDividend.sum_value;
    monthlyQuantitywfee.PNLwDiv = monthlyQuantitywfee.eom_PNL + ...
        monthlyQuantitywfee.totDiv;
    monthlyQuantitywfee.PNLwDivVal = monthlyQuantitywfee.PNLwDiv -  ...
        monthlyQuantitywfee.AccumQuantitywfee ;
    % return
    monthlyQuantitywfee.eom_return = monthlyQuantitywfee.PNLwDiv ./...
        monthlyQuantitywfee.AccumQuantitywfee - 1;
    monthlyQuantitywfee.eom_BenchReturn = monthlyQuantity.eom_BenchReturn;
    % aggregate portfolio analysis
    % total div
    totalDiv = sum(dividensTable.AccValDiv);
    % MWR  as IRR with div and fee
    cashflow2irrwfeendiv = [-table2array(monthlyQuantitywfee(:, 5));...
        totalDiv; monthlyQuantitywfee.PNLwDiv(end)];
    % portfolio TWR
    prtfTWRwfeendiv = calculatePortfTWR(AllRiskHistoricaltable,sumPnl0,...
        dailyinfoportfh, dividends, 'divnfee');
    % aggr. portfolio metrics
    aggrPortfTable = table();
    % returns
    % MWR
    aggrPortfTable.MWRportf = irr(cashflow2irr);
    aggrPortfTable.MWRportfwfeendiv = irr(cashflow2irrwfeendiv);
    aggrPortfTable.MWRportfbytick = sum(tablepricesMiddlePrice.PercByTick .*...
        tablepricesMiddlePrice.MWR);
    aggrPortfTable.MWRwdivnfeeportfbytick = sum(tablepricesMiddlePrice.PercByTick .*...
        tablepricesDivndFee.MWRtwd);
    % TWR
    aggrPortfTable.TWRportfwdivnfee = prtfTWRwfeendiv;
    aggrPortfTable.TWRportfbytick = sum(tablepricesMiddlePrice.PercByTick .*...
        tablepricesMiddlePrice.TWR);
    aggrPortfTable.TWRwdivnfeeportfbytick = sum(tablepricesMiddlePrice.PercByTick .*...
        tablepricesDivndFee.TWRtwd);
    aggrPortfTable.portfPercBytick = sum(tableprices0t.PercByTickt .*...
        tableprices0t.PercUsingPrice);
    aggrPortfTable.portfPercBytickMidPrice = sum(tablepricesMiddlePrice.PercByTick .* ...
        tablepricesMiddlePrice.PercByMidPr);
    % beta
    aggrPortfTable.Beta_w = sum(tablepricesMiddlePrice.PercByTick .* ...
        tableprices0t.Beta);
    aggrPortfTable.Betap = portfolioBeta(end);
    % Sharpe, Sortino, Informatio ratio,
    aggrPortfTable.alphavec = portalpha(pnlperc0,newBench,0,'xs');
    aggrPortfTable.SharpeMod = sharpemod(end);
    aggrPortfTable.Treynor = treynorRatio(end);
    aggrPortfTable.PortfReturn = pnlperc0(end);
    aggrPortfTable.BenchReturn = newBench(end);
    aggrPortfTable.TrackingError = TrackingError(end);
    aggrPortfTable.InfoRatio = inforatioc(end);
    % risk metrics
    aggrPortfTable.PortfES = -sqrt(sum((tableprices0t.PercByTickt(1:end-1)).^2 .*...
        (riskTablet.ExpShortFall(1:end-1)).^2));
    aggrPortfTable.PortfVaR = -sqrt(sum((tableprices0t.PercByTickt(1:end-1)).^2 .*...
        (riskTablet.VaR(1:end-1)).^2));
    aggrPortfTable.PortfVaRCF = -sqrt(sum((tableprices0t.PercByTickt(1:end-1)).^2 .* ...
        (riskTablet.VaRCF(1:end-1)).^2));
    aggrPortfTable.PortfESCF = -sqrt(sum((tableprices0t.PercByTickt(1:end-1)).^2 .* ...
        (riskTablet.ExpShortFallCF(1:end-1)).^2));
    aggrPortfTable2w = rows2vars(aggrPortfTable);
    aggrPortfTable2w.Properties.VariableNames = {'Metric', 'Value'};

    %% plot
    disp(['Is plot time!! at ', char(datetime('now','format','HH:mm'))])
    % plot  exposure pie
    g = figure('visible','off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
    subplot(1,2,1)
    pie(tableprices0t.PercByTickt(1:end))
    title('Percentage invested by position')
    legend(tableprices0t.UniqueSymbol(1:end),...
        'location', 'northOutside', 'Orientation','horizontal', ...
        'NumColumns',2)
    subplot(1,2,2)
    names = categorical(tablepricesMiddlePrice.Name(1:end));
    bar(names, [tablepricesMiddlePrice.PercByMidPr(1:end),...
        tablepricesDivndFee.TWRtwd(1:end), tablepricesDivndFee.MWRtwd(1:end),...
        riskTablet.ExpShortFall(1:end)])
    title('Main Stats by position')
    legend('PrByMidPr', 'TWR', 'MWR', 'ExpShortfall', 'location', 'south')
    ylabel('%'), xlabel('Position')
    saveas(g, char(strcat(pathfilewm, '\',...
        'Main_Stats(1)')), 'jpg')
    h = figure('visible','off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
    subplot(1,2,1)
    pie(tableprices0t.PercByTickT(1:end))
    title('Exposure percentage by position')
    legend(tableprices0t.UniqueSymbol(1:end),...
        'location', 'northOutside', 'Orientation','horizontal', ...
        'NumColumns',2)
    subplot(1,2,2)
    bar(monthlyQuantitywfee.month_Date,...
        [monthlyQuantitywfee.eom_return, monthlyQuantitywfee.eom_BenchReturn])
    title('Monthly Returns')
    legend('Portfolio','Benchmark', 'location', 'south')
    ylabel('%'), xlabel('Month')
    saveas(h, char(strcat(pathfilewm, '\',...
        'Main_Stats(2)')), 'jpg')
    % plot history portfolio analysis
    f = figure('visible','off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
    % plot aggregate returns and benchmark, legend, title...
    subplot(1,3,1)
    plotTimeSeries(RiskHistoricaltable.Time, [pnlperc0, newBench,...
        pnlpercbypos], ['PNL and Benchmark in ' ,...
        datestr(RiskHistoricaltable.Time(1), 'dd-mm-yyyy'), ' to ',...
        datestr(RiskHistoricaltable.Time(end), 'dd-mm-yyyy')],...
        {'Return','Benchmark', 'ReturnbyPosition'}, 'Date', '%')
    subplot(1,3,2)
    plotTimeSeries(RiskHistoricaltable.Time, [[0;diff(sumPnl{:,:}-sumPnl0)], ...
        portfvarval, portfvarwrval, portfvarwgval, ...
        portfesval, portfeswrval, portfeswgval],'Portfolio VaR and ES',...
        {'PNL', 't-VAR', 'tCopula-VAR', 'gCopula-VAR', 't-ES',...
        'tCopula-ES', 'gCopula-ES'},'Date', 'Eur' )
    subplot(1,3,3)
    plotTimeSeries(RiskHistoricaltable.Time,[[0;diff(sumPnl{:,:}-sumPnl0)],...
        hExpShortFallbyposval, hVaRbyposval, hVaRCFbyposval, ...
        hExpShortFallCFbyposval],...
        'Historical PortVaR and PortES',{'PNL',...
        'ES', 'VAR', 'VARcf', 'EScf'}, 'Date', 'Eur')
    saveas(f, char(strcat(pathfilewm, '\', 'PNLandBenchmark')), 'jpg')
    % plot returns distributions
    ff = figure('visible','off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
    % portfolio returns
    subplot(1,2,1)
    histogram(diff(pnlperc0), 40);
    title(['Returns from portfolio in ' ,...
        datestr(RiskHistoricaltable.Time(1), 'dd-mm-yyyy'), ' to ',...
        datestr(RiskHistoricaltable.Time(end), 'dd-mm-yyyy')]);
    xlabel('value'), ylabel('frecuency')
    % benchmark returns
    subplot(1,2,2)
    histogram(diff(newBench), 40);
    title(['Returns from Benchmark in ' ,...
        datestr(RiskHistoricaltable.Time(1), 'dd-mm-yyyy'), ' to ',...
        datestr(RiskHistoricaltable.Time(end), 'dd-mm-yyyy')]);
    xlabel('value'), ylabel('frecuency')
    saveas(ff, char(strcat(pathfilewm, '\', ...
        'PNLandBenchmarkDistribution')), 'jpg')
    % plot returns for each ticker
    g = figure('visible','off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
    % obtain a historical returns for each ticker
    names2filteridxret = startsWith(AllRiskHistoricaltable.Properties.VariableNames,...
        'VPNL0T_');
    names2filterret = AllRiskHistoricaltable.Properties.VariableNames(names2filteridxret);
    retallDataPort = AllRiskHistoricaltable{:, names2filterret};
    % plot, legend, title...
    plotTimeSeries(RiskHistoricaltable.Time, [retallDataPort, pnlperc0], ...
        'Returns By Middle price and Portfolio',...
        [tableprices0t.UniqueSymbol;'Portfolio'], 'Date', '%')
    saveas(g, char(strcat(pathfilewm, '\',...
        'IndividualReturnsByTicker')), 'jpg')
    % plot Returns with simulated VAR and ES with VG a T distribution
    hh = figure('visible','off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
    subplot(1,2,1)
    plotTimeSeries(RiskHistoricaltable.Time, [pnlperc2backtest,...
        portfSimulVGVaRndESinPerc], 'Daily returns, VAR and ES with VG model',...
        {'Daily pnl', 'PortfTVaRwVG', 'PortfTcopVaRwVG', 'PortfGcopVaRwVG',...
        'PortfTESwVG', 'PortfTcopESwVG', 'PortfGcopESwVG'}, 'Date', '%')
    subplot(1,2,2)
    plotTimeSeries(RiskHistoricaltable.Time,[pnlperc2backtest,...
        portfSimultVaRndESinPerc], 'Daily returns, VAR and ES with T model',...
        {'Daily pnl','PortfTVaRwT', 'PortfTcopVaRwT', 'PortfGcopVaRwT',...
        'PortfTESwT', 'PortfTcopESwT', 'PortfGcopESwT'}, 'Dates', '%')
    saveas(hh, char(strcat(pathfilewm, '\',...
        'VaRndESsimulatedModels')), 'jpg')
    % plot returns by category asset
    hg = figure('visible','off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
    subplot(3,3,1)
        plotTimeSeries(RiskHistoricaltable.Time, pnlPercByCategory{:,'Percpnl_Stock'}, ...
        'Returns In Stock', {'Return'}, 'Date', '%')
    subplot(3,3,4)
    plotTimeSeries(RiskHistoricaltable.Time, pnlPercByCategory{:,'Percpnl_ETF'}, ...
        'Returns In ETF', {'Return'}, 'Date', '%')
    subplot(3,3,7)
    plotTimeSeries(RiskHistoricaltable.Time, pnlPercByCategory{:,'Percpnl_Fund'}, ...
        'Returns In Funds', {'Return'}, 'Date', '%')
    subplot(3,3,[2,5,8])
    plotTimeSeries(RiskHistoricaltable.Time, pnlPercByCategory{:,'PercInNoCrypto'}, ...
        'Returns In Stocks, ETF and Funds', {'Return'}, 'Date', '%')   
    subplot(3,3,[3,6,9])
    plotTimeSeries(RiskHistoricaltable.Time, pnlPercByCategory{:,'Percpnl_Crypto'}, ...
        'Returns In Crypto', {'Return'}, 'Date', '%')    
    saveas(hg, char(strcat(pathfilewm, '\', 'ReturnsByCategory')), 'jpg')
        
    %% sold positions
    disp(['Creating sold positions table at ', char(datetime('now','format','HH:mm'))])
    % get sold positions
    hposSold = dailyinfoportf(dailyinfoportf.shares < 0, :);
    mean_Price_sold   = groupsummary(hposSold(:,{'Name', 'price'}),...
        {'Name'},'mean');
    sum_Quantity_sold = groupsummary(hposSold(:,{'Name', 'shares', ...
        'loc_quantity'}),{'Name'},'sum');
    sum_Quantity = groupsummary(dailyinfoportf(dailyinfoportf.shares > 0, ...
        {'Name', 'shares', 'loc_quantity'}),{'Name'},'sum');
    join_Quantity = outerjoin(sum_Quantity_sold(:,{'Name', 'sum_shares', ...
        'sum_loc_quantity'}), sum_Quantity(:,{'Name', 'sum_shares', ...
        'sum_loc_quantity'}),'type','left', 'keys',{'Name'}, ...
        'MergeKeys', true);
    join_Quantity.Properties.VariableNames = {'Name', 'Shares_Sold', ...
        'Loc_Value_Sold', 'Shares_Buyed', 'Loc_Value_Invested'};
    join_Quantity.Shares_Hold = join_Quantity.Shares_Buyed + ...
        join_Quantity.Shares_Sold ;
    join_Quantity.Loc_Value_Sold = -1 * join_Quantity.Loc_Value_Sold;
    join_Quantity.PnlValue = join_Quantity.Loc_Value_Sold - ...
        join_Quantity.Loc_Value_Invested;
    join_Quantity.RetPnlValue = (join_Quantity.Loc_Value_Sold ./...
        join_Quantity.Loc_Value_Invested) -1;
    join_Quantity.MiddlePrice =  join_Quantity.Loc_Value_Invested ./...
        join_Quantity.Shares_Buyed;
    join_Quantity.MeanPriceSold = mean_Price_sold.mean_price;
    join_Quantity.RetbyMiddlePrice = (join_Quantity.MeanPriceSold ./...
        join_Quantity.MiddlePrice) -1 ;
    soldPostable = table();
    % obtain profit/loss and returns
    for b = 1:height(join_Quantity)
        % get name
        name = join_Quantity.Name(b);
        soldPositiontable = table();
        soldPositiontable.Name = name;
        idxst = ismember(hposSold.Name, name);
        % obtain ticker
        tick = hposSold.symbol(idxst);
        % get prices and quantity
        idxp = ismember(dailyinfoportf.symbol, tick(1));
        prices = dailyinfoportf.price(idxp);
        quantity = dailyinfoportf.loc_quantity(idxp);
        % we search dividends only is entirely sold the positions
        [~, iddiv] = ismember(dividensTable.Name, tick(1));
        if sum(iddiv) > 0
            soldPositiontable.Valdiv = dividensTable.AccValDiv(iddiv > 0);
            if join_Quantity.Shares_Hold(b) == 0
                soldPositiontable.Valdiv = soldPositiontable.Valdiv;
            else
                soldPositiontable.Valdiv = 0;
            end
        else
            soldPositiontable.Valdiv = 0;
        end
        % TWR using price and MWR
        ret = price2ret(prices);
        periods = numel(ret);
        soldPositiontable.TWRuPrice = (prod(1+ ret))^(1/periods) - 1;
        soldPositiontable.MWR = irr([-quantity;soldPositiontable.Valdiv]);
        soldPostable = [soldPostable;soldPositiontable];
    end
    PosSoldTable = [join_Quantity,soldPostable(:,2:end)];
    
    %% Writetable
    disp(['Writing tables at ', char(datetime('now','format','HH:mm'))])
    % sheet 1
    % structure
    tbles1(1).tble = tableprices0t;
    tbles1(2).tble = tablepricesMiddlePrice;
    tbles1(3).tble = tablepricesDivndFee;
    tbles1(4).tble = timetable2table(tomorrowPriceticker);
    tbles1(5).tble = sumUpWorstDayTable;
    tbles1(6).tble = sumUpAverageValuesTable;
    tbles1(7).tble = PosSoldTable;
    
    % start writing the first table
    row = 1;
    columnLetters = char(xlscol(row));
    cellReference = sprintf('%s1', columnLetters);
    writetable(tableprices0t, ...
        char(strcat(pathfilewm, '\', 'Portfolio_PNL','_', ...
        datestr(today()) ,'.xlsx')), 'filetype',...
        'spreadsheet', 'sheet', 1, 'Range', cellReference)
    % write other tables
    row2 = 1;
    for k = 2:numel(tbles1)
        row2 = row2 + 2 + height(tbles1(k-1).tble);
        cellReference2 = strcat(columnLetters, num2str(row2));
        writetable(tbles1(k).tble, ...
            char(strcat(pathfilewm, '\', 'Portfolio_PNL','_', ...
            datestr(today()) ,'.xlsx')), 'filetype',...
            'spreadsheet', 'sheet', 1, 'Range', cellReference2)
    end
    
    % sheet 2
    % structure
    tbles2(1).tble  = riskTablet;
    tbles2(2).tble  = sumUpVaRTable;
    tbles2(3).tble  = sumUpESTable;
    tbles2(4).tble  = hVARaggsumUp;
    tbles2(5).tble  = hVARaggTLtable;
    tbles2(6).tble  = hVARaggPofTable;
    tbles2(7).tble  = hVARaggCCITable;
    tbles2(8).tble  = hVARaggTBFITable;
    tbles2(9).tble  = tableBacktestinghESaggWithId;
    tbles2(10).tble  = VARaggsumUp;
    tbles2(11).tble = VARaggTLtable;
    tbles2(12).tble = VARaggPofTable;
    tbles2(13).tble = VARaggCCITable;
    tbles2(14).tble = VARaggTBFITable;
    tbles2(15).tble = tableBacktestingESaggWithId;
    tbles2(16).tble = VARaggwVGsumUp;
    tbles2(17).tble = VARaggwVGTLtable;
    tbles2(18).tble = VARaggwVGPofTable;
    tbles2(19).tble = VARaggwVGCCITable;
    tbles2(20).tble = VARaggwVGTBFITable;
    tbles2(21).tble = tableBacktestingESwVGaggWithId;
    tbles2(22).tble = VARaggwTsumUp;
    tbles2(23).tble = VARaggwTTLtable;
    tbles2(24).tble = VARaggwTPofTable;
    tbles2(25).tble = VARaggwTCCITable;
    tbles2(26).tble = VARaggwTTBFITable;
    
    % write the first table
    row = 1;
    columnLetters2 = char(xlscol(row));
    cellReference2 = sprintf('%s1', columnLetters2);
    writetable(riskTablet, ...
        char(strcat(pathfilewm, '\','Portfolio_PNL','_',...
        datestr(today()),'.xlsx')), 'filetype',...
        'spreadsheet', 'sheet', 2, 'Range', cellReference2)
    % write other tables
    row3 = 1;
    for k = 2:numel(tbles2)
        row3 = row3 + 2 + height(tbles2(k-1).tble);
        cellReference3 = strcat(columnLetters, num2str(row3));
        writetable(tbles2(k).tble, ...
            char(strcat(pathfilewm, '\', 'Portfolio_PNL','_',...
            datestr(today()),'.xlsx')), 'filetype',...
            'spreadsheet', 'WriteRowNames', 1, 'sheet', 2, 'Range', ...
            cellReference3)
    end
    
    % sheet 5
    % structure
    tbles3(1).tble  = aggrPortfTable2w;
    tbles3(2).tble  = tblAssetType;
    tbles3(3).tble  = tblExch;
    tbles3(4).tble  = tblSector;
    tbles3(5).tble  = tblCountry;
    tbles3(6).tble  = monthlyQuantity;
    tbles3(7).tble  = monthlyQuantitywfee;
    
    % write the first table
    row = 1;
    columnLetters = char(xlscol(row));
    cellReference = sprintf('%s1', columnLetters);
    writetable(aggrPortfTable2w, ...
        char(strcat(pathfilewm, '\', 'Portfolio_PNL','_',...
        datestr(today()),'.xlsx')), 'filetype',...
        'spreadsheet', 'sheet', 5, 'Range', cellReference)
    % write other tables
    row4 = 1;
    for k = 2:numel(tbles3)
        row4 = 2 + row4 +  height(tbles3(k-1).tble);
        cellReference4 = strcat(columnLetters, num2str(row4));
        writetable(tbles3(k).tble, ...
            char(strcat(pathfilewm, '\', 'Portfolio_PNL','_',...
            datestr(today()),'.xlsx')), 'filetype',...
            'spreadsheet', 'sheet', 5, 'Range', cellReference4)
    end
    
    % write timetable
    % sheet 3
    writetimetable(RiskPortfolioAnalysis, ...
        char(strcat(pathfilewm, '\', 'Portfolio_PNL','_',...
        datestr(today()),'.xlsx')), 'filetype',...
        'spreadsheet', 'sheet', 3)
    % sheet 4
    writetimetable(VaRandESPortfolioAnalysis, ...
        char(strcat(pathfilewm, '\', 'Portfolio_PNL','_',...
        datestr(today()),'.xlsx')), 'filetype',...
        'spreadsheet', 'sheet', 4)
    % sheet 6
    writetimetable(simulatedMetricsTable, ...
        char(strcat(pathfilewm, '\', 'Portfolio_PNL','_',...
        datestr(today()),'.xlsx')), 'filetype',...
        'spreadsheet', 'sheet', 6)

    % sheet 7
    writetimetable(gnPnlByCategory, ...
        char(strcat(pathfilewm, '\', 'Portfolio_PNL','_',...
        datestr(today()),'.xlsx')), 'filetype',...
        'spreadsheet', 'sheet', 7)
    
    % open Activex server
    e = actxserver('Excel.Application');
    % open file (enter full path!)
    ewb = e.Workbooks.Open(char(strcat(pathfilewm, '\', 'Portfolio_PNL',...
        '_', datestr(today()),'.xlsx')));
    % rename sheets
    ewb.Worksheets.Item(1).Name  = 'Main_Values';
    ewb.Worksheets.Item(2).Name  = 'Risk_Analysis';
    ewb.Worksheets.Item(3).Name  = 'Portfolio_Benchmark';
    ewb.Worksheets.Item(4).Name  = 'VaR_and_ES_detail';
    ewb.Worksheets.Item(5).Name  = 'Portfolio_Values';
    ewb.Worksheets.Item(6).Name  = 'VaR_and_ES_Simulated_detail';
    ewb.Worksheets.Item(7).Name  = 'PNL_CategoryAsset'; 
    % save to the same file
    ewb.Save
    ewb.Close(true)
    e.Quit
    
else
    disp(['Portfolio_PNL file at ',...
        datestr(today()), ' already exist'])
end
% delete old portfolio_Pnl Analysis
if exist(char(strcat(pathfilewm, '\', 'Portfolio_PNL','_', ...
        datestr(today()-1) ,'.xlsx')), 'file')
    disp(['Deleting old portfolio file at ', char(datetime('now','format','HH:mm'))])
    delete(char(strcat(pathfilewm, '\', 'Portfolio_PNL','_', ...
        datestr(today()-1) ,'.xlsx')))
end

%% END
t = toc;
disp(['This is the end! at ', char(datetime('now','format','HH:mm'))])
disp(['This code is excecuted in ', num2str(t/60),' minutes'])
end