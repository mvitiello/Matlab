function [tablePPD, tableTechAnalysis] = Tecnicallyanalysis(pricestable,...
    period, path2save, ticker)
% Input
% pricetable =  table with High, Low, Close (weekly) and Open prices
% dailyClose = daily Closing prices;
% period = period of the analysis
% Output
% tablePPD = table with pivot points
% tableTechAnalysis = table with different technical indicators


% Weekly Closing, High and Low prices
H = pricestable.High;
L = pricestable.Low;
O = pricestable.Open;
try 
    C = pricestable{:, matlab.lang.makeValidName(ticker)};
catch ME
    C = pricestable{:,2};
end 

prices = C; 

prices2Candle = timetable(pricestable.Date, O, H, L, C);
prices2Candle.Properties.VariableNames = {'Open', 'High', 'Low','Close'};
% plot candlestick 
g = figure('visible', 'off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
candle(prices2Candle(end-31:end,:) , 'r');
title(['Candlestick chart for ', char(ticker)])
saveas(g, char(strcat(path2save, '\', ...
        'Candlestick chart')), 'jpg')

%% A) Pivot Points (PP)

% A.1) standard Pivot Points
PPs = 1/3*(H + L + C);
S1s = 2*PPs - H ;
R1s = 2*PPs - L ;
S2s = PPs - (H - L);
R2s = PPs +(H - L);
S3s = S2s;
S4s = S3s;
R3s = R2s;
R4s = R3s;
classicPP = [S4s  S3s S2s S1s PPs R1s R2s R3s R4s];
% A.2) Woodie's Pivot Points
PPw = 1/4*((H + L) + 2*C);
R2w = PPw + (H - L);
R1w = 2*PPw - L;
S1w = 2*PPw - H;
S2w = PPw - (H - L);
S3w = S2w;
S4w = S3w;
R3w = R2w;
R4w = R3w;
woodiePP = [S4w  S3w S2w S1w PPw R1w R2w R3w R4w];
% A.3) Camarilla Pivot points
PPc = PPs;
R4c = 0.55*(H - L) + C;
R3c = 0.275*(H - L) + C;
R2c = 0.183*(H - L) + C;
R1c = 0.0916*(H - L) + C;
S1c = C - 0.0916*(H - L);
S2c = C - 0.183*(H - L);
S3c = C - 0.275*(H - L);
S4c = C - 0.55*(H - L);
camarillaPP = [S4c  S3c S2c S1c PPc R1c R2c R3c R4c];
% A.4) Fibonacci Pivot Points
PPf = PPs;
R1f = PPf + (.382 * (H - L));
R2f = PPf + (.618 * (H - L));
R3f = PPf + (H - L);
S1f = PPf - (.382 * (H - L));
S2f = PPf - (.618 * (H - L));
S3f = PPf - (H - L);
S4f = S3f;
R4f = R3f;
fibonacciPP = [S4f  S3f S2f S1f PPf R1f R2f R3f R4f];
% A.5) Demark Pivot Points
if C(end) > O (end)
    X = 2*H + L + C;
elseif C(end) == O(end)
    X = H + L + 2*C;
else
    X = H + 2*L + C;
end
PPd = X./4;
R1d = X./2 - L;
S1d = X./2 - H;
S2d = S1d;
S3d = S2d;
S4d = S3d;
R2d = R1d;
R3d = R2d;
R4d = R3d;
demarkPP = [S4d  S3d S2d S1d PPd R1d R2d R3d R4d];
% table with Pivot points
tablePP = array2table([classicPP(end,:);woodiePP(end, :); camarillaPP(end,:); ...
    fibonacciPP(end,:); demarkPP(end,:)]);
tablePP.Properties.VariableNames = {'S4', 'S3', 'S2', 'S1', 'PP', 'R1',...
    'R2', 'R3', 'R4'};
tablePP.Properties.RowNames = {'Classic', 'Woodie', 'Camarilla',...
    'Fibonacci', 'Demark'};

% "Quick" technical investment descision
for k = 1:height(tablePP)
    if prices(end) > tablePP{k,1} || prices(end) > tablePP{k,2} ||...
            prices(end) > tablePP{k,3} || prices(end) > tablePP{k,4}
        PPinvDes(k,1)  = {'Hold/Buy'};
        PPinvDes(k,2) = PPinvDes(k,1);
        PPinvDes(k,3) = PPinvDes(k,1);
        PPinvDes(k,4) = PPinvDes(k,1);
    elseif prices(end) < tablePP{k,1} || prices(end) < tablePP{k,2} ||...
            prices(end) < tablePP{k,3} || prices(end) < tablePP{k,4}
        PPinvDes(k,1)  = {'Sell'};
        PPinvDes(k,2) = PPinvDes(k,1);
        PPinvDes(k,3) = PPinvDes(k,1);
        PPinvDes(k,4) = PPinvDes(k,1);
    elseif prices(end) > tablePP{k,6} || prices(end) > tablePP{k,7} ||...
            prices(end) > tablePP{k,8} || prices(end) > tablePP{k,9}
        PPinvDes(k,6)  = {'Hold/Sell'};
        PPinvDes(k,7) = PPinvDes(k,6);
        PPinvDes(k,8) = PPinvDes(k,6);
        PPinvDes(k,9) = PPinvDes(k,6);
    else prices(end) < tablePP{k,6} || prices(end) < tablePP{k,7} ||...
            prices(end) < tablePP{k,8} || prices(end) < tablePP{k,9};
        PPinvDes(k,6)  = {'Buy'};
        PPinvDes(k,7) = PPinvDes(k,6);
        PPinvDes(k,8) = PPinvDes(k,6);
        PPinvDes(k,9) = PPinvDes(k,6);
    end
end
PPinvDes = reshape(PPinvDes(~cellfun('isempty',PPinvDes)),5,4);
for k = 1:height(tablePP)
[sPP,~,jPP] = unique(PPinvDes(k,:));
mostcommonPP(k) = sPP(mode(jPP));
end 
tableD = cell2table(mostcommonPP');
tableD.Properties.VariableNames = {'Action'};
closeT = table(C(end)*ones(5,1));
closeT.Properties.VariableNames = {'ClosePrice'};
tablePPD = [tablePP, tableD];

tablePP = [tablePP, closeT];

g1 = figure('visible', 'off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
X = categorical({'Classic', 'Woodie', 'Camarilla',...
    'Fibonacci', 'Demark'});
bar(X, tablePP{:,:})
legend('S4', 'S3', 'S2', 'S1', 'PP', 'R1',...
    'R2', 'R3', 'R4', 'ClosePrice', 'Location','northeastOutside')
title(['Pivot Prices - ',char(ticker)])
saveas(g1, char(strcat(path2save, '\', ...
        'PivotPrices')), 'jpg')
%% B) Indicators

% B.1) RSI - Relative Strength Index
% absolute gain and losses
diffdata    = diff(prices);
indwe = (diffdata == 0);
diffdata(indwe) = [];
pricechange = abs(diffdata);
advances = pricechange;
declines = pricechange;
advances (diffdata < 0) = 0;
declines (diffdata >= 0) = 0;

for k = period:size(diffdata, 1)
    % total Gains/losses
    totalGain = sum(advances((k - (period-1)):k));
    totalLoss = sum(declines((k - (period-1)):k));
    % Calculate RSI
    rs = totalGain ./ totalLoss;
    trsi(k) = 100 - (100 / (1+rs));
end

% technical investment descision - RSI
for k = 1:length(trsi)
    if trsi(k) < 30 && trsi(k) >  0
        RSIinvDes(k) = {'Sell'};
    elseif trsi(k) > 70
        RSIinvDes(k) = {'Sell'};
    elseif trsi(k) == 0
        RSIinvDes(k) = {'Neutral'};
    else
        RSIinvDes(k) = {'Buy'};
    end
end

% get the dates
DateReturns = pricestable.Date; 
% first date
startDate   = datenum(DateReturns(1));
% last date
endDate     = datenum(DateReturns(size(DateReturns, 1)));
% axes with dates
xData       = linspace(startDate, endDate, size(DateReturns, 1));             
g2 = figure('visible', 'off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
% plot, legend, title... 
plot(xData(~(indwe)),[diffdata, trsi']);
legend('absReturns','RSI')
xlim([startDate endDate])
datetick('x', 'ddmmyyyy', 'keeplimits');
title(['Relative Strength Index - ',char(ticker)])
saveas(g2, char(strcat(path2save, '\', ...
        'RSI')), 'jpg')

tableRSI = createTableTechIndicator({'RSI'}, trsi, RSIinvDes);

% B.2) Stochastic Oscillator(k, d) and Stochastic RSI

ko = 9;
d = 6;
osch = (prices - min(prices(end-ko:end)))/(max(prices(end-ko:end))- ...
    min(prices(end-ko:end)));
kdOsch = movmean(osch, d);
signOsch = osch - kdOsch;

% technical investment descision - Stochastic Oscillator
for k = 1:length(signOsch)
    if signOsch(k) > 0
        oschinvDes(k) = {'Buy'};
    else
        oschinvDes(k) = {'Sell'};
    end
end

% get the dates
DateReturns = pricestable.Date; 
% first date
startDate   = datenum(DateReturns(1));
% last date
endDate     = datenum(DateReturns(size(DateReturns, 1)));
% axes with dates
xData       = linspace(startDate, endDate, size(DateReturns, 1));             
% plot, legend, title... 
g3 = figure('visible', 'off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);

plot(xData, signOsch); 
xlim([startDate endDate])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title(['Stochastic Oscillator(9, 6) - ',char(ticker)])
saveas(g3, char(strcat(path2save, '\', ...
        'StochasticOscillator')), 'jpg')

tableOsch = createTableTechIndicator({'Osch'}, signOsch, oschinvDes);

stochRSI = (trsi - min(trsi(end-ko:end)))/(max(trsi(end-ko:end))- ...
    min(trsi(end-ko:end)));

% technical investment descision - Stochastic RSI
for k = 1:length(stochRSI)
    if stochRSI(k) < 2
        RSIStochinvDes(k) = {'Sell'};
    elseif stochRSI(k) > 3
        RSIStochinvDes(k) = {'Sell'};
    else
        RSIStochinvDes(k) = {'Buy'};
    end
end

tableRSIOsch = createTableTechIndicator({'RSIOsch'}, stochRSI, RSIStochinvDes);

% B.3) MACD and signal

MACD = tmovavg(prices', 'e', 12) -  tmovavg(prices', 'e', 26);
MACD(1:25) = 0;
signal = tmovavg(MACD, 'e', 9);
MACDh = MACD - signal;

% technical investment descision - MACDh
for k = 1:length(MACDh)
    if MACDh(k) > 0
        macdhinvDes(k) = {'Buy'};
    else
        macdhinvDes(k) = {'Sell'};
    end
end

tableMacdh = createTableTechIndicator({'Macdh'}, MACDh, macdhinvDes);

% MA standard
MAs20 = cummeanv(prices, 20);
MAs100 = cummeanv(prices, 100);
signMAs = MAs20 - MAs100;

% technical investment descision - MA
for k = 1:length(signMAs)
    if signMAs(k) > 0
        MAsinvDes(k) = {'Buy'};
    else
        MAsinvDes(k) = {'Sell'};
    end
end

tableMAs = createTableTechIndicator({'MAs'}, signMAs, MAsinvDes);

% MA exponential
MAe20 = tmovavg(prices', 'e', 20);
MAe100 = tmovavg(prices', 'e', 100);
signMAe = MAe20 - MAe100;

% technical investment descision - MAE
for k = 1:length(signMAe)
    if signMAe(k) > 0
        MAeinvDes(k) = {'Buy'};
    else
        MAeinvDes(k) = {'Sell'};
    end
end

tableMAe = createTableTechIndicator({'MAe'}, signMAe, MAeinvDes);

% HMA
HMA20 = cummeanv(2*movmean(prices, 10) - cummeanv(prices, 20),...
    round(sqrt(20),0));
HMA100 = cummeanv(2*cummeanv(prices, 50) - cummeanv(prices, 100),...
    round(sqrt(100),0));
signHMA = HMA20 - HMA100;

% technical investment descision - HMA
for k = 1:length(signHMA)
    if signHMA(k) > 0
        HMAinvDes(k) = {'Buy'};
    else
        HMAinvDes(k) = {'Sell'};
    end
end

tableHMA = createTableTechIndicator({'HMA'}, signHMA, HMAinvDes);

% plot, legend, title... 
g4 = figure('visible', 'off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
subplot(2,1,1)
plot(xData, [prices,signMAs,signMAe', signHMA]);
legend('Price', 'MAs', 'MAe', 'HullMA', 'location','best')
xlim([startDate endDate])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title(['Prices and Signal - ',char(ticker)])
subplot(2,1,2)
bar(xData, MACDh);
xlim([startDate endDate])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title('MACDh')
saveas(g4, char(strcat(path2save, '\', ...
        'MACD_MA')), 'jpg')

%% OHL analysis 

% B.4) ADX() - Average Directional Index,  Wilder's DMI (ADX)
observ = height(pricestable);
% true range
hml = H-L;                                  % high - low
hmc = [0;abs(H(2:observ)- C(1:observ-1))];  % abs(high - close)
lmc = [0;abs(L(2:observ)- C(1:observ-1))];  % abs(low - close)
tr = max([hml,hmc,lmc],[],2);               % true range

% Directional Movement
hmh = H(2:observ)- H(1:observ-1);           % high - high
lml = L(1:observ-1)- L(2:observ);           % low - low
pdm1  = zeros(observ-1,1);                  % preallocate pdm1
maxh = max(hmh,0);
pdm1(hmh > lml) = maxh(hmh > lml);          % plus
mdm1  = zeros(observ-1,1);                  % preallocate mdm1
maxl = max(lml,0);
mdm1(lml > hmh) = maxl(lml > hmh);          % minus
pdm1 = [nan;pdm1];
mdm1 = [nan;mdm1];

% Preallocate 14 period tr, pdm, mdm, adx
trperiod  = zeros(observ,1);  % period true range
pdmperiod = trperiod;       % period plus directional movement
mdm1period = trperiod;      % period minus directional movement
adx   = trperiod;           % average directional index

% Calculate ADX
trperiod(period+1)  = sum(tr(period+1-period+1:period+1), 'omitnan');
pdmperiod(period+1) = sum(pdm1(period+1-period+1:period+1), 'omitnan');
mdm1period(period+1) = sum(mdm1(period+1-period+1:period+1), 'omitnan');
for k = period+2:observ
    trperiod(k)  = trperiod(k-1)-trperiod(k-1)/period+tr(k);
    pdmperiod(k) = pdmperiod(k-1)-pdmperiod(k-1)/period+pdm1(k);
    mdm1period(k) = mdm1period(k-1)-mdm1period(k-1)/period+mdm1(k);
end
pdiperiod = 100*pdmperiod./trperiod;   % period plus directional indicator
mdiperiod = 100*mdm1period./trperiod;  % period minus directional indicator
% directional movement index
dmx   = 100*abs(pdiperiod-mdiperiod)./(pdiperiod+mdiperiod);
% Average Directional Index
adx(2*period) = sum(dmx(period+1:2*period))/(2*period-period-1);
for k = 2*period+1:observ
    adx(k) = (adx(k-1)*(period-1)+dmx(k))/period;
end

% technical investment descision - ADX
for k = 1:length(adx)
    if adx(k) < 30
        ADXinvDes(k) = {'Sell'};
    elseif adx(k) > 50
        ADXinvDes(k) = {'Sell'};
    else
        ADXinvDes(k) = {'Buy'};
    end
end
% get the dates
DateReturnsOHL = pricestable.Date; 
% first date
startDateOHL   = datenum(DateReturnsOHL(1));
% last date
endDateOHL     = datenum(DateReturnsOHL(size(DateReturnsOHL, 1)));
% axes with dates
xDataOHL       = linspace(startDateOHL, endDateOHL, size(DateReturnsOHL, 1));             

% plot, legend, title...
g5 = figure('visible', 'off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
subplot(2,3,1)
plot(xDataOHL, adx); 
xlim([startDateOHL endDateOHL])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title(['Average Directional Index - ',char(ticker)])

tableADX = createTableTechIndicator({'ADX'}, adx, ADXinvDes);

% B.5)  % William's %R
% Highest High and Lowest Low
llv = zeros(observ,1);                % preallocate lowest low
llv(1:period) = min(L(1:period));     % lowest low of first kperiods
for k = period:observ                 % cycle through rest of data
    llv(k) = min(L(k-period+1:k));    % lowest low of previous kperiods
end
hhv = zeros(observ,1);                % preallocate highest high
hhv(1:period) = max(H(1:period));     % highest high of first kperiods
for k = period:observ                 % cycle through rest of data
    hhv(k) = max(H(k-period+1:k));    % highest high of previous kperiods
end

% Williams %R
wpctr = nan(observ,1);
nzero = find((hhv-llv) ~= 0);
wpctr(nzero) = ((hhv(nzero)-C(nzero))./(hhv(nzero)-llv(nzero))) * -100;

% technical investment descision - Williams'R%
for k = 1:length(wpctr)
    if abs(wpctr(k)) < 20
        WRinvDes(k) = {'Sell'};
    elseif abs(wpctr(k)) > 80
        WRinvDes(k) = {'Sell'};
    else
        WRinvDes(k) = {'Buy'};
    end
end

subplot(2,3,2)
plot(xDataOHL, wpctr); 
xlim([startDateOHL endDateOHL])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title(['Williams %R - ',char(ticker)])

tableWR = createTableTechIndicator({'WR'}, wpctr, WRinvDes);

% B.6) CCI, Commodity Channel Index
% Typical Price and a constant
tp = (H + L + C)/3;
const  = 0.015;
% Simple moving average of typical price
smatp = movmean(tp, period);
% Sum of the mean absolute deviation
smad = nan(observ,1);
cci  = smad;    % preallocate cci
for k = period:observ
    smad(k) = sum(abs(smatp(k)-tp(k-period+1:k)));
end

% Commodity Channel Index
k = period:observ;
cci(k) = (tp(k)-smatp(k))./(const*smad(k)/period);

% technical investment descision - CCI
for k = 1:length(cci)
    if cci(k) > 50
        CCIinvDes(k) = {'Buy'};
    else
        CCIinvDes(k) = {'Sell'};
    end
end

subplot(2,3,3)
plot(xDataOHL, cci); 
xlim([startDateOHL endDateOHL])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title(['Commodity Channel Index - ',char(ticker)])

tableCCI = createTableTechIndicator({'CCI'}, cci, CCIinvDes);

% B.7) Average True Range ATR and volatility ratio
ATR = tmovavg(tr', 'e', period);
vr = tr./ATR';
averageVol = mean(vr, 'omitnan');

% technical investment descision - Vol Ratio
for k = 1:length(vr)
    if vr(k) > averageVol
        VRinvDes(k) = {'Sell'};
    else
        VRinvDes(k) = {'Buy'};
    end
end

subplot(2,3,4)
plot(xDataOHL, vr); 
xlim([startDateOHL endDateOHL])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title(['Volatility Ratio - ',char(ticker)])

tableVolRatio = createTableTechIndicator({'VolRatio'}, vr, VRinvDes);

% B.8 ) ROC
roc = nan(observ,1);
roc(period+1:observ) = ((C(period+1:observ)- ...
    C(1:observ-period))./C(1:observ-period))*100;

minRoc = min(roc);
% technical investment descision - ROC
for k = 1:length(roc)
    if roc(k) >= minRoc
        rocinvDes(k) = {'Buy'};
    else
        rocinvDes(k) = {'Sell'};
    end
end

subplot(2,3,5)
plot(xDataOHL, roc); 
xlim([startDateOHL endDateOHL])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title(['Rate Of Change - ',char(ticker)])

tableROC = createTableTechIndicator({'ROC'}, roc, rocinvDes);

% B.9) HHLL, Highest High, Lowest Low

% Lowest Low
llv = nan(observ,1);                      % preallocate lowest low
llv(1:period) = min(L(1:period));         % lowest low of first kperiods
for k = period:observ                     % cycle through rest of data
    llv(k) = min(L(k-period+1:k));        % lowest low of previous kperiods
end

% Highest High
hhv = nan(observ,1);                      % preallocate highest high
hhv(1:period) = max(H(1:period));         % highest high of first kperiods
for k = period:observ                     % cycle through rest of data
    hhv(k) = max(H(k-period+1:k));        % highest high of previous kperiods
end

% Midpoint
mp = (hhv+llv)/2;
% Format Output
HHL = [hhv C llv];

% technical investment descision - HHLL
for k = 1:length(HHL)
    if C(k) > hhv(k)
        HHLinvDes(k) = {'Sell'};
    elseif C(k) < llv(k)
        HHLinvDes(k) = {'Buy'};
    else
        HHLinvDes(k) = {'Neutral'};
    end
end

tableHHLL = createTableTechIndicator({'HHLL'}, mp, HHLinvDes);

% B.10) BB power, Bulls and Bears Power
bearsPower = L - tmovavg(C', 'e', period)';
bullPower = H - tmovavg(C', 'e', period)';
BBP = bullPower./bearsPower ;

% technical investment descision - BBP
for k = 1:length(BBP)
    if BBP(k) > 0
        BBPinvDes(k) = {'Buy'};
    else
        BBPinvDes(k) = {'Sell'};
    end
end

subplot(2,3,6)
plot(xDataOHL, BBP); 
xlim([startDateOHL endDateOHL])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title(['Bulls and Bears Power - ',char(ticker)])
saveas(g5, char(strcat(path2save, '\', ...
        'SomeIndicators')), 'jpg')

tableBBP = createTableTechIndicator({'BBP'}, BBP, BBPinvDes);

% B.11) BB, Bollinger Bands

weight = 0;
nstd   = 2;
% Create output vectors.
mid  = nan(size(C, 1), 1);
uppr = mid;
lowr = mid;

% Create weight vector.
wtsvec = ((1:period).^weight) ./ (sum((1:period).^weight));

% Save the original data and remove NaN's from the data to be processed.
nnandata = C(~isnan(C));

% Calculate middle band moving average using convolution.
cmid    = conv(nnandata, wtsvec);
nnanmid = cmid(period:length(nnandata));

% Calculate shift for the upper and lower bands. The shift is a
% moving standard deviation of the data.
mstd = nnandata(period:end); % Pre-allocate
for k = period:length(nnandata)
    mstd(k-period+1, :) = std(nnandata(k-period+1:k));
end

% Calculate the upper and lower bands.
nnanuppr = nnanmid + nstd.*mstd;
nnanlowr = nnanmid - nstd.*mstd;

% Return the values.
nanVec = nan(period-1,1);
mid(~isnan(C))  = [nanVec; nnanmid];
uppr(~isnan(C)) = [nanVec; nnanuppr];
lowr(~isnan(C)) = [nanVec; nnanlowr];

% out BB
BB = [uppr, mid, C ,lowr];

% technical investment descision - BB
for k = 1:length(HHL)
    if C(k) > uppr(k)
        BBinvDes(k) = {'Sell'};
    elseif C(k) < lowr(k)
        BBinvDes(k) = {'Buy'};
    else
        BBinvDes(k) = {'Neutral'};
    end  
end

g6 = figure('visible', 'off', 'units', 'normalized',...
        'outerposition', [0 0 1 1]);
plot(xDataOHL, BB);
legend('Upper', 'Middle', 'Close', 'Lower')
xlim([startDateOHL endDateOHL])
datetick('x', 'ddmmyyyy', 'keeplimits'); 
title(['Bollinger Bands - ',char(ticker)])
saveas(g6, char(strcat(path2save, '\', ...
        'Bollinger Bands')), 'jpg')

tableBB = createTableTechIndicator({'BB'}, mid, BBinvDes);

% Final table
tableTechAnalysis = [tableRSI; tableOsch; tableRSIOsch; tableMacdh; ...
    tableMAs; tableMAe; tableHMA; tableADX; tableWR; tableCCI;... 
    tableVolRatio; tableROC; tableHHLL; tableBBP; tableBB];

end

function tableout = createTableTechIndicator(Indicator,techval, techDes)
Val = techval(end);
Des = techDes(end);
[s,~,j] = unique(Des);
mostCommon = s(mode(j));
tableout = table(Indicator, Val, Des, mostCommon);
end