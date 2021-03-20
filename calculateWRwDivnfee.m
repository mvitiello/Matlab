function [TWRtwd, MWRtwd] = calculateWRwDivnfee(tickers, PnlTable,...
    dividendsTable, tableAcc, lastPrice2r)

% this function calculate the MWR and TWR of each asset in the portfolio
% input:
% tickers = list of the assets
% PnlTable = each asset pnl
% tableAcc = historical purchases
% dividendsTable = table with all dividends
% lastprice2r = the last price for each asset
% output = TWR and MWR for each assets (include dividends and fees)

% TWR
% get names
names2filteridxPnl = startsWith(PnlTable.Properties.VariableNames,...
    'PNL_');
names2filteridxPnl0 = startsWith(PnlTable.Properties.VariableNames,...
    'PNL0_');
names2filterPnl = PnlTable.Properties.VariableNames(names2filteridxPnl);
names2filterPnl(strcmp(names2filterPnl, 'PNL_USD_EUR_USDollarEuro')) = [];
names2filterPnl0 = PnlTable.Properties.VariableNames(names2filteridxPnl0);
names2filterPnl0(strcmp(names2filterPnl, 'PNL0_USD_EUR_USDollarEuro')) = [];

TWRtwd = table();
for k = 1: numel(names2filterPnl)
    % pnl for each asset
    Pnl = timetable(PnlTable.Time, PnlTable{:, names2filterPnl(k)});
    allpnl0 = PnlTable{:, names2filterPnl0(k)};
    % cut the data in each purchase
    cutoff = find(diff(allpnl0) > 0);
    PnlCutedt0 = Pnl(cutoff, :);
    PnlCutedt1 = [Pnl(cutoff+1, :);Pnl(end, :)];
    % join the pieces
    S = sortrows([PnlCutedt0;PnlCutedt1], 'Time');
    S(1,:) = [];
    name = split(names2filterPnl(k), '_');
    if numel(name) == 2
        iticker = name(2);
    elseif numel(name) == 3
        iticker = strcat(name(2), '.', name(3));
    else
        iticker = strcat(name(2), '/', name(3));
    end
    % add fees
    q = startsWith(tableAcc.symbol, iticker);
    tableAccf = tableAcc(q > 0, :);
    if isempty(tableAccf)
        iticker = char(iticker);
        q = startsWith(tableAcc.symbol, iticker(2:end));
        tableAccf = tableAcc(q > 0, :);
        iticker = cellstr(iticker);
    end
    tablefee = table2timetable(tableAccf(:,{'Date','fee'}));
    TT = synchronize(S,tablefee);
    idx = ismissing(TT(:,{'fee'}));
    TT{:,{'fee'}}(idx) = 0;
    S = timetable(TT.Time, sum(TT.Var1 + TT.fee,2, 'omitnan'));
    % add dividend if they have
    [~, p] = ismember(dividendsTable.symbol, iticker);
    if sum(p) > 0
        datedividend = dividendsTable.Date(p > 0);
        datevaluem = dividendsTable.Valuelm(p > 0);
        for g = 1: sum(p)
            finddate = S.Time - datedividend(g);
            [~, mindatep] = min(abs(finddate));
            S.Var1(mindatep) = S.Var1 (mindatep) + datevaluem(g);
        end
    end
    % find zeros and delete it
    [ii, jj] = find(~S.Var1);
    if jj
        S(ii,:) = [];
    end
    % calculate TWR
    if height(S) < 2
        rettwr = 0;
    else
        rettwr = price2ret(S.Var1);
    end
    ret2twr = rettwr(1:2:end) + 1;
    periods = numel(ret2twr);
    twr = prod(ret2twr)^(1/periods) - 1;
    twrtable = table(iticker,twr);
    TWRtwd = [TWRtwd;twrtable];
end

% MWR
purchases = tableAcc(:,{'Date','symbol', 'price', 'shares', 'fee'});
% calculate MWR
MWRtwd = table();
for k = 1:length(tickers)
    tick = tickers(k);
    tickIdx = ismember(purchases{:,2}, tick);
    % change the sign of the purchases
    pricesPurchasesTick = -purchases{tickIdx,3};
    sharesbyPurchases = purchases{tickIdx,4};
    tickIdxp = ismember(lastPrice2r{:,1}, tick);
    pricesPurchasesTickp = lastPrice2r{tickIdxp,2} * sum(sharesbyPurchases);
    % add the last pnl of the asset and the fees
    prices2irr = [(pricesPurchasesTick .* sharesbyPurchases) - ...
        purchases{tickIdx,5} ; pricesPurchasesTickp];
    iticker = tickers(k);
    [~, p] = ismember(dividendsTable.symbol, iticker);
    % add dividends if they have
    if sum(p) > 0
        [GroupsDividens, Gnames] = findgroups(dividendsTable.symbol);
        accdiv = splitapply(@sum, dividendsTable.Valuelm,...
            GroupsDividens);
        tablediv = table(Gnames, accdiv);
        [~, tick] = ismember(tablediv.Gnames, iticker);
        accdiv2irr =  tablediv.accdiv(tick > 0);
        prices2irr = [prices2irr; accdiv2irr];
    end
    % calculate MWR as IRR for each asset
    mwr = irr(prices2irr);
    tableMWRtick = table(iticker, mwr);
    MWRtwd = [MWRtwd; tableMWRtick];
end

end