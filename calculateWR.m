function [TWRt, MWRt] = calculateWR(tickers, PnlTable, tableAcc,...
    lastPrice2r)

% this function calculate the TWR and MWR for each asset in the portfolio
% input:
% tickers = list of the assets
% PnlTable = each asset pnl
% tableAcc = historical purchases
% lastprice2r = the last price for each asset
% output = TWR and MWR for each assets

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

TWRt = table();
for k = 1: numel(names2filterPnl)
    % pnl for each asset
    Pnl = timetable(PnlTable.Time, PnlTable{:, names2filterPnl(k)});
    allpnl0 = PnlTable{:, names2filterPnl0(k)};
    % cut the data in each purchase
    cutoff = find(diff(allpnl0) > 0);
    PnlCutedt0 = Pnl(cutoff, :);
    PnlCutedt1 = [Pnl(cutoff+1, :); Pnl(end, :)];
    % join the pieces
    S = sortrows([PnlCutedt0;PnlCutedt1], 'Time');
    S(1,:) = [];
    name = split(names2filterPnl(k), '_');
    if numel(name) == 2
        iticker = name(2);
    elseif numel(name) == 3
        iticker = strcat(name(2), '.', name(3));
    else
        iticker = strcat(name(2), '/', name(3), ' - ',name(4));
    end
    % calculate the TWR
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
    TWRt = [TWRt;twrtable ];
end

purchases = tableAcc(:,{'Date','symbol', 'price', 'shares'});%

% calculate MWR
MWRt = table();
for k = 1:length(tickers)
    tick = tickers(k);
    tickIdx = ismember(purchases{:,2}, tick);
    % change the sign of the purchases
    pricesPurchasesTick = -purchases{tickIdx,3};
    sharesbyPurchases = purchases{tickIdx,4};
    tickIdxp = ismember(lastPrice2r{:,1}, tick);
    pricesPurchasesTickp = lastPrice2r{tickIdxp,2} * sum(sharesbyPurchases);
    % add the last pnl of the asset
    prices2irr = [pricesPurchasesTick .* sharesbyPurchases ;...
        pricesPurchasesTickp];
    % Calculate MWR as IRR
    mwr = irr(prices2irr);
    tableMWRtick = table(tick, mwr);
    MWRt = [MWRt; tableMWRtick];
end

end