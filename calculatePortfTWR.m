function portfTWR = calculatePortfTWR(PnlTable, totPnl0, ...
    tableAcc, dividendsTable, varargin)

% this function calculate the TWR of the portfolio
% input:
% PnlTable = each asset pnl
% totPnl0 = total value invested
% tableAcc = historical purchases
% dividendsTable = table with all dividends
% varargin = 'divnfee' include dividends and fee in the TWR
% output = portfolio's TWR

% Get total portfolio's PNL adding each asset pnl
names2filteridxPnl = startsWith(PnlTable.Properties.VariableNames,...
    'PNL_');
names2filterPnl = PnlTable.Properties.VariableNames(names2filteridxPnl);
PnlallData = PnlTable{:, names2filterPnl};
sumPnl =  timetable(PnlTable.Time,sum(PnlallData, 2, 'omitnan'));
% cut the data in each purchase
cutoff = find(diff(totPnl0) > 0);
PnlCutedt0 = sumPnl(cutoff, :);
PnlCutedt1 = [sumPnl(cutoff+1, :);sumPnl(end, :)];
% join the pieces
S = sortrows([PnlCutedt0;PnlCutedt1], 'Time');
S = unique(S);
if strcmp(varargin, 'divnfee')
    % add fees
    tablefee = table2timetable(tableAcc(:,{'Date','fee'}));
    TT = synchronize(S,tablefee);
    idx = ismissing(TT(:,{'fee'}));
    TT{:,{'fee'}}(idx) = 0;
    S = timetable(TT.Time, sum(TT.Var1 + TT.fee,2,'omitnan'));
    % add dividends
    datedividend = dividendsTable.Date;
    datevaluem = dividendsTable.value;
    for g = 1: numel(datedividend)
        finddate = S.Time - datedividend(g);
        [~, mindatep] = min(abs(finddate));
        S.Var1(mindatep) = S.Var1 (mindatep) + datevaluem(g);
    end
end
% find zeros
[ii, jj] = find(~S.Var1);
if jj
    S(ii,:) = [];
end
% calculate TWR
rettwr = price2ret(S.Var1);
ret2twr = rettwr(1:2:end) + 1;
periods = numel(ret2twr);
portfTWR = prod(ret2twr)^(1/periods) - 1;
end