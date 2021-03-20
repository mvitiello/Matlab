% code to update  funds positions
% required: .xls from R4
pathfilewm = cd();
xlsfileacc = 'Account.xlsx';
xlsfileaccfund = 'fondos.xls';
xlsfileaccpath = fullfile(pathfilewm, xlsfileacc);
xlsfileaccfundpath = fullfile(pathfilewm, xlsfileaccfund);
% read tables
fundData = readtable(xlsfileaccfundpath, 'filetype',...
    'spreadsheet');
% delete first 6 lines (empty innecessary data)
fundData(1:6,:) = [];
% variables names adding 'Tikcker'
titles = fundData{1,:};
% find n of account finding entire missing lines
skipacc = find(ismissing(fundData(:,1)));
nAccounts = size(skipacc,1) + 1 ;
% first account
firstAcc = fundData{2,1};
% get info first account
firstAccValues = fundData(3:skipacc(1)-1,:);
firstAccValuestrtt = obtainTable(firstAccValues,titles, firstAcc);
% second account
secondAcc = fundData{skipacc(1)+1,1};
% get info
secondAccValues = fundData(skipacc(1)+2:skipacc(2)-1,:);
secondAccValuestrtt = obtainTable(secondAccValues,titles, secondAcc);
% third account
thirdAcc = fundData{skipacc(2)+1,1};
thirdAccValues = fundData(skipacc(2)+2:end,:);
thirdAccValuestrtt = obtainTable(thirdAccValues,titles, thirdAcc);
% join tables
AccountsJoined = sortrows([firstAccValuestrtt; secondAccValuestrtt;...
    thirdAccValuestrtt], 'AccValuesD', 'ascend');
% get info for account
account = readtable(xlsfileaccpath, 'filetype', ...
    'spreadsheet', 'sheet', 1);
idxF = strcmp(account.Category1, 'Fund');
fundAccInvsAll = account(idxF,:);
idxFR4 = strcmp(fundAccInvsAll.Orig, 'RENTA4');
fundAccInvs = fundAccInvsAll(idxFR4,:);

% check new positions
maxDateInFund = max(AccountsJoined.AccValuesD);
maxDateFundAcc = max(fundAccInvs.Date);
if (maxDateInFund > maxDateFundAcc)
    postDates = AccountsJoined.AccValuesD > maxDateFundAcc;
    newPosInFunds = AccountsJoined(postDates, :);
    % check that new positions belong to my portfolio (exclude MSS GLOBAL BRANDS)
    fundsInAccount = unique(fundAccInvs{:,'Name'});
    % identify positions
    identifynewPosInAcc = ismember(newPosInFunds.Ticker, fundsInAccount);
    newPosInacc = newPosInFunds(identifynewPosInAcc,:);
    newPosInaccUpdate = table();
    % add positions
    for k = 1:height(newPosInacc)
        % get ticker 
        tick = newPosInacc.Ticker(k);
        % get tickers in fund
        posFundInacc = fundAccInvs(strcmp(tick, fundAccInvs.Name),:);
        % get a positions in Account as example
        posFundInaccEx = posFundInacc(end,:);
        % get position in Fund
        newPosInaccToUpdate = newPosInacc(k,:);
        % change inputs
        posFundInaccEx.Date = datetime(newPosInaccToUpdate.AccValuesD,...
            'format','MM/dd/yyyy');
        posFundInaccEx.loc_quantity = newPosInaccToUpdate.ImporteNETO;
        posFundInaccEx.quantity = newPosInaccToUpdate.ImporteBrutoDiv_;
        posFundInaccEx.shares = newPosInaccToUpdate.Participaciones;
        posFundInaccEx.price = posFundInaccEx.quantity/posFundInaccEx.shares;
        posFundInaccEx.fx = (posFundInaccEx.loc_quantity - ...
            newPosInaccToUpdate.Comisi_n_)/posFundInaccEx.quantity;
        posFundInaccEx.fee = newPosInaccToUpdate.Comisi_n_;
        posFundInaccEx.fee___ = 100*(posFundInaccEx.fee / posFundInaccEx.loc_quantity);
        newPosInaccUpdate = [newPosInaccUpdate;posFundInaccEx];
    end
    % writetable if is not a empty table 
    if ~isempty(newPosInaccUpdate)
        accountUpdated = [account;newPosInaccUpdate];
        writetable(accountUpdated, xlsfileaccpath, 'filetype',...
            'spreadsheet', 'sheet', 1);
        disp('Updated new positions in Funds')
    end 
else
    disp('There is no new positions in Funds')
end

% INTERNAL FUNCTION

function ttv = obtainTable(AccValues,titles, ticker)
% change to numerical
AccValuesN = str2double(AccValues{:,3:8});
% change datetime
AccValuesD = datetime(AccValues{:,1},'InputFormat', 'dd/MM/yyyy');
% convert to table
AccValuestr = array2table(AccValuesN);
AccValuestr.Properties.VariableNames = matlab.lang.makeValidName(cellstr(titles(3:8)));
% timetable
AccValuestrtt = table2timetable([table(AccValuesD),AccValuestr]);
% get number of operations to add a ticker identifier
nopsfirstAcc = height(AccValuestr);
% array with a account identifier
firstaccountticker = repmat(ticker,nopsfirstAcc,1);
% table with first acc
AccValuestrtt.Ticker = firstaccountticker;
ttv = AccValuestrtt;
end
