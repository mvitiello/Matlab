% code to update  stok and ETF positions
% required: .xls from DEGIRO
% code to update  funds positions
% required: .xls from R4
pathfilewm = cd();
xlsfileacc = 'Account.xlsx';
xlsfileaccStockndETF = 'Transactions.xls';
xlsfileaccpath = fullfile(pathfilewm, xlsfileacc);
xlsfileaccStockndETFpath = fullfile(pathfilewm, xlsfileaccStockndETF);
% read tables
StockndETFData = readtable(xlsfileaccStockndETFpath, 'filetype',...
    'spreadsheet');
% replace '.' by ',' and change to numerical
AccValuesN = strrep(StockndETFData{:,{'Precio', 'ValorLocal', ...
    'Valor', 'CostesDeTransacci_n', 'Total'}}, '.' , '');
AccValuesN = strrep(AccValuesN, ',' , '.');
% change to double and change the sign: (for us) buys are positive
% and sells are negative
AccValuesN = str2double(AccValuesN);
AccValuesN(:,2:end) = -AccValuesN(:,2:end);
% add number of shares
AccValuesN = [str2double(StockndETFData{:,'N_mero'}),AccValuesN];
% change datetime
AccValuesD = datetime(StockndETFData{:,1},'InputFormat', 'dd-MM-yyyy');
% convert a table
AccValuestr = array2table(AccValuesN);
% add variable names
AccValuestr.Properties.VariableNames = {'Shares', 'Price', 'Loc_quantity',...
    'quantity', 'fee', 'Total'};
% timetable
AccValuestrtt = table2timetable([table(AccValuesD),AccValuestr]);
AccValuestrtt.ISIN = StockndETFData.ISIN;
% get info for account
account = readtable(xlsfileaccpath, 'filetype', ...
    'spreadsheet', 'sheet', 1);
idxS = strcmp(account.Category1, 'Stock');
idxF = strcmp(account.Category1, 'ETF');
idx = logical(idxS + idxF);
AccInvs = account(idx,:);
% check new positions
maxDateInBroker = max(AccValuestrtt.AccValuesD);
maxDateAcc = max(AccInvs.Date);
if maxDateInBroker > maxDateAcc
    postDates = AccValuestrtt.AccValuesD > maxDateAcc;
    newPos = AccValuestrtt(postDates, :);
    % check that new positions belong to my portfolio
    InstrInAccount = unique(AccInvs{:,'ISIN'});
    % identify positions
    identifynewPosInAcc = ismember(newPos.ISIN, InstrInAccount);
    if any(~identifynewPosInAcc)
        newPosNID = newPos(~identifynewPosInAcc,:);
        productNID = StockndETFData.Producto(ismember(StockndETFData.ISIN,...
            newPosNID{:,'ISIN'}));
        disp([repmat('The following position is new: ', ...
            size(productNID, 1), 1),char(newPosNID{:,...
            'ISIN'}),repmat('(',size(productNID,1),1), char(productNID),...
            repmat(')',size(productNID,1),1)])
    end
    newPosInacc = sortrows(newPos(identifynewPosInAcc,:), 'AccValuesD',...
        'ascend');
    % if this table is empty this means that the Isin is new, so
    % we can not update because we need more info
    if isempty(newPosInacc)
        error('There are only new positions and is needed more info')
    end 
    newPosInaccUpdate = table();
    for k = 1:height(newPosInacc)
        % get ticker
        tick = newPosInacc.ISIN(k);
        % get tickers in fund
        posInacc = AccInvs(strcmp(tick, AccInvs.ISIN),:);
        % get a positions in Account as example
        posInaccEx = posInacc(end,:);
        % get position in Fund
        newPosInaccToUpdate = newPosInacc(k,:);
        % change inputs
        posInaccEx.Date = datetime(newPosInaccToUpdate.AccValuesD,...
            'format','MM/dd/yyyy');
        posInaccEx.loc_quantity = newPosInaccToUpdate.Total;
        posInaccEx.quantity = newPosInaccToUpdate.Loc_quantity;
        posInaccEx.shares = newPosInaccToUpdate.Shares;
        posInaccEx.price = newPosInaccToUpdate.Price;
        posInaccEx.fx = (posInaccEx.loc_quantity - newPosInaccToUpdate.fee)...
            / posInaccEx.quantity; %
        posInaccEx.fee = newPosInaccToUpdate.fee;
        posInaccEx.fee___ = 100*(posInaccEx.fee / posInaccEx.loc_quantity);
        newPosInaccUpdate = [newPosInaccUpdate;posInaccEx];
    end
    % writetable if is not a empty table
    if ~isempty(newPosInaccUpdate)
        accountUpdated = [account;newPosInaccUpdate];
        writetable(accountUpdated, xlsfileaccpath, 'filetype',...
            'spreadsheet', 'sheet', 1)
        newPosID = newPos(identifynewPosInAcc,:);
        productID = unique(StockndETFData.Producto(ismember(StockndETFData.ISIN,...
            newPosID{:,'ISIN'})), 'stable');
        disp([repmat('The following position is updated: ', ...
            size(productID, 1), 1),char(newPosID{:,...
            'ISIN'}),repmat('(',size(productID,1),1), char(productID),...
            repmat(')',size(productID,1),1)])
    end
else
    disp('There is no new positions in Stocks or ETFs')
end




