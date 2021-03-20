function TableVarES = CFPVARNT(ret, DateReturns, ti, startDate,...
    endDate, tickname, ftick, model, rdix, path)

% TableVarES =  CFPVARESNT(ret, DateReturns, ti, startDate,...
%     endDate, tickname, ftick, rdix, path)
% calculate the parametric dynamic VaR and ES (i.e, using a time
% horizon window) the code calculate Normal and T distribution
% parametric VaR with 95 and 99 confidence level (CL).
% Using the Cornish Fisher Expansion. 
% the inputs are:
% ret, returns variable
% tickname,ftick, are the ticker of the variable, e.g.,'^GSPC',
% and ftick is the list with the tckers and with long names
% of the ticker,  e.g., ftick = {'^GSPC', 'SP500';...
% '^IBEX', 'IBEX35';...
% '^FCHI', 'CAC40';...
% '^FTSE', 'FTSE100'};
% DateReturns, dates downloaded from yahoo
% startDate and endDate are the point in time date that we
% want to start and finish with in the analysis.
% ti, a year, usually one year more than the startDate. This
% variable is used to calculate the horizon window and
% is assumed a 250 days of window
% model is the name of the model of the simulated returns, e.g.,
% 'BrowG' refers to the geometric brownian motion. This variable is
% only needed in the name of the plot (only to identify the model)
% rdix, indicates if ret is a simulated returns (rdix = 1) or
% real returns (rdix = 0).
% path is the folder path where the files are saved
% the outputs are:
% NV95, Normal parametric dynamic Var with 95 CL
% TV95, t parametric dynamic VaR with 95 CL
% all variables in a table.
% Try:
% ftick = {'^GSPC', 'SP500';...
% '^IBEX', 'IBEX35';...
% '^FCHI', 'CAC40';...
% '^FTSE', 'FTSE100'};
% inidate = '01012014';
% enddate = '01102018';
% DateReturns = linspace(datetime(inidate,'InputFormat', 'ddMMyyyy'),...
%     datetime(enddate,'InputFormat', 'ddMMyyyy'), ...
%     datenum(datestr(datetime(enddate,'InputFormat', 'ddMMyyyy')))-...
%     datenum(datestr(datetime(inidate,'InputFormat', 'ddMMyyyy'))))';
% tableVRES =  CFPVARESNT(randn(1000,1), DateReturns, 2015, inidate,...
%     enddate, '^FCHI', ftick,'BrowG', 0, ...
%     'C:\Users\Miguel\Desktop\Matlab\');

% check folder or create it 
if ~exist(path, 'dir')
    mkdir(path)
end 

% check the list
idx = ismember(ftick, tickname);
rtick = ftick(idx(:, 1), :);
% Size and window
sampleSize = length(ret);
testWindowStart = find(year(DateReturns)== ti, 1);
if testWindowStart <= 250
    testWindowStart = 251;
end
testWindow = testWindowStart : sampleSize;
% is assumed a window of 250 days
estimationWindowSize = 250;
% confidence level
pVaR = 0.95;
% first distribution - Normal distribution
dist1 = mle(ret, 'distribution', 'norm');
nmu = dist1(1);
% second distribution - t location scale
dist2 = mle(ret, 'distribution', 'tlocationscale');
n1 = dist2(3);
% if the degree of freedom (dof) is less than 2 is a problem
% with the variance scale (cannot be scale it by dof because
% must be higher than 2 - see formula below)
if n1 <= 2
    n1 = 2.01;
end
mu = dist2(1);
% initial variables
NV95  = zeros(size(testWindow, 2), 1);
TV95  = zeros(size(testWindow, 2), 1);


% calculate a "dynamic" Var, i.e., with time horizon window
for t = testWindow
    % get the window
    i = t - testWindowStart + 1;
    estimationWindow = t-estimationWindowSize:t-1;
    % we obtain the standard desviation, skew and kurtosis
    % in the time horizon window
    sigma = std(ret(estimationWindow));
    skewnd = skewness(ret(estimationWindow), 0);
    kurtnd = kurtosis(ret(estimationWindow), 0) - 3;
    % use parametric Normal VaR with 95 and 99 % of confidence level
    % using Cornish Fisher (CF)
    NV95(i)  = VARCornishFisher(-nmu, sigma, skewnd,kurtnd, ...
        pVaR(1), 0, 1);
end
for t = testWindow
    % get the window
    i = t - testWindowStart + 1;
    estimationWindow = t-estimationWindowSize : t - 1;
    % we obtain the standard desviation in the time horizon window
    sigmat = std(ret(estimationWindow));
    skewtd = skewness(ret(estimationWindow));
    kurttd = kurtosis(ret(estimationWindow), 0) - 3;
    % use parametric t VaR with 95 and 99 % of confidence level
    % using a CF method
    TV95(i)  = VARCornishFisher(-mu, sigmat, skewtd,kurttd,...
    pVaR(1), n1, 3);
end
% create the dates (abcises)
startDate2 = datenum(DateReturns(find(year(DateReturns) == ti, 1)));
enddate = datenum(DateReturns(size(DateReturns, 1)));
xData2 = linspace(startDate2, enddate, size(NV95, 1));
m = size(ret, 1)- size(NV95, 1);
ret2 = ret(m + 1 : sampleSize);
% plot
subplot(1, 2, 1)
plot(xData2, [NV95 TV95]);
xlim([startDate2 enddate])
legend('VN95', 't95')
xlabel('Dates'), ylabel('VAR'),
if rdix == 1
    title(['Parametric CF Normal and t VaR, model CF-',...
        char(string(model))])
else
    title('Parametric CF Normal and t VaR')
end
datetick('x', 'yyyy', 'keeplimits')
subplot(1, 2, 2)
plot(xData2, [ret2 NV95 TV95]); ...
    xlim([startDate2 enddate]),
datetick('x', 'yyyy', 'keeplimits')
if rdix == 1
title(['Parametric CF VaR and ES from ', char(string(rtick(1, 2))), ...
    ' , CF-',char(string(model)), ' Model']),
legend('Simulated returns', 'VN95','t95')
else
    title(['Parametric CF VaR and ES from ',...
        char(string(rtick(1, 2)))])
legend('Real returns', 'VN95','t95')
end 
xlabel('Dates'), ylabel('returns and VAR')
Names = {'CFVN95', 'CFVT95'};
TableVarES = table(NV95,TV95, 'VariableNames', Names);

% write the table
if rdix == 1
    writetable(TableVarES, [path, char(string(rtick(1, 2))), ...
        '_VARESRealSimul_',char(string(model)),'_',...
        datestr(datetime(startDate,'InputFormat',...
        'ddMMyyyy'), 'ddmmyyyy'), ' to ',...
        datestr(datetime(endDate,...
        'InputFormat', 'ddMMyyyy'), 'ddmmyyyy'), '.xls'], 'Sheet', 9)
else
    writetable(TableVarES, [path, char(string(rtick(1, 2))),...
        '_VARESRealSimul_',char(string(model)),'_',...
        datestr(datetime(startDate,'InputFormat',...
        'ddMMyyyy'), 'ddmmyyyy'), ' to ',...
        datestr(datetime(endDate,...
        'InputFormat', 'ddMMyyyy'), 'ddmmyyyy'), '.xls'], 'Sheet', 8)
end

% open Activex server
e = actxserver('Excel.Application');
% open file (enter full path!)
ewb = e.Workbooks.Open([path, ...
    char(string(rtick(1, 2))),'_VARESRealSimul_',...
    char(string(model)),'_', datestr(datetime(startDate,...
    'InputFormat', 'ddMMyyyy'), 'ddmmyyyy'), ' to ',...
    datestr(datetime(endDate,...
    'InputFormat', 'ddMMyyyy'), 'ddmmyyyy'), '.xls']);
% rename 3st sheet
if rdix == 1
    ewb.Worksheets.Item(9).Name = [char(string(rtick(1, 2))), ...
        '_CFPVARSimul',char(string(model))] ;
else
    % rename 2st sheet
    ewb.Worksheets.Item(8).Name = [char(string(rtick(1, 2))),...
        '_CFPVARReal',char(string(model))] ;
end
% save the file, close and quit
ewb.Save
ewb.Close(false)
e.Quit
end
