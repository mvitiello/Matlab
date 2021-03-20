function newtableprices = addweekendprices(pricestable)
% this function add week end prices (prev value) or prices in no working days
% only needed in new positions or with missed data
% input = historical prices table
% output = historical prices table with prices in week end or in
% any  no working day

% check dates
startdate = pricestable.Date(1);
enddate = pricestable.Date(end);
adates = linspace(startdate, enddate, days(enddate - startdate)+1)';
dates2fill = ismember(adates, pricestable.Date);
% create and insert prices
prices2fill = nan(sum(dates2fill == 0), width(pricestable) - 1);
weekendprices = timetable(adates(~dates2fill), prices2fill(:,1),...
    prices2fill(:,2), prices2fill(:,3), prices2fill(:,4));
weekendprices.Properties.VariableNames = pricestable.Properties.VariableNames(2:5);
% newtable
newtableprices = fillmissing(sortrows([table2timetable(pricestable); ...
    weekendprices], 'Date'), 'previous');
end