function plotTimeSeries(Time,Serie, titleprov, legendprov, xlabelprov, ...
    ylabelprov)
% first date
startDate = datenum(Time(1));
% last date
endDate = datenum(Time(end));
% axes with dates
xData = linspace(startDate, endDate, length(Time));
%plot
plot(xData, Serie)
xlim([startDate endDate])
datetick('x', 'ddmmyyyy', 'keeplimits');
% title and legend
title(titleprov)
legend(legendprov, 'Location','best')
xlabel(xlabelprov)
ylabel(ylabelprov)
end

