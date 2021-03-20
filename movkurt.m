function mkurt = movkurt(data, varargin)

% function that calculate moving kurtosis using a window determined by the
% variable varargin. If varargin is empty a accumulated kurtosis moving 
% one day is obtained 
% input:
% data = variable to calculate the metric
% varargin = window to calculate the metric
% output:
% mkurt = moving kurtosis

if isempty(varargin)
for k = 1 : length(data)
    % moving kurtosis by accumulating data, converges to kurtosis(data)
    datat = data(1:k);
    kurtcumulativefunc(k) = kurtosis(datat, 0) - 3;
    mkurt = kurtcumulativefunc';
    mkurt(isnan(mkurt)) = 0;
end
else
    % moving skewness by sliding a window of length s along data
    moving = varargin{1};
    for k = 1 : length(data)- moving
        kurtcumulativefunc(k)= kurtosis(data(k:(k+moving)), 0) - 3;
        mkurt = kurtcumulativefunc';
        mkurt = [zeros(moving,1); kurtcumulativefunc'];
        mkurt(isnan(mkurt)) = 0;
    end
end
end

