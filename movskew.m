function mskew = movskew(data, varargin)

% function that calculate moving skew using a window determined by the
% variable varargin. If varargin is empty a accumulated skew moving one day
% is obtained
% input:
% data = variable to calculate the metric
% varargin = window to calculate the metric
% output:
% mskew = moving skew

if isempty(varargin)
    % moving skewness by accumulating data, converges to skewness(data)
    for k = 1 : length(data)
        datat = data(1:k);
        skewcumulativefunc(k) = skewness(datat);
        mskew = skewcumulativefunc';
        mskew(isnan(mskew)) = 0;
    end
else
    % moving skewness by sliding a window of length s along data
    s = varargin{1};
    for k = 1 : length(data)- s
        skewcumulativefunc(k)= skewness(data(k:(k+s)));
        mskew = [zeros(s,1); skewcumulativefunc'];
        mskew(isnan(mskew)) = 0;
    end
end
end
