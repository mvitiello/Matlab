function mquantile = movingquantile(data, p, ndays)
if  ndays > 0
    for k = 1 : size(data, 1) - ndays
        datatquantile = data(k:(k + ndays));
        movingqunatilefunc(k) = quantile(datatquantile, p);
        mquantile = movingqunatilefunc';
    end
else
    mquantile = quantile(data, p);
end
end