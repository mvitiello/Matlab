function ES = calcParamES(data, method, al, folder2save, varargin)

% function to calculate the parametric ES
% inputs: 
% data = returns of the historical serie
% method = model to calculate the ES, admit Normal, t and GHS distribution, 
% al = significance level
% folder2save = folder path to save info 
% varargin = indicate the period in the moving volatility to obtain the ES
% Output = ES,  expected shortfall with the model decided 

% Initial parameters 
% constant parameters  
if isempty(varargin)
    mu = mean(data);
    sigma = std(data);
else
    % moving parameters 
    moving = varargin{1};
    mu = movmean(data, moving);
    sigma = movstd(data, moving);   
end

% Normal ES 
if  isempty(method)|| strcmp(method, 'Normal')
    ES = -sigma ./(sqrt(2.*pi).*(al).*exp(erfinv(2.*al-1).^2));
% ES with t distribution     
elseif strcmp(method, 't')
    pt = mle(data, 'distribution', 'tlocationscale');
    nu = pt(3);
    if nu < 2
        nu = 2.65; % value from "literature"
    end
    ES = (-sigma.* (tpdf((tinv(al, nu)), nu)/(al)).* ...
        ((nu +(tinv(al, nu)).^2)./(nu - 1)));
    
% ES with GHS distribution  
elseif strcmp(method, 'GHS')
    
    ES =  +mu -2/pi*-sigma.*log(tan(pi/2*al))-2/(al*pi^2)*...
        -sigma.*1i*(dilog(-1i*tan(pi/2*al))-dilog(1i*tan(pi/2*al)));

% ES with Cornish Fisher Normal distribution  
elseif strcmp(method, 'CFn')
    skew = skewness(data);
    kurt = kurtosis(data, 0) - 3;
    VaR = VARCornishFisher(0,-sigma, skew, kurt, 1-al, 0, 1);
    ES = CFisherCVAR(VaR, 0, -sigma, 1-al, skew, kurt, 0, 1, folder2save);
    
 % ES with Cornish Fisher t distribution      
elseif strcmp(method, 'CFt')
    pt = mle(datat, 'distribution', 'tlocationscale');
    nu = pt(3);
    if nu < 2
        nu = 2.65; % value from "literature"
    end
    skew = skewness(data);
    kurt = kurtosis(data, 0) - 3;
    VaR = VARCornishFisher(0,-sigma, skew, kurt, 1-al, nu, 3);
    ES = CFisherCVAR(VaR, 0, -sigma, 1-al, skew, kurt, nu, 3, folder2save);
end

% avoid non-sense in CF method
if mean(ES) > 0 
    ES = - ES; 
end
ES(ES < -0.90) = NaN;
ES = fillmissing(ES,'pchip');  
end

