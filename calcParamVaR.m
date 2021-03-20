function VaR = calcParamVaR(data, method, al, varargin)

% function to calculate the parametric VaR
% inputs: 
% data = returns of the historical serie
% method = model to calculate the VaR, admit Normal, t and GHS distribution, 
% al = significance level
% varargin = indicate the period in the moving volatility to obtain the VaR
% Output = VaR,  Value at risk with the model decided 

% Initial parameters 
% constant parameters  
if isempty(varargin)
    mu = mean(data);
    sigma = std(data);
    skew = skewness(data);
    kurt = kurtosis(data, 0) - 3;
else
    % moving parameters 
    moving = varargin{1};
    mu = movmean(data, moving);
    sigma = movstd(data, moving);
    skew = movskew(data, moving);
    kurt = movkurt(data, moving);
end

% Normal VaR
if  isempty(method)|| strcmp(method, 'Normal')
    VaR = - mu + sigma.*norminv(al);
    
% t student VaR    
elseif strcmp(method, 't')
    pt = mle(data, 'distribution', 'tlocationscale');
    nu = pt(3);
    if nu < 2
        nu = 2.65; % value from "literature"
    end
    
    VaR = tinv(al, nu)* sigma .* ...
        sqrt((nu - 2)/nu);

% VaR with Generalised Hyperbolic Secant or NEF-GHS distribution    
elseif strcmp(method, 'GHS') 
    
    VaR = (2/pi)*log(tan((pi/2)*(al)));

% VaR with Cornish-Fisher - Normal distribution    
elseif strcmp(method, 'CFn')
    
    VaR = VARCornishFisher(0,-sigma, skew, kurt, 1-al, 0, 1);

%VaR with Cornish-Fisher - t distribution      
elseif strcmp(method, 'CFt')
    pt = mle(datat, 'distribution', 'tlocationscale');
    nu = pt(3);
    if nu < 2
        nu = 2.65; % value from "literature"
    end
    
    VaR = VARCornishFisher(0,-sigma, skew, kurt, 1-al, nu, 2);
    
end
% sometimes with Cornish Fisher a VaR higher than 0.90 is obtained.
% In the "real world" we need to put a floor in the VAR if the value is
% higher than 90% 
VaR(VaR < -0.90) = NaN;
VaR = fillmissing(VaR,'pchip');  
end

