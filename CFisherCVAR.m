function CFCVAR = CFisherCVAR(VAR, mu, sigma, alph, skew, kurt,...
    dof, typeofdist, folder)

% This function calculate the ES with the Cornish-Fisher Expansion
% we have different distributions and we use the same 'methodology'
% used in the ES (most conservative!)
% we change the distribution with the input typeofdist
% this code also calculate the CF parameters to simulate with the
% CF distribution (see Maillard, 2013).  
% In this code is included the range where the expansion is valid. 
% the INPUTS are:
% VAR, Value at risk using CF expansion
% mu, location parameter
% sigma, scale parameter
% skew, skewness parameter
% kurt, kurtosis parameter
% alph, percentile
% dof, degree of freedom, for t-sudent
% typeofdist, 1 for Normal, 2 for Lognormal, 3 for t distribution,
% 4 for Weibull
% folder = folder to save the parameters
% The OUTPUT is the ES with the model selected 


% type of distribution 
n = typeofdist;

% we need "estimated" parameters of the distribution
if abs(skew) >= 1.5 | abs(kurt)>= 2 % "rule of thumb"
    f = @(x)calculateCFKS(x,skew,kurt);
    x0 = [0,0];
    s = fsolve(f,x0);
    k1 = s(1);
    k2 = s(2);
    % we have the following restriction for this transformation
    if (6*(sqrt(2)-1) >= abs(k1) || abs(k1) >= 6*(sqrt(2)+1)) && ...
            ((1+k1^2 - sqrt(k1^4-6*k1^2 +1))/6 >= abs(k2) ||...
            abs(k2) >= (1+k1^2 + sqrt(k1^4-6*k1^2 +1))/6)
        disp('Cornish Fisher inequality is satisfied')
    else
        warning('Cornish Fisher inequality is not satisfied')
    end
else
    k1 = skew;
    k2 = kurt;
end

% parameters to simulate
K = k1/24;
S = k2/6;
a1 = K - 2*S^2;
a2 = -S/(a1);
p = ((1 - 3*K + 5*S^2)/a1) - 1/3*(S^2/(a1)^2);
q = a2 - 1/3*(S*((1 - 3*K + 5*S^2)/a1^2)) - 2/27*(S/a1)^3;
z = a2/3 + ((-q + sqrt(q^2+4/27*p^3))/2)^(1/3) + ...
    ((-q + sqrt(q^2*4/27*p^3))/2)^(1/3);

% write table
tableParameters = table(real(k1),real(k2), real(z));
tableParameters.Properties.VariableNames = {'K', 'S', 'z'};
tableParameters = rows2vars(tableParameters);
tableParameters.Properties.VariableNames(2) = {'Params'};
writetable(tableParameters, char(strcat(folder, '\',...
    'CornishFisherParam', '.xlsx')), 'filetype', 'spreadsheet', ...
    'sheet', 1)

% calculate ES 
switch n
    case 1 % Normal Distribution
        
        CFCVARnk = + mu + (sigma ./ (1-alph)) .* normpdf(((VAR + mu)...
            ./ sigma),0,1).* (1 + ((k1./6) .* (((VAR + mu)./sigma).^3))+...
            k2./24 .* ((((VAR + mu)./sigma)).^4 - 2 .* (((VAR + mu)...
            ./sigma)).^2 - 1));
        
        CFCVARk = + mu + (sigma ./ (1-alph)) .* normpdf(((VAR + mu)...
            ./ sigma),0,1).* (1 + ((skew./6) .* (((VAR + mu)./sigma).^3))+...
            kurt./24 .* ((((VAR + mu)./sigma)).^4 - 2 .* (((VAR + mu)...
            ./sigma)).^2 - 1));
        
    case 2 % LogNormal Distribution
        
        CFCVARnk = 1-(exp(- mu + (0.5* sigma^2 ))./ (1-alph)) .*...
            normpdf(((log(VAR + mu)/sigma)),0,1).* (1 + ((k1./6).*...
            normcdf(log((VAR + mu + sigma^2)/sigma)).^3 +...
            k2./24 .* normcdf(log((VAR + mu +sigma^2)/sigma)).^4 ...
            - 2 .* normcdf(log((VAR + mu +sigma^2)./sigma)).^2 - 1));
        
        CFCVARk = 1-(exp(- mu + (0.5* sigma^2 ))./ (1-alph)) .*...
            normpdf(((log(VAR + mu)/sigma)),0,1).* (1 + ((skew./6).*...
            normcdf(log((VAR + mu + sigma^2)/sigma)).^3 +...
            kurt./24 .* normcdf(log((VAR + mu +sigma^2)/sigma)).^4 ...
            - 2 .* normcdf(log((VAR + mu +sigma^2)./sigma)).^2 - 1));
        
    case 3 % t Distribution
        pd = makedist('tLocationScale','mu',max(mu),'sigma',max(sigma),...
            'nu',dof);
        y = cdf(pd,VAR);
        
        CFCVARnk = + mu + (sigma ./ (1-alph))... %.*(sqrt((dof - 2)/dof))...
            .* y .* (1 + ((k1./6).* (tpdf(tinv(VAR, dof),dof).^3)) + ...
            k2./24 .*(tpdf(tinv(VAR, dof),dof).^4 -...
            2 .* (tpdf(tinv(VAR, dof),dof).^2 - 1)));
        
        CFCVARk = + mu + (sigma ./ (1-alph))... %.*(sqrt((dof - 2)/dof))...
            .* y .* (1 + ((skew./6).* (tpdf(tinv(VAR, dof),dof).^3)) + ...
            kurt./24 .*(tpdf(tinv(VAR, dof),dof).^4 -...
            2 .* (tpdf(tinv(VAR, dof),dof).^2 - 1)));
    case 4 % Weibull Distribution
        % We use the same 'metodology' for ES to be conservative!
        CFCVARnk =  (mu./(1- alph)) .* wblcdf(VAR,mu,sigma).*...
            (1 + ((k1./6).* ((gammainc(1 + 1./sigma,-log(alph),...
            'upper')) .^3)) + k2./24 .*((gammainc(1 + 1./sigma,...
            -log(alph),'upper')).^4 - 2 .* (gammainc(1 + 1./sigma,...
            -log(alph),'upper')).^2 - 1));
        
        CFCVARk =  (mu./(1- alph)) .* wblcdf(VAR,mu,sigma).*...
            (1 + ((skew./6).* ((gammainc(1 + 1./sigma,-log(alph),...
            'upper')) .^3)) + kurt./24 .*((gammainc(1 + 1./sigma,...
            -log(alph),'upper')).^4 - 2 .* (gammainc(1 + 1./sigma,...
            -log(alph),'upper')).^2 - 1));
        
end
% select the highest CVAR
if mean(CFCVARk) > mean(CFCVARnk)
    CFCVAR = CFCVARk;
else
    CFCVAR = CFCVARnk;
end

end

% internal function 
function F = calculateCFKS(x, skew, kurt)
F(1) = ((x(1) - 76/216 .*x(1).^3 + 85/1296.*x(1).^5 +...
    1/4 .*x(1).*x(2) - 13/144 .* x(1)^.3.*x(2) + ...
    1/32 .* x(2).^2.*x(1))./ ((1 + 1/96.*x(2).^2 + ...
    25/1296.*x(1).^4 - 1/36.*x(2).*x(1).^2).^3/2))- skew;
F(2) = (( 3 + x(2) +7/16 .*x(2).^2 + 3/32.*x(2).^3 +...
    31/3072 .*x(2).^4 - 7/216 .* x(1).^4 - 25/486 .*x(1).^6 +...
    21605/559872 .* x(1).^8 - 7/12.*x(2).*x(1).^2 +...
    113/452 .*x(2).*x(1).^4 - 5155/46656 .*x(2).*x(1).^6 - ...
    7/24 .*x(2).^2.*x(1).^2 + 2455/20736.*x(2).^2.*x(1).^4 -...
    65/1152 .*x(2)^3.*x(1).^2 )./((1 + 1/96.*x(2).^2 + ...
    25/1296*x(1).^4 - 1/36.*x(2).*x(1).^2).^2)) - kurt;
end