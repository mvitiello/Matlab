function VARCF = VARCornishFisher(mu,sigma, skew, kurt, alph,...
    dof, typeofdist)
% this function calculate the VaR with the Cornish-Fisher Expansion
% we have different distributions and we only change the percentile
% (most conservative!)
% we change the distribution with the input typeofdist
% the inputs are:
% mu, location parameter
% sigma, scale parameter
% skew, skewness parameter
% kurt, kurtosis parameter
% alph, percentile
% dof, degree of freedom, for t-sudent
% typeofdist, 1 for Normal, 2 for Lognormal, 3 for t distribution,
% 4 for Weibull Distribution
% we use the Hermite polynomials

n = typeofdist;

switch n
    case 1 % Normal Distribution
        p = norminv(alph);
    case 2 % LogNormal Distribution
        p = logninv(alph,mu,sigma);
    case 3 % t distribution
        p = tinv(alph,dof);
    case 4 % Weibull distribution
        p = wblinv(alph,mu,sigma);
end

% using the Hermite Polynomials,
% Note: the matlab function HermiteH obtain the
% "physicists' Hermite polynomials" (Hn(n,x)) and we need the
% "probabilists' Hermite polynomials" (He(n,x)), so we have
% to rescaling: He(n,x) = 2^-n/2* Hn(n,x/sqrt(2)) where n is the
% order of the polinomy

% we need "estimated" parameters of the distribution
if abs(skew) >= 1.5 | abs(kurt)>= 2 % "rule of thumb"
    f = @(x)calculateCFKS(x,skew,kurt);
    x0 = [0,0];
    s = fsolve(f,x0);
    k1 = s(1);
    k2 = s(2);
    % we have the following restriction for this transformation
    if   (6*(sqrt(2)-1) >= abs(k1) || abs(k1) >= 6*(sqrt(2)+1)) && ...
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
h1 = (2^(-2/2)*hermiteH(2, p/sqrt(2)))/6;
h2 = (2^(-3/2)*hermiteH(3, p/sqrt(2)))/24;
h11 = - (2*(2^(-3/2)*hermiteH(3, p/sqrt(2))) + ...
    (2^(-1/2)*hermiteH(1, p/sqrt(2))))/36;

perc = p + (k1 .* h1) + (k2 .* h2 + (k1.^2) .* h11);
VARCF =  mu + sigma .* perc;

end

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