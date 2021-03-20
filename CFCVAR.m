function CFCVAR = CFCVAR(VAR, mu, sigma, alph, skew, kurt,dof,typeofdist)
n = typeofdist;
switch n
    case 1 % Normal Distribution
        CFCVAR = - mu + (sigma ./ alph) .* (VAR).* (1 + (skew./...
            (6.* (VAR).^3)) + kurt./24 .*((VAR).^4 - 2 .* (VAR).^2 - 1));
    case 2 % LogNormal Distribution
        CFCVAR = - mu + (sigma ./ alph) .* normcdf(log(VAR)).* (1 +...
            (skew./ (6.* (VAR).^3)) + kurt./24 .*((VAR).^4 - 2 .* (VAR).^2 - 1));
    case 3 % t Distribution
        pd = makedist('tLocationScale','mu',mu,'sigma',sigma,'nu',dof);
        y = pdf(pd,VAR);
        y = normcdf(y,mu,sigma);
        CFCVAR = - mu + (sigma ./ alph) .* y .* (1 + (skew./ (6.* (VAR).^3)) + ...
            kurt./24 .*((VAR).^4 - 2 .* (VAR).^2 - 1));
    case 4 % Pareto Distribution
        CFCVAR = - mu + (sigma ./ alph) .* gppdf(VAR,mu, sigma,0).* ...
            (1 + (skew./ (6.* (VAR).^3)) + kurt./24 .*((VAR).^4 - 2 .* (VAR).^2 - 1));
end
end