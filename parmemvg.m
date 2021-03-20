% Approximation of VG pdf
function [a, m, sig, m0] = parmemvg(x, m0_in, m_in, beta, ...
    sig_in, a_in, u, w, passi) 
a = [a_in zeros(1, passi)];
m = [m_in zeros(1, passi)];
m0 = [m0_in zeros(1, passi)];
sig = [sig_in zeros(1, passi)];
phi = [psi(a_in) zeros(1, passi)];
T = length(x);
for h = 2 : passi
    phi(1, h) = sum(log(u)' * postdistrib(x, m0(1, h-1), m(1, h-1), ...
        beta, sig(1, h-1), a(1,h-1), u, w)) / T;
    a(1, h) = real(invpsi(phi(1, h)));
    m0(1, h) = (((u.^-1)' * postdistrib(x, m0(1, h-1), m(1, h-1),...
        beta, sig(1, h-1), a(1, h-1), u, w)* x) - ...
        T * sum(x)/sum(u' * postdistrib(x, m0(1, h-1), m(1, h-1), ...
        beta, sig(1, h-1), a(1, h-1), u, w)))/...
        (sum((u.^-1)' * postdistrib(x, m0(1, h-1), m(1, h-1), ...
        beta, sig(1, h-1), a(1, h-1), u, w))- ...
        T^2/sum(u' * postdistrib(x, m0(1, h-1), m(1, h-1),...
        beta, sig(1, h-1), a(1, h-1), u, w)));
    m(1, h) = (sum(x) - T*m0(1, h)) * beta/sum(u' * postdistrib(x, ...
        m0(1, h-1), m(1, h-1), beta, sig(1, h-1), a(1, h-1), u, w));
    sig(1, h) = sqrt(1/T*sum(sum(((x * ones(1, length(u)))' - ...
        (u * ones(1, length(x)))* m(1, h)/beta - m0(1,h)).^2 .* ...
        postdistrib(x, m0(1, h-1), m(1, h-1), beta, ...
        sig(1, h-1), a(1, h-1), u, w))));
end
m = m(passi);
sig = sig(passi);
a = a(passi);
m0 = m0(passi);
end

function Y = invpsi(X)
% Y = INVPSI(X)
% Inverse digamma (psi) function.  The digamma function is the
% derivative of the log gamma function.  This calculates the value
% Y > 0 for a value X such that digamma(Y) = X.
% This algorithm is from Paul Fackler:
% http://www4.ncsu.edu/~pfackler/
  L = 1;
  Y = exp(X);
 while L > 10e-8 
    Y = Y + L*sign(X-psi(real(Y)));
    L = L / 2;
 end
end

function y = postdistrib(x, d, m, beta, sig, a, u, w)
pesi = approxgammadens(a, u, w);
num = zeros(length(u), length(x));
y = zeros(length(u), length(x));
for i=1:length(x)
    num(:,i) = normpdf(x(i), d + m*u/beta, ...
        sig*sqrt(u/beta)).*pesi;
    y(:,i) = num(:,i)/sum(num(:,i));
end
end

function y = approxgammadens(a, u, w)
 y_num = u.^(a-1).*w;
 den = sum(y_num);
 y = y_num/den;
end
