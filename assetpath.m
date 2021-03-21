function stockpath = assetpath(S0, m, s, t, M, N,...
    method, simMethod, plotg)
% The function assetpath() simulate N price paths with M points (days)
% With a geometric brownian motion. The following INPUTS are needed:
% m   = process's drift
% s   = process's diffusion (volatility), constant
% S0  = initial price of the path
% t   = period of time, usually daily, i.e t = 1
% M   = time horizon: points, days, usually 252 days or more
% N   = numbers of paths
% method, if the method is 'direct' the path os obtaned with a written code
% if the method is 'GBM', is used  a MATLAB GBM object
% simMethod, only works with GBM method is used, the simMethod are
% 'byEuler' and 'bySolution'
% plotg, to plot the graphs
% OUTPUT
% stockpath = prices paths (Matrix M,N)
% try
%  S = assetpath(100, 0.004, 0.08, 1, 4*252, 2, 'direct', '', 1)
%  S = assetpath(100, 0.004, 0.08, 1, 4*252, 2, 'GBM', 'default',1)
%  S = assetpath(100, 0.004, 0.08, 1, 4*252, 2, 'GBM', 'ByEuler',1)
%  S = assetpath(100, 0.004, 0.08, 1, 4*252, 2, 'GBM', 'BySolution',1)

% First, the brownian motion process is created
% The brownian motion is esentialy constant difussion
% and a random component

if contains(method, 'direct')
    % time step
    dt   = t; % t/M % determinist part
    drift = (m - 0.5 * s ^ 2)* dt * ones(M, N);
    % stochastic part
    difussion = (s*sqrt(dt))*randn(M, N);
    % both part are added and exponentialized
    brownianProcess  = exp(drift + difussion);
    % Next, the initial price is added
    % The spot price is needed
    spotPrice   = S0 * ones(1, N);
    simulation  = [spotPrice; brownianProcess];
    stockpath = cumprod(simulation);
elseif contains(method, 'GBM')
    % simulation using a GBM object
    GBMobj = gbm(m, s, 'StartState', S0);
    rng(142857,'twister')
    if contains(simMethod, 'default')
        % simulation using simulate 
        stockpath = simulate(GBMobj, M, 'DeltaTime', t, ...
            'nTrials', N);
    elseif contains(simMethod, 'ByEuler')
        % simByEuler uses the Euler approach to approximate
        % continuous-time stochastic processes.
        stockpath = simByEuler(GBMobj, M, 'DeltaTime', t, ...
            'nTrials', N);
    elseif contains(simMethod, 'BySolution')
        % The simBySolution function simulates the state vector 
        % using an approximation of the closed-form solution of diagonal 
        % drift Hull-White/Vasicek Gaussian Diffusion (HWV) processes
        stockpath = simBySolution(GBMobj, M, 'DeltaTime', t, ...
            'nTrials', N);
    end
end
if plotg
    p1 = figure('visible', 'on', 'units', 'normalized', ...
        'outerposition', [0 0 1 1]);
    if contains(method, 'direct')
        plot(stockpath)
        xlim([1 M])
        title('Simulations with Geometric Brownian Motion')
    else
        for k = 2:N
            plot(stockpath(:,:,1))
            hold on
            plot (stockpath(:,:,k))
            xlim([1 M])
            title(['Simulations with Geometric Brownian Motion',...
                ' with method ', simMethod])
        end
    end
    xlabel('days'), ylabel('Prices')
    saveas(p1, 'Simulations with Geometric Brownian Motion', 'jpg')
end
end