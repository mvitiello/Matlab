function tableWithBacktestingES = ...
    bcktstingexepshrtfall2pnl(ret, prices, ...
    Var, tailvar, model, p, folder)

% this function obtains a table with the results
% from a ES's backtesting:
% we calculate the test2 following the Acerbi Paper:
% https://www.msci.com/documents/10199/22aa9922-f874-4060-b77a-0f0e267a489b
% the inputs are:
% ret = the returns
% prices = historical prices 
% Var = value at risk 
% tailvar = Expected Shortfall(ES) value
% model =  model used to calculate the VaR: normal, t, cornish-fisher,
% variance gamma....
% p = significance level,
% folder = folder path to get info and find the ticker  
% the output is table with a backtesting table using a unconditional test
% and conditional test 
% with normal and t student distributon both test are calculated
% with cf and vg distribution only unconditional test is calculated. 

% Do not take account weekends
retv = ret(ret~=0);
if ~(size(Var,1) == 1)
    Var = Var(ret~=0, :);
    tailvar = tailvar(ret~=0, :);
end
% fill 0 values with the mean 
if any(any(tailvar == 0))
    tailvar = standardizeMissing(tailvar,0);
    tailvar = fillmissing(tailvar, 'movmean', 5);
end
% starting index
n = size(retv,1);
m = size(tailvar,2);
q = size(tailvar,1);
I = zeros(n,m);

% this is in case the ES is a single value
if m == 1
    tailvar = tailvar .* ones(1, m);
end

% change the dimensions 
if q == 1
    tailvar = tailvar .* ones(n,m);
    Var = Var .* ones(n,m);
end

% obtaining exceptions
for i = 1:size(I, 1)
    for j = 1:m
        I(i,j) = retv(i,:) < Var(1,j);
    end
end

I2 = (retv .* I)./ -tailvar;

% Index to test2
test2 = 1 + sum(I2) /(p * n);


tf = table();

% unconditional and conditional backtesting test 
for k = 1:numel(model)
    if strcmp(model(k), 't')|| strcmp(model(k), 'normal')
        Mu = movmean(retv, 5);
        Sigma = movstd(retv, 5);
        % get the dof modelling the ret return as a t student 
        param = fitdist(ret, 'tlocationscale');
        % a entire number and higher than 3 of df is needed 
        DoF = ceil(param.nu); 
        if DoF < 3
            DoF = ceil(3 + param.nu);
        end
        % IDs
        if strcmp(model(k), 't')
            IDs = "t(dof) 95%";
            IDs = strrep(IDs,"dof",num2str(DoF));
        else
            IDs = "N95%";
        end
        % condfidence level
        VaRLevel = 1 - p ;
        % Ticker
        findticker = split(folder, '\');
        ticker = findticker(end);
        % with matlab internal functions
        obj = esbacktestbysim(retv,abs(-Var(:,k)),abs(-tailvar(:,k)),string(model(k)),...
            'DegreesOfFreedom',DoF,...
            'Location',Mu,...
            'Scale',Sigma,...
            'PortfolioID',string(ticker),...
            'VaRID',IDs,...
            'VaRLevel',VaRLevel);
        % conditional test 
        t1 = conditional(obj);
        t1.TypeTest = strcat('Conditional', '-', t1.VaRTest);
        t1.Conditional = [];
        t1.VaRTest = [];
        t1.VaRTestResult = [];
        t1.VaRTestPValue = [];
        t1.Properties.VariableNames('ConditionalOnly') = {'TestDesicion'};
        % unconditonal test 
        % matlab internal functions simulate the path with the normal
        % distribution (geometric brownian motion) and t student 
        t2 = unconditional(obj);
        t2.Properties.VariableNames(4) = {'TestDesicion'};
        t2.TypeTest = 'Unconditional' ;
        tf = [tf; [t1;t2]];
    else
        % Unconditional test with variance gamma model and Cornins-Fisher
        % distribution 
        % Ticker
        findticker = split(folder, '\');
        ticker = repmat(findticker(end),m,1);
        % type test
        TypeTest = repmat({'Unconditional'}, m,1);
        % VaR level
        VaRLevel = (1 - p) .*ones(m,1);
        TestLevel = VaRLevel;
        Observations = n .*ones(m,1);
        TestStatistic = test2;
        % Calculate p value with simulation
        Scenarios = 1000;
        simul = zeros(Scenarios, size(model,2));
        if strcmp(model(k), 'vg')
            ID(k,:) = {'VG95%'};
            % obtain parameteres 
            paramsVG = estimateVGParametersGen(prices, folder);
            % simulate a path with variance gamma 
            simulp = price2ret(simulaVG(prices(1), paramsVG.Params.params(1),...
                paramsVG.Params.params(3), paramsVG.Params.params(2), ...
                paramsVG.Params.params(4), Scenarios, 1, 0));
        elseif strcmp(model(k), 'cf')
            ID(k,:) = {'CF95%'};
            % obtain CF parameters, obtained in CFisherCVAR.m
            paramsCF = readtable([char(folder), '\', 'CornishFisherParam.xlsx']);
            z = paramsCF.Params(3);
            K = paramsCF.Params(1);
            S = paramsCF.Params(2);
            % simulate a path with CF distribution: geometric  
            % brownian motion path adjusted with kurtosis and skew  
            mu = mean(ret(ret~=0));
            s = std(ret(ret~=0));
            dt = 1;
            drift = (mu - 0.5 * s ^ 2)* dt * ones(Scenarios, 1);
            difussion = (s*sqrt(dt))*randn(Scenarios, 1);
            brownianProcess  = exp(drift + difussion);
            spotPrice   = prices(1) * ones(1, 1);
            simulation  = [spotPrice; brownianProcess];
            stockpath = cumprod(simulation);
            retsimul = price2ret(stockpath);
            % adjust with curtosis and skew
            simulp = retsimul ./ (z^2*(K/8 - S^2/6)+S*z/3+1-K/8+5*S^2/36);
        end
        simul(:,k) = simulp;
        % ratio to calculate pvalue
        pvalue = sum(simul<test2)./Scenarios;
        % critical value and test Desicion
        SimTestStats = sort(simul);
        CutOff = Scenarios*(1-TestLevel);
        Ind = max(ceil(CutOff),1);
        CriticalValue = SimTestStats(Ind, k);
        testDesicion = abs(TestStatistic') <  abs(CriticalValue);
        testDesicion = categorical(testDesicion,[true;false],{'accept';'reject'});
    end
end

% unconditional test table for vg and cf distributions 
if strcmp(model(k), 'vg')|| strcmp(model(k), 'cf')
 Names = {'PortfolioID','VaRID','VaRLevel','TestDesicion','PValue',...
            'TestStatistic','CriticalValue','Observations','Scenarios',...
            'TestLevel', 'TypeTest'};
        t3 = table(ticker, ...
            ID, VaRLevel, testDesicion , pvalue', TestStatistic',...
            CriticalValue, Observations, Scenarios.*ones(m,1),TestLevel, ...
            TypeTest, 'VariableNames', Names);
        t3b = table2cell(t3);
        t3b = t3b(~any(cellfun('isempty', t3b), 2), :);
        t3b = cell2table(t3b);
        t3b.Properties.VariableNames = t3.Properties.VariableNames;
        t3 = t3b;
end

if ~exist('t3', 'var')
    t3 = table();
end

% join tables 
tableWithBacktestingES = [tf;t3];
%tableWithBacktestingES = rmmissing(tableWithBacktestingES);
end