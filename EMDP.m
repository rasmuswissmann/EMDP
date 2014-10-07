function [res] = EMDP()

format long
warning off all

res = [];

%% Example: European Arithmetic Basket in a Log-Normal Stock Model

% N = 10;
% T = 1;
% sigma = 0.2;
% rho = AssetModel.buildSimpleConstantCorrelationMatrix(N,0.7);
% assetModel = LogConstCoeffModel(N,sigma,rho,0);
% % assetModel = assetModel.setRanges([-Inf,Inf,log(10),log(1000)]);
% 
% omega = ones(N,1)/N;
% K = 100;
% derivative = ArithmeticBasketExp(T,omega,K);
% 
% S = log(ones(N,1)*100);
% t = 0;
% dt = 0.1;
% 
% J = 500;
% order = 1;
% pricer = PDEExpansionPricer(J,dt,order);
% VPDE = pricer.price(S,t,assetModel,derivative);
% VPDE
% 
% numPaths = 10^7;
% pricer = MCPricer(numPaths,dt);
% [VMC,sMC] = pricer.price(S,t,assetModel,derivative);
% VMC,sMC
% 
% res = [VMC,sMC,VPDE];

%% Example: European Arithmetic and Geometric Basket in an Exponential Volatility Model

% N = 5;
% T = 1;
% assetModel = LogExpVolModel(N);
% 
% omega = [ones(N,1)/N;zeros(N,1)/N];
% K = 100;
% derivative = ArithmeticBasketExp(T,omega,K);
% %derivative = GeometricBasketExp(T,omega,K);
% 
% S0 = log(ones(N,1)*100);
% Y0 = -1.5*ones(N,1);
% S = [S0;Y0];
% t = 0;
% dt = 0.01;
% 
% J = 300;
% order = 1;
% pricer = PDEExpansionPricer(J,dt,order,'HV');
% VPDE = pricer.price(S,t,assetModel,derivative);
% VPDE
% 
% numPaths = 10^6;
% pricer = MCPricer(numPaths,dt);
% [VMC,sMC] = pricer.price(S,t,assetModel,derivative);
% VMC,sMC
% 
% res = [VMC,sMC,VPDE];

%% European Geometric Basket in a Variable Coefficient Stock Model

% N = 10; 
% T = 1;
% sigma = @(S,t) (0.1*ones(N,1)+.1*(0:N-1)'/(N-1)) * (1+t/T)/2;
% rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.25);
% mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
% assetModel = VarCoeffModel(N,mu,sigma,rho);
% assetModel = assetModel.setIsConstCorrelationModel(true);
% assetModel = assetModel.setIsVarSpaceModel(false);
% 
% K = 100;
% omega = [2*ones(N/2,1);ones(N/2,1)]/(N*3/2);
% derivative = GeometricBasketExp(T,omega,K);
% 
% S = log(ones(N,1)*100);
% t = 0;
% dt = 0.01;
% 
% assetModel = assetModel.precompute(zeros(N,1),linspace(0,T,T/dt+1));
% 
% J = 300;
% order = 1;
% pricer = PDEExpansionPricer(J,dt,order,'ADI');
% VPDE = pricer.price(S,t,assetModel,derivative);
% 
% numPaths = 10^6;
% pricer = MCPricer(numPaths,dt);
% [VMC,sMC] = pricer.price(S,t,assetModel,derivative);
% 
% res = [VMC,sMC,VPDE];

%% Example: European Arithmetic Baskets in Variable Coefficient Stock Models

N = 10;
T = 1;

assetModel = 0;
model = 'time-dependent volatilities';
if strcmp(model,'const')
    % Log Const - tested against first example above
    sigma = @(S,t) 0.2*ones(N,1);
    rho = @(S,t) AssetModel.buildSimpleConstantCorrelationMatrix(N,0.7);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'const_1')
    % Log Const - tested against first example above
    sigma = @(S,t) 0.2*ones(N,1);
    rho = @(S,t) AssetModel.buildSimpleConstantCorrelationMatrix(N,0.8 - 0.2*4*(0/T-.5)^2);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'const_2')
    % Log Const - tested against first example above
    sigma = @(S,t) 0.2*ones(N,1);
    rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.25 - 0.15*4*(0/T-.5)^2);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'const_3')
    % Log Const - tested against first example above
    sigma = @(S,t) (0.1*ones(N,1)+.1*(0:N-1)'/(N-1)) * (1+0/T);
    rho = @(S,t) AssetModel.buildSimpleConstantCorrelationMatrix(N,0.7);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'const_3_2')
    % Log Const - tested against first example above
    sigma = @(S,t) (0.1*ones(N,1)+.1*(0:N-1)'/(N-1)) * (1+.5*T/T);
    rho = @(S,t) AssetModel.buildSimpleConstantCorrelationMatrix(N,0.7);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'const_3_3')
    % Log Const - tested against first example above
    sigma = @(S,t) (0.1*ones(N,1)+.1*(0:N-1)'/(N-1)) * (1+0);
    rho = @(S,t) AssetModel.buildSimpleConstantCorrelationMatrix(N,0.7);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'const_4')
    % Log Const - tested against first example above
    sigma = @(S,t) (0.1*ones(N,1)+.1*(0:N-1)'/(N-1)) * (1+0/T);
    rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.25 - 0.15*4*(0/T-.5)^2);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'const_4_2')
    % Log Const - tested against first example above
    sigma = @(S,t) (0.1*ones(N,1)+.1*(0:N-1)'/(N-1)) * (1+.5/T);
    rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.25 - 0.15*4*(.5*T/T-.5)^2);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model, 'const_5')
    sigma = @(S,t) 0.2*ones(N,1);
    S_0 = log(ones(N,1)*100);
    rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.1 + 0.1*exp( -sum((S_0-100).^2/(20)^2)/N ));
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'time-dependent volatilities')
    % Time-dependent volatilities
    sigma = @(S,t) (0.1*ones(N,1)+.1*(0:N-1)'/(N-1)) * (1+t/T);
    rho = @(S,t) AssetModel.buildSimpleConstantCorrelationMatrix(N,0.7);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstCorrelationModel(true);
elseif strcmp(model,'time-dependent simple correlation')
    % Time-dependent correlation
    sigma = @(S,t) 0.2*ones(N,1);
    rho = @(S,t) AssetModel.buildSimpleConstantCorrelationMatrix(N,0.8 - 0.2*4*(t/T-.5)^2);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'time-dependent exp correlation')
    % Time-dependent correlation
    sigma = @(S,t) 0.2*ones(N,1);
    rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.25 - 0.15*4*(t/T-.5)^2);
    %rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.2 - 0.1*4*(t/T-.5)^2);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
elseif strcmp(model,'mixed case')
    sigma = @(S,t) (0.1*ones(N,1)+.1*(0:N-1)'/(N-1)) * (1+t/T);
    rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.25 - 0.15*4*(t/T-.5)^2);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(false);
elseif strcmp(model, 'var space')
    sigma = @(S,t) (0.1*ones(N,1) + 0.2*exp(-(S-100).^2/25));
    rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.2);
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsConstCorrelationModel(true);
elseif strcmp(model, 'space corr')
    sigma = @(S,t) 0.2*ones(N,1);
    rho = @(S,t) AssetModel.buildSimpleExponentialCorrelationMatrix(N,0.1 + 0.1*exp( -sum((S-100).^2/(20)^2)/N ));
    mu = @(S,t) zeros(N,1) - .5*sigma(S,t).^2;
    assetModel = VarCoeffModel(N,mu,sigma,rho);
    assetModel = assetModel.setIsVarSpaceModel(true);
    assetModel = assetModel.setIsConstCorrelationModel(false);
    assetModel = assetModel.setIsConstDriftModel(true);
    assetModel = assetModel.setIsConstVolatilityModel(true);
end

omega = ones(N,1)/N;
%omega = [2*ones(N/2,1);ones(N/2,1)]/(N*3/2);
%omega = [ones(7,1);-ones(3,1)]/4;

K = 100;
derivative = ArithmeticBasketExp(T,omega,K);

S = log(ones(N,1)*100);
t = 0;
dt = 0.1;

assetModel = assetModel.precompute(zeros(N,1),linspace(0,T,T/dt+1)); % precompute with Srange being irrelevant here (just needs to indicate N in first size)

numPaths = 10^4;
%pricer = MCPricer(numPaths,dt);
%[VMC,sMC] = pricer.price(S,t,assetModel,derivative);
%VMC,sMC

J = 100;
J3 = 50;
order = 2;
pricer = PDEExpansionPricer(J,dt,order,'ADI');
VPDE_ADI = pricer.price(S,t,assetModel,derivative);
VPDE_ADI

%pricer = PDEExpansionPricer(J,dt,order,'HV');
%VPDE_HV = pricer.price(S,t,assetModel,derivative);
%VPDE_HV

%res = [VMC,sMC,VPDE_ADI,VPDE_HV];

%% Example: European Arithmetic and Geometric Basket in the Heston-Model

% N = 3;
% T = 1;
% assetModel = LogHestonModel(N);
% 
% omega = [ones(N,1)/N;zeros(N,1)/N];
% K = 1;
% %derivative = ArithmeticBasketExp(T,omega,K);
% derivative = GeometricBasketExp(T,omega,K);
% 
% S0 = log(ones(N,1)*1);
% v0 = log([.2723,.2536,.1539]');
% S = [S0;v0];
% t = 0;
% dt = 0.05;
% 
% J = 100;
% order = 1;
% pricer = PDEExpansionPricer(J,dt,order,'ADI');
% VPDE = pricer.price(S,t,assetModel,derivative);
% VPDE
% 
% numPaths = 10^5;
% pricer = MCPricer(numPaths,dt);
% [VMC,sMC] = pricer.price(S,t,assetModel,derivative);
% VMC,sMC
% 
% res = [VMC,sMC,VPDE];

end

