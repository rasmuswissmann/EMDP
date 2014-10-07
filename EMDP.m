function [res] = EMDP()

format long
warning off all

res = [];

%% Example: European Arithmetic Basket in a Log-Normal Stock Model

N = 10;
T = 1;
sigma = 0.2;
rho = AssetModel.buildSimpleConstantCorrelationMatrix(N,0.7);
assetModel = LogConstCoeffLogNormalModel(N,sigma,rho,0);
% assetModel = assetModel.setRanges([-Inf,Inf,log(10),log(1000)]);

omega = ones(N,1)/N;
K = 100;
derivative = ArithmeticBasketExp(T,omega,K);

S = log(ones(N,1)*100);
t = 0;
dt = 0.1;

J = 500;
order = 1;
pricer = PDEExpansionPricer(J,dt,order);
VPDE = pricer.price(S,t,assetModel,derivative);
VPDE

numPaths = 10^7;
pricer = MCPricer(numPaths,dt);
[VMC,sMC] = pricer.price(S,t,assetModel,derivative);
VMC,sMC

res = [VMC,sMC,VPDE];

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
% VPDE
% 
% numPaths = 10^4;
% pricer = MCPricer(numPaths,dt);
% [VMC,sMC] = pricer.price(S,t,assetModel,derivative);
% VMC,sMC
% 
% res = [VMC,sMC,VPDE];

end

