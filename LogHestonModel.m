classdef LogHestonModel < AssetModel
    % both S and vol in log coordinates
    
    %%% Using this model successfully will require previously modifying the
    %%% current MCPricer to guarantee positivity of volatilites and the
    %%% PDEExpansionPrice to deal with boundary conditions.
    
    properties
        r = [0,0,0]';
        q = [0,0,0]';
        kappa = [1.0121,0.5217,0.5764]';
        eta = [.7627,.4611,.3736]';
        p = [-.9137,-.9322,-.4835]';
        %rho = [1,-0.913700000000000,0.0234000000000000,-0.0219000000000000,0.0625000000000000,-0.0302000000000000;-0.913700000000000,1,-0.0214000000000000,0.0200000000000000,-0.0571000000000000,0.0276000000000000;0.0234000000000000,-0.0214000000000000,1,-0.932200000000000,0.693400000000000,-0.335200000000000;-0.0219000000000000,0.0200000000000000,-0.932200000000000,1,-0.646400000000000,0.312500000000000;0.0625000000000000,-0.0571000000000000,0.693400000000000,-0.646400000000000,1,-0.483500000000000;-0.0302000000000000,0.0276000000000000,-0.335200000000000,0.312500000000000,-0.483500000000000,1;];
        rho = [1.0000    ,0.0234  ,  0.0625  , -0.9137 ,  -0.0219 ,  -0.0302;    0.0234 ,   1.0000,    0.6934 ,  -0.0214  , -0.9322 ,  -0.3352;    0.0625  ,  0.6934  ,  1.0000 ,  -0.0571  , -0.6464 ,  -0.4835;   -0.9137  , -0.0214  , -0.0571 ,   1.0000  ,  0.0200 ,   0.0276;   -0.0219  , -0.9322  , -0.6464 ,   0.0200  ,  1.0000 ,   0.3125;   -0.0302  , -0.3352  , -0.4835 ,   0.0276  ,  0.3125 ,   1.0000];
        theta = [.2874,.2038,.1211]';
    end
    
    methods
        function obj = LogHestonModel(N)
            obj = obj@AssetModel(2*N); % N stocks, N volatilities
            obj.ranges = [ones(N,1)*[-Inf,Inf,log(10),log(1000)]; [ones(N,1)*[-Inf,Inf],log(obj.theta/2),log(obj.theta*2)]];
            obj.isConstCorrelation = true;
            obj.isConstVolatility = false;
            obj.isConstDrift = false;
            obj.isVarSpace = true;
            obj.isDriftTimesS = false;
            obj.isVolTimesS = false;
        end
        
        function sigma = calculateVolatility(obj,S,~)
            N = obj.N;
            eta = obj.eta;
            vec = ones(1,size(S,2));
            sigma = zeros(size(S));
            sigma(1:(N/2),:) = exp(S((1+(N/2)):N,:)/2);
            sigma((1+(N/2)):N,:) = (eta*vec)./sigma(1:(N/2),:);
            sigma(sigma == Inf) = 10^16; %vol of vol can go to inf for vol -> 0
        end
        
        function rho = calculateCorrelation(obj,~,~)
            rho = obj.rho;
        end
        
        function mu = calculateDrift(obj,S,~)
            r = obj.r;
            q = obj.q;
            N = obj.N;
            kappa = obj.kappa;
            theta = obj.theta;
            eta = obj.eta;
            vec = ones(1,size(S,2));
            mu = (r-q)*vec-1/2*S((1+N/2):N,:);
            expMinusLogVol = exp(-S((1+N/2):N,:));
            mu = [mu;diag(kappa)*((theta*vec).*expMinusLogVol-1)-1/2*((eta.^2)*vec).*expMinusLogVol];
        end
        
    end
    
end

