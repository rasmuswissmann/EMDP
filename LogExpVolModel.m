classdef LogExpVolModel < AssetModel
    % both S and vol in log coordinates
    
    properties
        xi;
        ybar;
        theta;
        rho;
    end
    
    methods
        function obj = LogExpVolModel(N)
            obj = obj@AssetModel(2*N); % N stocks, N volatilities
            
	          obj.xi = 0.6*ones(N,1);
	          obj.ybar = -1.5*ones(N,1);
	          obj.theta = 0.5*ones(N,1);
            rho1 = AssetModel.buildSimpleConstantCorrelationMatrix(N,0.85);
            rho2 = eye(N);
            obj.rho = [rho1,zeros(N,N);zeros(N,N),rho2];
            
            obj.isConstCorrelation = true;
            obj.isConstVolatility = false;
            obj.isConstDrift = false;
            obj.isVarSpace = true;
            obj.isDriftTimesS = false;
            obj.isVolTimesS = false;
        end
        
        function sigma = calculateVolatility(obj,S,~)
            N = obj.N;
            vec = ones(1,size(S,2));
            sigma = zeros(size(S));
            sigma(1:(N/2),:) = exp(S((1+(N/2)):N,:));
            sigma((1+(N/2)):N,:) = obj.theta*vec;
        end
        
        function rho = calculateCorrelation(obj,~,~)
            rho = obj.rho;
        end
        
        function mu = calculateDrift(obj,S,~)
            N = obj.N;
            vec = ones(1,size(S,2));
            mu1 = -1/2*exp(2*S((1+(N/2)):N,:));
            mu2 = diag(obj.xi)*(obj.ybar*vec-S((1+(N/2)):N,:));
            mu = [mu1;mu2];
        end
        
    end
    
end

