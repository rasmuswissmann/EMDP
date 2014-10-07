classdef ConstCoeffModel < AssetModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigma;
        rho;
        mu;
        Sigma;
    end
    
    methods
        function obj = ConstCoeffModel(N,sigma,rho,mu)
            obj = obj@AssetModel(N);
            if nargin > 0
                if (length(sigma) == 1)
                    sigma = sigma*ones(obj.N,1);
                end
                if (length(rho) == 1)
                    rho = AssetModel.buildSimpleConstantCorrelationMatrix(obj.N,rho);
                end
                if (length(mu) == 1)
                    mu = mu*ones(obj.N,1);
                end
                obj.sigma = sigma;
                obj.rho = rho;
                obj.mu = mu;
                obj.isConstCorrelation = true;
                obj.isConstVolatility = true;
                obj.isConstDrift = true;
                obj.isVarSpace = false;
                obj.isDriftTimesS = true;
                obj.isVolTimesS = true;
            end
            obj.Sigma = AssetModel.buildCovarianceMatrix(obj.sigma,obj.rho);
        end
        
        function sigma = calculateVolatility(obj,~,~)
            sigma = obj.sigma;
        end
        
        function rho = calculateCorrelation(obj,~,~)
            rho = obj.rho;
        end
        
        function Sigma = getCovariance(obj,~,~)
            Sigma = obj.Sigma;
        end
        
        function mu = calculateDrift(obj,~,~)
            mu = obj.mu;
        end
    end
    
end

