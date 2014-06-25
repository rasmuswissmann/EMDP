classdef VarCoeffModel < AssetModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigma;
        rho;
        mu;
    end
    
    methods
        function obj = VarCoeffModel(N,mu,sigma,rho)
            obj = obj@AssetModel(N);
            obj.ranges = ones(N,1)*[-Inf,Inf,log(1),log(1000)];
            obj.isDriftTimesS = false;
            obj.isVolTimesS = false;
            obj.isConstCorrelation = false;
            obj.isConstVolatility = false;
            obj.isConstDrift = false;
            
            % These are functions of the form f(S,t)
            obj.mu = mu;
            obj.sigma = sigma;
            obj.rho = rho;
        end
        
        function sigma = calculateVolatility(obj,S,t)
             if size(S,2)==1
                sigma = obj.sigma(S,t);
             else
                sigma = zeros(size(S));
                for i=1:size(S,2)
                    sigma(:,i) = obj.sigma(S(:,i),t);
                end
             end
        end
        
        function rho = calculateCorrelation(obj,S,t)
            if size(S,2)==1
                rho = obj.rho(S,t);
            else
                rho = zeros(size(S,1),size(S,1),size(S,2));
                for i=1:size(S,2)
                    rho(:,:,i) = obj.rho(S(:,i),t);
                end
            end
        end
        
        function mu = calculateDrift(obj,S,t)
            if size(S,2)==1
                mu = obj.mu(S,t);
            else
                mu = zeros(size(S));
                for i=1:size(S,2)
                    mu(:,i) = obj.mu(S(:,i),t);
                end  
            end
        end
    end
    
end

