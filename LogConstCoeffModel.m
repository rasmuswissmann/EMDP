classdef LogConstCoeffModel < ConstCoeffModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = LogConstCoeffModel(N,sigma,rho,mu)
            obj = obj@ConstCoeffModel(N,sigma,rho,mu-sigma.^2/2);
            obj.isDriftTimesS = false;
            obj.isVolTimesS = false;
        end
        
    end
    
end

