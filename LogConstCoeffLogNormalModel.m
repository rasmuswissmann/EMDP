classdef LogConstCoeffLogNormalModel < ConstCoeffLogNormalModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = LogConstCoeffLogNormalModel(N,sigma,rho,mu)
            obj = obj@ConstCoeffLogNormalModel(N,sigma,rho,mu-sigma.^2/2);
            obj.isDriftTimesS = false;
            obj.isVolTimesS = false;
        end
        
    end
    
end

