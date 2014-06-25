classdef Derivative
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T; % Terminal time
        omega; % Vector of asset weights (or approximation of it)
    end
    
    
    methods
        function obj = Derivative(T)
            if nargin > 0
                obj.T = T;
            end
        end
        
        % should be vectorized
        getPayout(obj,S,t);
        
        function T = getTerminalTime(obj)
            T = obj.T;
        end
        
        function omega = getAssetWeights(obj)
            omega = obj.omega;
        end
        
    end
    
end