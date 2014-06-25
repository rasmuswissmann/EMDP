classdef ArithmeticBasket < Derivative
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Strike;
    end
    
    methods
        function obj = ArithmeticBasket(T,omega,Strike)
            obj = obj@Derivative(T);
            if nargin > 0
                obj.omega = omega;
                obj.Strike = Strike;
            end
        end
        
        function v = getPayout(obj,S,~)
            vec = ones(1,size(S,2));
            v = max(obj.omega'*S-obj.Strike*vec,0*vec);
        end
    end
    
end

