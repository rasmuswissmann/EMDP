classdef ArithmeticBasketExp < ArithmeticBasket
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = ArithmeticBasketExp(T,omega,Strike)
            obj = obj@ArithmeticBasket(T,omega,Strike);
        end
        
        function v = getPayout(obj,S,~)
            v = getPayout@ArithmeticBasket(obj,exp(S));
            %vec = ones(1,size(S,2));
            %v = max(obj.mu'*S-obj.Strike*vec,0*vec);
        end
    end
    
end

