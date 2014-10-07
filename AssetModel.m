classdef AssetModel
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N;
        % Optimizations to speed up calculations
        isConstCorrelation = false;
        isConstVolatility = false;
        isConstDrift = false;
        isVarSpace = true;
        isDriftTimesS = false;
        isVolTimesS = false;
        ranges = [];
        rhopre;
        sigmapre;
        mupre;
        S_range;
        t_range;
        isPreloaded = false;
        MAXVOL = 10000; % Cap max volatility to avoid Inf/NaN in covariance matrix
    end
    
    methods
        function obj = AssetModel(N)
            if nargin > 0
                obj.N = N;
                obj.ranges = ones(N,1)*[0,Inf,1,1000]; % Specifies for each asset: min, max, min_focus_range, max_focus_range
            end
        end
        
        % If not const, then these need to be vectorised in the sense of being able to handle multiple S values at once
        calculateVolatility(obj,S,t);
        calculateCorrelation(obj,S,t);
        calculateDrift(obj,S,t);
        
        function sigma = getVolatility(obj,S,t)
            if obj.isPreloaded
                if obj.isConstVolatility % const 
                    sigma = obj.sigmapre;
                elseif (obj.isVarSpaceModel()==false) % time-dependent only
                    [~, i] = min(abs(obj.t_range - t));
                    sigma = obj.sigmapre(:,:,i);
                else % time- and space-dependent 
                    % TODO
                    sigma = obj.calculateVolatility(S,t);
                end
            else
                sigma = obj.calculateVolatility(S,t);
            end
            sigma = min(sigma,obj.MAXVOL);
        end
        
        function rho = getCorrelation(obj,S,t)
            if obj.isPreloaded
                if obj.isConstCorrelation % const 
                    rho = obj.rhopre;
                elseif (obj.isVarSpaceModel()==false) % time-dependent only
                    [~, i] = min(abs(obj.t_range - t));
                    rho = obj.rhopre(:,:,:,i);
                    %rho = squeeze(obj.rhopre(:,:,:,i));
                else % time- and space-dependent
                    % TODO
                    rho = obj.calculateCorrelation(S,t);
                end
            else
                rho = obj.calculateCorrelation(S,t);
            end
        end
        
        function mu = getDrift(obj,S,t)
            if obj.isPreloaded
                if obj.isConstDrift % const 
                    mu = obj.sigmapre;
                elseif (obj.isVarSpaceModel()==false) % time-dependent only
                    [~, i] = min(abs(obj.t_range - t));
                    mu = obj.mupre(:,:,i);
                else % time- and space-dependent
                    % TODO
                    mu = obj.calculateDrift(S,t);
                end
            else
                mu = obj.calculateDrift(S,t);
            end
        end
        
        
        function Sigma = getCovariance(obj,S,t)
            Sigma = AssetModel.buildCovarianceMatrix(obj.getVolatility(S,t),obj.getCorrelation(S,t));
        end
        
        function N = getNumberOfAssets(obj)
            N = obj.N;
        end
        
        function b = isConstCorrelationModel(obj)
            b = obj.isConstCorrelation;
        end
        
        function b = isConstVolatilityModel(obj)
            b = obj.isConstVolatility;
        end
        
        function b = isConstDriftModel(obj)
            b = obj.isConstDrift;
        end
        
        function b = isVarSpaceModel(obj)
            b = obj.isVarSpace;
        end
        
        function b = isConstCovarrianceModel(obj)
            b = obj.isConstCorrelationModel() * obj.isConstVolatilityModel();
        end
        
        function b = isDriftTimesSModel(obj)
            b = obj.isDriftTimesS;
        end
        
        function b = isVolTimesSModel(obj)
            b = obj.isVolTimesS;
        end
        
        function obj = setIsVarSpaceModel(obj,b)
            obj.isVarSpace = b;
        end
        
        function obj = setIsConstCorrelationModel(obj,b)
            obj.isConstCorrelation = b;
        end
        
        function obj = setIsConstVolatilityModel(obj,b)
            obj.isConstVolatility = b;
        end
        
        function obj = setIsConstDriftModel(obj,b)
            obj.isConstDrift = b;
        end
        
        function obj = set_ranges(obj,range)
            if (size(range,1) == 1)
                obj.ranges = ones(obj.N,1)*range;
            else
                obj.ranges = range;
            end
        end
        
        function ranges = get_ranges(obj)
            ranges = obj.ranges;
        end
        
        function obj = precompute(obj,S_range,t_range)
            if obj.isConstDriftModel()
                obj.mupre = obj.getDrift(S_range,0);
            else
                obj.mupre = zeros(size(S_range,1),size(S_range,2),length(t_range));
            end
            if obj.isConstVolatilityModel()
                obj.sigmapre = obj.getVolatility(S_range,0);
            else
                obj.sigmapre = zeros(size(S_range,1),size(S_range,2),length(t_range));
            end
            if obj.isConstCorrelationModel()
                obj.rhopre = obj.getCorrelation(S_range,0);
            else
                obj.rhopre = zeros(size(S_range,1),size(S_range,1),size(S_range,2),length(t_range));
            end
            for i=1:length(t_range)
                t = t_range(i);
                if ~obj.isConstDriftModel() 
                    obj.mupre(:,:,i) = obj.getDrift(S_range,t); 
                end
                if ~obj.isConstVolatilityModel() 
                    obj.sigmapre(:,:,i) = obj.getVolatility(S_range,t);
                end
                if ~obj.isConstCorrelationModel() 
                    obj.rhopre(:,:,:,i) = obj.getCorrelation(S_range,t);
                end
            end
            obj.S_range = S_range;
            obj.t_range = t_range;
            obj.isPreloaded = true;
        end
        
    end
    
    methods(Static = true)
        function Sigma = buildCovarianceMatrix(sigma,rho)
            vol = diag(sigma);
            Sigma = vol*rho*vol;
        end
        
        function rho = buildSimpleConstantCorrelationMatrix(N,r)
            rho = r*ones(N,N) + (1-r)*diag(ones(N,1));
        end
        
        function rho = buildSimpleExponentialCorrelationMatrix(N,r)
            rho = exp(-r*abs(ones(N,N)*diag(1:N)-diag(1:N)*ones(N,N)));
        end
        
        function gamma = gamma_d(rho,d)
            [Q,L] = eig(rho);
            L = diag(L);
            [L,IX] = sort(L,'descend');
            gamma = Q(:,IX(1:d))*diag(sqrt(L(1:d)));
            for n=1:size(gamma,1)
                gamma(n,:) = gamma(n,:)/norm(gamma(n,:));
            end
            gamma = gamma';
            %rho = gamma'*gamma;
        end
    end
    
end

