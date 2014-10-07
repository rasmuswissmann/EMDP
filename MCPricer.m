classdef MCPricer < Pricer
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        numPaths = 10^5;
        dt = 0.1;
        M_par = 10^5;
        MAX_PAYOUT = 10000;
    end
    
    methods
        function obj = MCPricer(numPaths,dt)
            obj = obj@Pricer();
            if nargin >= 2
                obj.dt = dt;
            end
            if nargin >= 1
                obj.numPaths = numPaths;
            end
        end
        
        function [Value,s] = price(obj,S,t,assetModel,derivative)
            tic
            % PARAMETERS
            dt = obj.dt;
            N = assetModel.getNumberOfAssets();
            T = derivative.getTerminalTime();
            L_initial = S(1:N)';
            M = T/dt;
            if (M~=ceil(M))
                M = ceil(M);
                dt = T/ceil(M);
                disp(strcat('Warning: Decreased dt from '),num2str(obj.dt),' to ',num2str(dt),' for computation to make T/dt integer.');
            end
            numTimeSteps = M+1;
            numPaths = obj.numPaths;
            M_par = obj.M_par;
            if round(numPaths/M_par)~=(numPaths/M_par)
                while (round(numPaths/M_par)~=(numPaths/M_par))
                    M_par = M_par/10;
                    if (M_par == 100)
                        break;
                    end
                end
                if (M_par == 100)
                    M_par = 1000;
                    numPaths = ceil(numPaths/M_par) * M_par;
                    disp('Warning: Increased numPaths from ', num2str(obj.numPaths) , ' to ', num2str(numPaths) ,' for computation to make it a multiple of 1000.');
                end
            end
            
            stream = RandStream('mt19937ar','seed',sum(100*clock));
            RandStream.setGlobalStream(stream);
            
            %sigm = zeros(N,numTimeSteps);
            %gamma = zeros(N,N,numTimeSteps);
            %d = N;
            %for j=1:numTimeSteps
            %    sigm(:,j) = assetModel.getVolatility(S,(j-1)*dt);
            %    rho = assetModel.getCorrelation(S,(j-1)*dt);
            %    gamma(:,:,j) = MCPricer.gamma_d(rho,d); % TODO in model?
            %end
            
            L_initial = repmat(L_initial,M_par,1);
            
            Value = 0;
            s = 0;
            %% MC Loop
            for m=1:(numPaths/M_par)
                L = MCPricer.MC(L_initial,T,dt,assetModel);
                Values = min(derivative.getPayout(L(:,:,numTimeSteps)',0),obj.MAX_PAYOUT);
                Value = Value + sum(Values);
                s = s + sum(Values.^2);
                toc
            end
            
            % Final Calculations
            Value = Value/numPaths;
            s = sqrt((s-numPaths*Value.^2)/(numPaths*(numPaths-1)));
            
            toc
        end
    end
    
    methods(Static = true) % Internal subfunctions
        
        function [L] = MC(L_initial,T,dt,assetModel)
            
            N = size(L_initial,2);
            M_par = size(L_initial,1);
            numTimeSteps = T/dt;
            L = zeros(M_par,N,numTimeSteps);
            L(:,:,1) = L_initial;
            sdt = sqrt(dt);
            isConstCorrelation = assetModel.isConstCorrelationModel();
            if (isConstCorrelation)
                cholRho = sdt*chol(assetModel.getCorrelation(0,0))';
            end
            isConstVolatility = assetModel.isConstVolatilityModel();
            if (isConstVolatility)
                sigma = assetModel.getVolatility(0,0);
                sigma = sigma*ones(1,M_par);
            end
            isConstDrift = assetModel.isConstDriftModel();
            if (isConstDrift)
                drift = dt*assetModel.getDrift(0,0);
                drift = drift*ones(1,M_par);
                drift = drift';
            end
            
            %if (isConstDrift && isConstColatility && isConstCorrelation)
            %TODO: one step computation
            %else
            % Time iteration
            
            % The following takes properties of the assetModel into
            % account. There is no need to make those 'ifs' here, but it
            % speed up computations tremendously.
            vec2 = ones(M_par,N);
            for j=1:numTimeSteps
                
                % Drift
                if (~isConstDrift)
                    if (assetModel.isVarSpaceModel)
                        drift = dt*assetModel.getDrift(L(:,:,j)',(j-1)*dt);
                    else
                        drift = dt*assetModel.getDrift(L(1,:,j)',(j-1)*dt);
                        drift = repmat(drift,1,M_par);
                    end
                    drift = drift';
                end
            
                % Brownian motion
                dW = randn(M_par,N);
                if (~isConstCorrelation)
                    if (assetModel.isVarSpaceModel)
                        for i=1:M_par
                            cholRho = sdt*chol(assetModel.getCorrelation(L(i,:,j)',(j-1)*dt))';
                            dW(i,:) = dW(i,:)*(cholRho');
                        end
                    else
                        cholRho = sdt*chol(assetModel.getCorrelation(zeros(N,1),(j-1)*dt))';
                        dW = dW*(cholRho');
                    end
                else
                    dW = dW*(cholRho');
                end
                
                % Volatility
                if (~isConstVolatility)
                    if (assetModel.isVarSpaceModel)
                        sigma = assetModel.getVolatility(L(:,:,j)',(j-1)*dt);
                    else
                        sigma = assetModel.getVolatility(L(1,:,j)',(j-1)*dt);
                        sigma = repmat(sigma,1,M_par);
                    end
                end
                
                if (assetModel.isDriftTimesS() && assetModel.isVolTimesS())
                    L(:,:,j+1) = L(:,:,j).*(vec2 + drift + sigma'.*dW);
                elseif (assetModel.isDriftTimesS() && ~assetModel.isVolTimesS())
                    L(:,:,j+1) = L(:,:,j).*(vec2 + drift) + sigma'.*dW;
                elseif (~assetModel.isDriftTimesS() && assetModel.isVolTimesS())
                    L(:,:,j+1) = L(:,:,j).*(vec2 + sigma'.*dW) + drift;
                else
                    L(:,:,j+1) = L(:,:,j) + drift + sigma'.*dW;
                end
                
                % Assure positivty? 
                % L(:,:,j+1) = max(L(:,:,j+1),0);
                
            end
            %end
            
        end
        
        function [gamma] = gamma_d(rho,d)
            [Q,L] = eig(rho);
            L = diag(L);
            [L,IX] = sort(L,'descend');
            gamma = Q(:,IX(1:d))*diag(sqrt(L(1:d)));
            for n=1:size(gamma,1)
                gamma(n,:) = gamma(n,:)/norm(gamma(n,:));
            end
            gamma = gamma';
        end
        
    end
    
end

