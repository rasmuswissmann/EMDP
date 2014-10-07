classdef PDEExpansionPricer < Pricer
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        J = 400;
        J3 = 400;
        dt = 0.1;
        scheme = 'ADI'; % ADI or HV
        ADIMethod3D = 'Brian'; % Douglas or Brian ADI for 3D problems
        MAX_PAYOUT = 10000; % Maximal payout used to truncate solution towards infinity, if necessary
        order = 1;
        methodForQ = 'PCA'; % PCA or LT (LT not implemented yet)
    end
    
    methods
        function obj = PDEExpansionPricer(J,dt,order,scheme,J3,ADIMethod3D,MAX_PAYOUT)
            obj = obj@Pricer();
            if nargin >= 1
                obj.J = J;
                obj.J3 = J;
            end
            if nargin >= 2
                obj.dt = dt;
            end
            if nargin >= 3
                obj.order = order;
            end
            if nargin >= 4
                obj.scheme = scheme;
            end
            if nargin >= 5
                obj.J3 = J3;
            end
            if nargin >= 6
                obj.ADIMethod3D = ADIMethod3D;
            end
            if nargin >= 7
                obj.MAX_PAYOUT = MAX_PAYOUT;
            end
        end
        
        function value = price(obj,S,~,assetModel,derivative)
            tic
            % PARAMETERS
            scheme = obj.scheme;
            methodForQ = obj.methodForQ;
            dt = obj.dt;
            J1 = obj.J;
            J2 = obj.J;
            J3 = obj.J3;
            N = assetModel.getNumberOfAssets();
            T = derivative.getTerminalTime();
            L = S(1:N);
            M = T/dt;
            if (M~=ceil(M))
                M = ceil(M);
                dt = T/ceil(M);
                disp(strcat('Decreased dt from '),num2str(obj.dt),' to ',num2str(dt),' for computation to make T/dt integer.');
            end
            omega = derivative.getAssetWeights();
            %TODO: implement t
            
            % PROBLEM SET UP AND COORDINATE TRANSFORMATIONS
            Q = PDEExpansionPricer.principalComponentTrafo(N,L,T,assetModel,omega,methodForQ);
            [beta,z] = PDEExpansionPricer.precomputeDrift(Q,L,T,M,assetModel);
            % TODO: this needs to be via ranges or similar, so as to not be
            % tied to exp/log models
            Lplus = log(exp(S)*5); Lminus = log(exp(S)/5);
            v0 = z + beta(:,M+1);
            [gamma,y1,y2] = PDEExpansionPricer.arctanTrafoOld(N,M,beta,z,Lplus,Lminus,Q,J1,J2); %TODO: parameters, and for J3!
            
            % 1-DIM CASE for first EV
            toc
            disp('1-dim computation');
            A_1 = PDEExpansionPricer.matrix1Dim(J1,gamma);
            L_grid = PDEExpansionPricer.L_grid1Dim(N,J1,Q,y1,beta,v0,1);
            lambda = zeros(N,N);
            if (assetModel.isVarSpaceModel())
                lambda = zeros(N,N,size(L_grid,1));
            end
            
            u1 = PDEExpansionPricer.initialCondition(L_grid,derivative,obj.MAX_PAYOUT);
            for ct=1:M
                if ((ct==1) || ~assetModel.isConstCovarrianceModel())
                    LRange = L';
                    if (assetModel.isVarSpaceModel())
                        LRange = L_grid;
                    end
                    lambda = PDEExpansionPricer.updateLambda(lambda,assetModel,LRange,T,ct/M*T,Q,1);
                end
                u1 = PDEExpansionPricer.crankNicolson(u1,J1+1,lambda,dt,A_1);
            end
            u1d = PDEExpansionPricer.u_interpL1D(u1,v0,y1);
            u_com = u1d;
            clear u_1d; clear A_1;
            disp(u_com);
            
            if (obj.order >= 1)
                % 2-DIM CASE for 1 and k
                [A1,Ak,A1k] = PDEExpansionPricer.matrices2Dim(J2);
                L_grid2D_1 = PDEExpansionPricer.L_grid2Dim1(N,J2,Q,y2);
                lambda = zeros(N,N);
                if (assetModel.isVarSpaceModel())
                    lambda = zeros(N,N,size(L_grid,1));
                end
                
                u2d = zeros(N,1);
                for k=2:N
                    toc
                    disp(strcat('2-dim computation for k=',num2str(k)));
                    
                    L_grid2D_k = PDEExpansionPricer.L_grid2Dimk(N,J2,Q,y2,k);
                    L_grid = PDEExpansionPricer.L_grid2Dim(N,beta,Q,v0,k,L_grid2D_1+L_grid2D_k,1);
                    u2 = PDEExpansionPricer.initialCondition(L_grid,derivative,obj.MAX_PAYOUT);
                    
                    if strcmp(scheme,'HV')
                        LRange = L';
                        if (assetModel.isVarSpaceModel())
                            LRange = L_grid;
                        end
                        lambda = PDEExpansionPricer.updateLambda(lambda,assetModel,LRange,T,0,Q,[1,k]);
                        [F1,Fk,F0] = PDEExpansionPricer.HVMatrices(lambda,gamma,A1,Ak,A1k,1,k);
                    end
                    for ct=1:M
                        if ((ct==1) || ~assetModel.isConstCovarrianceModel())
                            LRange = L';
                            if (assetModel.isVarSpaceModel())
                                LRange = L_grid;
                            end
                            lambda = PDEExpansionPricer.updateLambda(lambda,assetModel,LRange,T,ct/M*T,Q,[1,k]);
                        end
                        if strcmp(scheme,'ADI')
                            u2 = PDEExpansionPricer.ADI2D(u2,(J2+1)^2,lambda,gamma,dt,A1,Ak,1,k);
                        elseif strcmp(scheme,'HV')
                            [u2,F1,Fk,F0] = PDEExpansionPricer.HV(u2,(J2+1)^2,lambda,gamma,dt,A1,Ak,A1k,F1,Fk,F0,1,k);
                        end
                    end
                    
                    u2d(k) = PDEExpansionPricer.u_interpL2D(u2,k,v0,J2,y2)
                    u_com = u_com + u2d(k) - u1d
                end
                disp(u_com);
            end
            
            J3
            
            if (obj.order >= 2)
                % 3-DIM CASE for 1, k, m
                % This does so far not support spatially varying lambda
                [A1,Ak,Am] = PDEExpansionPricer.matrices3Dim(J3);
                L_grid3D_1 = PDEExpansionPricer.L_grid3Dim1(N,J3,Q,y2);
                lambda = zeros(N,N);
                u3d = zeros(N,N);
                for k=2:N
                    
                    L_grid3D_1k = L_grid3D_1 + PDEExpansionPricer.L_grid3Dimk(N,J3,Q,y2,k);
                    
                    for m=k+1:N
                        toc
                        disp(strcat('3-dim computation for k=',num2str(k),', m=',num2str(m)));
                        
                        L_grid3D_m = PDEExpansionPricer.L_grid3Dimm(N,J3,Q,y2,m);
                        L_grid = PDEExpansionPricer.L_grid3Dim(N,beta,Q,v0,k,m,L_grid3D_1k+L_grid3D_m,1);
                        u3 = PDEExpansionPricer.initialCondition(L_grid,derivative,obj.MAX_PAYOUT);
                        
                        LRange = L';
                        if (assetModel.isConstCovarrianceModel())
                            lambda = PDEExpansionPricer.updateLambda(lambda,assetModel,LRange,T,0/M*T,Q,[1,k,m]);
                            % all M time steps
                            u3 = PDEExpansionPricer.ADI3D(u3,(J3+1)^3,lambda,gamma,dt,M,A1,Ak,Am,1,k,m,obj.ADIMethod3D);
                        else
                            for ct=1:M
                                lambda = PDEExpansionPricer.updateLambda(lambda,assetModel,LRange,T,ct/M*T,Q,[1,k,m]);
                                % just one time step
                                u3 = PDEExpansionPricer.ADI3D(u3,(J3+1)^3,lambda,gamma,dt,1,A1,Ak,Am,1,k,m,obj.ADIMethod3D);
                            end
                        end
                        
                        u3d(k,m) = PDEExpansionPricer.u_interpL3D(u3,k,m,v0,J3,y2)
                        u_com = u_com + u3d(k,m) - u2d(k) - u2d(m) + u1d
                    end
                end
                disp(u_com);
            end
            
            toc
            value = u_com;
            
        end
    end
    
    
    methods(Static = true) % Internal subfunctions
        
        function [u,u2] = initialCondition(L_grid,derivative,MAX_PAYOUT)
            u = zeros(size(L_grid,1),1);
            for i=1:length(u)
                u(i) = derivative.getPayout(L_grid(i,:)',0);
            end
            u = min(u,MAX_PAYOUT);
            % works but unfortunately is even slower:
            %u2 = cell2mat(cellfun(@(x) derivative.getPayout(x',0), mat2cell(L_grid(:,:), ones(1,size(L_grid,1))), 'UniformOutput', false));
            %u2 = min(u2,MAX_PAYOUT);
        end
        
        function [Q] = principalComponentTrafo(N,L,T,assetModel,omega,methodForQ)
            Sigma = assetModel.getCovariance(L,0);
            if strcmp(methodForQ,'PCA')
                [Q,l] = eig(Sigma);
                l = diag(l);
                [l,IX] = sort(l,'descend');
                Q = Q(:,IX)'; % bad idea to use Q'? confuses notation
            elseif strcmp(methodForQ,'LT')
                Cpc = principalComponentTrafo(N,T,L,omega,'PCA')';
                dsl = diag(sqrt(diag(Cpc'*Sigma*Cpc)));
                Cpc = Cpc*dsl;
                %Cpc = chol(Sigma);
                A = orthogonalTrafoLT(N,T,L,Cpc,sigma,omega);
                Q = A'*(Cpc^-1);
            else
                Q = eye(N);
            end
        end
        
        function [beta,z] = precomputeDrift(Q,L,T,M,assetModel)
            % Pre-Computation of Drift Term
            N = assetModel.getNumberOfAssets();
            dt = T/M;
            beta = zeros(N,M+1);
            Int = zeros(N,1);
            for tau_c=1:M
                % TODO more exact integration, e.g. middle?
                mu = assetModel.getDrift(L,T-tau_c*dt);
                beta(:,tau_c+1) = beta(:,tau_c) + dt*Q*mu;
            end
            z = Q*L;
        end
        
        function [A] = orthogonalTrafoLT(N,T,L,Q,sigma,omega)
            L = exp(L);
            r = 0.0;
            mu1 = omega.*L;
            mu2 = (r*ones(N,1)-sigma.^2/2)*T;
            s = N;
            Cch = Q;
            A = zeros(N,N);
            Cbar = diag(mu1.*exp(mu2))*Cch; %to avoid imaginary numbers
            A(:,1) = sum(Cbar,1)';
            A(:,1) = A(:,1)/norm(A(:,1));
            for i=2:N
                mu3 = mu2 + Cch*A(:,i-1);
                Cbar = diag(mu1.*exp(mu3))*Cch;
                A(:,i) = sum(Cbar,1)';
                A(:,i) = A(:,i)/norm(A(:,i));
                while (abs(sum(A(:,i)'*A(:,1:i-1)))>10^-15)
                    A(:,i) = A(:,i);
                    for j=1:i-1
                        A(:,i) = A(:,i) - (A(:,i)'*A(:,j))*A(:,j);
                    end
                    A(:,i) = A(:,i)/norm(A(:,i));
                end
            end
        end
        
        % 1-D grid
        
        function [L_grid] = L_grid1Dim(N,J,Q,y1,beta,v0,tau_c)
            Qinv = PDEExpansionPricer.inverseQ(Q);
            griddim = J+1;
            L_grid_1 = zeros(griddim,N);
            for ii = 1:N
                L_grid_1(:,ii) = Qinv(ii,1)*y1;
            end
            val = -Qinv*beta(:,tau_c);
            val = val+(Qinv(:,2:N)*v0(2:N));
            L_grid = L_grid_1 + repmat(val',griddim,1);
        end
        
        % 2-D grid
        
        function [L_grid_1] = L_grid2Dim1(N,J,Q,y)
            Qinv = PDEExpansionPricer.inverseQ(Q);
            L_grid_1 = zeros((J+1)^2,N);
            for j=1:J+1
                L_grid_1(j:J+1:(J+1)^2,:) = repmat((Qinv(:,1)*y(1,j))',J+1,1);
            end
        end
        
        function [L_grid_k] = L_grid2Dimk(N,J,Q,y,k)
            Qinv = PDEExpansionPricer.inverseQ(Q);
            L_grid_k = zeros((J+1)^2,N);
            for j=1:J+1
                L_grid_k((j-1)*(J+1)+1:j*(J+1),:) = repmat((Qinv(:,k)*y(k,j))',J+1,1);
            end
        end
        
        function [L_grid] = L_grid2Dim(N,beta,Q,v0,k,L_grid,tau_c)
            Qinv = PDEExpansionPricer.inverseQ(Q);
            val = -Qinv*beta(:,tau_c);
            if (k>2) val = val +(Qinv(:,2:k-1)*v0(2:k-1)); end
            if (k<N) val = val + (Qinv(:,k+1:N)*v0(k+1:N)); end
            L_grid = L_grid + repmat(val',size(L_grid,1),1);
        end
        
        % 3-D grid
        
        function [L_grid_1] = L_grid3Dim1(N,J,Q,y)
            Qinv =PDEExpansionPricer.inverseQ(Q);
            L_grid_1 = zeros((J+1)^3,N);
            for j=1:J+1
                ind = j:(J+1):(J+1)^3;
                L_grid_1(ind,:) = repmat((Qinv(:,1)*y(1,j))',(J+1)^2,1);
            end
        end
        
        function [L_grid_k] = L_grid3Dimk(N,J,Q,y,k)
            Qinv =PDEExpansionPricer.inverseQ(Q);
            L_grid_k = zeros((J+1)^3,N);
            for j=1:J+1
                ind0 = (j-1)*(J+1)+1:j*(J+1);
                ind = zeros((J+1)^2,1);
                for m=0:J
                    ind((1:(J+1))+m*(J+1)) = ind0 + m*(J+1)^2;
                end
                L_grid_k(ind,:) = repmat((Qinv(:,k)*y(k,j))',(J+1)^2,1);
            end
        end
        
        function [L_grid_m] = L_grid3Dimm(N,J,Q,y,m)
            Qinv =PDEExpansionPricer.inverseQ(Q);
            L_grid_m = zeros((J+1)^3,N);
            for j=1:J+1
                ind = (j-1)*(J+1)^2+1:j*(J+1)^2;
                L_grid_m(ind,:) = repmat((Qinv(:,m)*y(m,j))',(J+1)^2,1);
            end
        end
        
        function [L_grid] = L_grid3Dim(N,beta,Q,v0,k,m,L_grid,tau_c)
            Qinv =PDEExpansionPricer.inverseQ(Q);
            val = -Qinv*beta(:,tau_c);
            ind = [2:k-1,k+1:m-1,m+1:N];
            val = val+(Qinv(:,ind)*v0(ind));
            L_grid = L_grid+repmat(val',size(L_grid,1),1);
        end
        
        function [lambda] = updateLambda(lambda,assetModel,Ls,T,tau,Q,ind)
            %TODO: without loop?
            for i=1:size(Ls,1)
                L = Ls(i,:)';
                Sigma = assetModel.getCovariance(L,T-tau);
                if (sum(imag(Sigma))~=0)
                    Sigma
                end
                lambda(ind,ind,i) = Q(ind,:)*Sigma*Q(ind,:)';
            end
            lambda(abs(lambda) < 10^-12) = 0;
        end
        
        function [u] = ADI2D(u,dim,lambda,gamma,dt,Aj,Ak,j,k)
            % ADI scheme % could modify by using both prev and cur lambda to handle
            % time-var better, if necessary
            Fjp = speye(dim)+dt/2*lambda(j,j)*(gamma(j)^2)*Aj;
            Fjm = speye(dim)-dt/2*lambda(j,j)*(gamma(j)^2)*Aj;
            [FjmL,FjmU] = lu(Fjm,0);
            Fkp = speye(dim)+dt/2*lambda(k,k)*(gamma(k)^2)*Ak;
            Fkm = speye(dim)-dt/2*lambda(k,k)*(gamma(k)^2)*Ak;
            [FkmL,FkmU] = lu(Fkm,0);
            u = (Fkp * u);
            u = (FjmL \ u);
            u = (FjmU \ u);
            u = (Fjp * u);
            u = (FkmL \ u);
            u = (FkmU \ u);
        end
        
        function [u] = ADI3D(u,dim,lambda,gamma,dt,M,Aj,Ak,Am,j,k,m,ADIMethod3D)
            % ADI scheme
            
            timefac = dt/2;
            
            Fj = timefac*lambda(j,j)*(gamma(j)^2)*Aj;
            Fk = timefac*lambda(k,k)*(gamma(k)^2)*Ak;
            Fm = timefac*lambda(m,m)*(gamma(m)^2)*Am;
            [FjL,FjU] = lu(speye(dim) - Fj,0);
            [FkL,FkU] = lu(speye(dim) - Fk,0);
            [FmL,FmU] = lu(speye(dim) - Fm,0);
            
            for ct=1:M
                if strcmp(ADIMethod3D,'Brian')
                    U = (speye(dim)+Fk+Fm)*u;
                    U = (FjL \ U);
                    U = (FjU \ U);
                    V = (speye(dim)+Fm)*u + Fj*U;
                    V = (FkL \ V);
                    V = (FkU \ V);
                    u = (speye(dim)+Fk)*V + Fj*U;
                    u = (FmL \ u);
                    u = (FmU \ u);
                elseif strcmp(ADIMethod3D,'Douglas')
                    U = (speye(dim)+Fj+2*Fk+2*Fm)*u;
                    U = (FjL \ U);
                    U = (FjU \ U);
                    V = (speye(dim)+Fj+Fk+2*Fm)*u + Fj*U;
                    V = (FkL \ V);
                    V = (FkU \ V);
                    u = (speye(dim)+Fj+Fk+Fm)*u + Fj*U + Fk*V;
                    u = (FmL \ u);
                    u = (FmU \ u);
                end
            end
            
        end
        
        function [u2,Fj,Fk,F0] = HV(u2,dim,lambda,gamma,dt,Aj,Ak,Ajk,Fjp,Fkp,F0p,j,k)
            % HV scheme
            theta = 0.5+sqrt(3)/6;
            mu = 0.5;
            [Fj,Fk,F0] = PDEExpansionPricer.HVMatrices(lambda,gamma,Aj,Ak,Ajk,j,k);
            %Fjp = Fj; Fkp = Fk; F0p = 0*F0;
            y0 = u2 + dt*(Fjp+F0p+Fkp)*u2;
            [Lmat,Umat] = lu((speye(dim)-theta*dt*Fj),0);
            y1 =  Umat \ (Lmat \ (y0-theta*dt*Fjp*u2));
            [Lmat,Umat] = lu((speye(dim)-theta*dt*Fk),0);
            yk = Umat \ (Lmat \ (y1-theta*dt*Fkp*u2));
            y0 = y0 + mu*dt*((Fj+Fk+F0)*yk - (Fjp+Fkp+F0p)*u2);
            [Lmat,Umat] = lu((speye(dim)-theta*dt*Fj),0);
            y1 = Umat \ (Lmat \ (y0-theta*dt*Fj*yk));
            [Lmat,Umat] = lu((speye(dim)-theta*dt*Fk),0);
            u2 = Umat \ (Lmat \ (y1-theta*dt*Fk*yk));
        end
        
        function [Fj,Fk,F0] = HVMatrices(lambda,gamma,Aj,Ak,Ajk,j,k)
            if size(lambda,3) == 1
                Alambdaj = lambda(j,j,1);
                Alambdak = lambda(k,k,1);
                Alambdajk = lambda(j,k,1);
            else
                Alambdaj = spdiags(squeeze(lambda(j,j,:)),0,size(Aj,1),size(Aj,2));
                Alambdak = spdiags(squeeze(lambda(k,k,:)),0,size(Ak,1),size(Ak,2));
                Alambdajk = spdiags(squeeze(lambda(j,k,:)),0,size(Ajk,1),size(Ajk,2));
            end
            Fj = (gamma(j)^2)*Alambdaj*Aj;
            Fk = (gamma(k)^2)*Alambdak*Ak;
            F0 = gamma(j)*gamma(k)*Alambdajk*Ajk;
        end
        
        function u_1d = crankNicolson(u_1d,dim,lambda,dt,A_1)
            % Crank-Nicolson scheme
            if size(lambda,3) == 1
                A_1lambda = lambda(1,1,1);
            else
                A_1lambda = spdiags(squeeze(lambda(1,1,:)),0,size(A_1,1),size(A_1,2));
            end
            A_1lambda = A_1lambda*A_1;
            F1m = speye(dim)-dt/2*A_1lambda;
            F1p = speye(dim)+dt/2*A_1lambda;
            [F1mL,F1mU] = lu(F1m,0);
            u_1d = (F1p * u_1d);
            u_1d = (F1mL \ u_1d);
            u_1d = (F1mU \ u_1d);
        end
        
        function [u_out] = u_interpL1D(u_in,z,vspace)
            u_out = spline(vspace(1,:),u_in(:),z(1));
        end
        
        function [u_out] = u_interpL2D(u_in,k,z,J,vspace)
            utemp = reshape(u_in(:),J+1,J+1);
            [X,Y] = meshgrid(vspace(k,:),vspace(1,:));
            u_out = interp2(X,Y,utemp,z(k),z(1),'spline');
            % X and Y seem inverted because of how interp2 works
        end
        
        function [u_out] = u_interpL3D(u_in,k,m,z,J,vspace)
            utemp = reshape(u_in(:),J+1,J+1,J+1);
            [X,Y,Z] = meshgrid(vspace(k,:),vspace(1,:),vspace(m,:));
            u_out = interp3(X,Y,Z,utemp,z(k),z(1),z(m),'linear');
            % X and Y seem inverted because of how interp3 works
        end
        
        function [Qinv] = inverseQ(Q)
            Qinv =  Q';
            
            % the next operations work for a Q that has
            % orthogonal vectors of arbitrary length
            % use this to avoid calling Q^-1
            %vectornorms = sum(Q.^2,2);
            %Qinv = diag(vectornorms.^-1)*Q;
            %Qinv = Qinv';
            
            % if nothing else applies
            %Qinv = Q^-1;
        end
        
        function [i] = index3D(J,j,k,l)
            i = j + (k-1)*(J+1) + (l-1)*(J+1)^2;
        end
        
        function [j,k,l] = coord3D(J,i)
            j = mod(i,J+1);
            k = (mod(i,(J+1)^2)-j)/(J+1) + 1;
            l = (i-j-(k-1)*(J+1))/(J+1)^2 + 1;
        end
        
        function [trafoParams] = getFiniteTransformationParameters(ranges)
            % adjust ranges
            trafoParams = arctanTrafo(N,M,beta,z,lplus,lminus,Q,J1,J2);
        end
        
        function [gamma,y1,y2] = arctanTrafo(N,M,beta,z,lplus,lminus,Q,J1,J2)
            tau_c = M+1;
            QlnLp = sum(max(Q,0),2).*lplus - sum(max(-Q,0),2).*lminus;
            QlnLm = sum(max(Q,0),2).*lminus - sum(max(-Q,0),2).*lplus;
            wplus = ones(N,1)*0.9;
            wminus = ones(N,1)*0.1;
            gamma = (tan(pi*(wplus-1/2))-tan(pi*(wminus-1/2)))./(QlnLp-QlnLm);
            c = tan(pi*(wplus-1/2)) - gamma.*QlnLp;
            y1 = linspace(0,1,J1+1);
            y1 = 1/gamma(1)*(tan(pi*(y1-0.5))-c(1)*ones(1,J1+1));
            y2 = repmat(linspace(0,1,J2+1),N,1);
            y2 = diag(1./gamma)*(tan(pi*(y2-0.5))-c*ones(1,J2+1));
        end
        
        function [gamma,y1,y2,v0] = arctanTrafoOld(N,M,beta,z,lplus,lminus,Q,J1,J2)
            tau_c = M+1;
            v0 = z + beta(:,tau_c);
            QlnLp = sum(max(Q,0),2).*lplus - sum(max(-Q,0),2).*lminus;
            QlnLp = max(QlnLp,QlnLp + beta(:,tau_c));
            QlnLm = sum(max(Q,0),2).*lminus - sum(max(-Q,0),2).*lplus;
            QlnLm = min(QlnLm,QlnLm + beta(:,tau_c));
            wplus = ones(N,1)*0.9;
            wminus = ones(N,1)*0.1;
            gamma = (tan(pi*(wplus-1/2))-tan(pi*(wminus-1/2)))./(QlnLp-QlnLm);
            c = tan(pi*(wplus-1/2)) - gamma.*QlnLp;
            space_min = 0;
            space_max = 1 - space_min;
            y1 = linspace(space_min,space_max,J1+1);
            y1 = 1/gamma(1)*(tan(pi*(y1-0.5))-c(1)*ones(1,J1+1));
            y2 = repmat(linspace(space_min,space_max,J2+1),N,1);
            y2 = diag(1./gamma)*(tan(pi*(y2-0.5))-c*ones(1,J2+1));
        end
        
        function [A_1] = matrix1Dim(J1,gamma)
            h1 = 1/J1;
            A_1_1 = spdiags(repmat([-1/(2*h1) 0 1/(2*h1)],(J1+1),1),[-1 0 1],(J1+1),(J1+1));
            A_1_2 = spdiags(repmat([1/(h1^2) -2/(h1^2) 1/(h1^2)],(J1+1),1),[-1 0 1],(J1+1),(J1+1));
            vcos = zeros(J1+1,1);
            vsin = zeros(J1+1,1);
            for i=1:J1+1
                vcos(i) = (cos(pi*((i-1)*h1-1/2)));
                vsin(i) = (sin(pi*((i-1)*h1-1/2)));
            end
            Dsin = spdiags(vsin,0,(J1+1),(J1+1));
            Dcos = spdiags(vcos,0,(J1+1),(J1+1));
            A_1 = 1/2*(gamma(1)^2)/pi*(Dcos.^3)*(-2*Dsin*A_1_1+1/pi*Dcos*A_1_2);
        end
        
        function [A1,Ak,A1k] = matrices2Dim(J)
            h = 1/J;
            % Matrix for 2-dim case, 1st direction
            col1 = repmat([-1/(2*h) 0 +1/(2*h)],(J+1),1);
            col1(J+1,1) = 0;
            col1(1,3) = 0;
            % The zeros look like they are in the wrong spots, but this is correct due
            % to how spdiags works
            A1_1 = spdiags(repmat(col1,(J+1),1),[-1 0 1],(J+1)^2,(J+1)^2);
            col2 = repmat([1/(h^2) -2/(h^2) 1/(h^2)],(J+1),1);
            col2(J+1,1) = 0;
            col2(1,3) = 0;
            % The zeros look like they are in the wrong spots, but this is correct due
            % to how spdiags works
            A1_2 = spdiags(repmat(col2,J+1,1),[-1 0 1],(J+1)^2,(J+1)^2);
            vcos = zeros(J+1,1);
            vsin = zeros(J+1,1);
            for i=1:J+1
                vcos(i) = (cos(pi*((i-1)*h-1/2)));
                vsin(i) = (sin(pi*((i-1)*h-1/2)));
            end
            vcos = repmat(vcos,J+1,1);
            vsin = repmat(vsin,J+1,1);
            Dcos = spdiags(vcos,0,(J+1)^2,(J+1)^2);
            Dsin = spdiags(vsin,0,(J+1)^2,(J+1)^2);
            A1 = 1/(2*pi)*(Dcos.^3)*(-2*Dsin*A1_1+1/pi*Dcos*A1_2);
            clear A1_1; clear A1_2;
            Dcos_1 = Dcos;
            
            % Matrix for 2-dim case, direction k
            % TODO: check that border is recpested correctly like for col1, col2 in A1
            col1 = repmat([-1/(2*h) 0 +1/(2*h)],J+1,1);
            Ak_1 = spdiags(repmat(col1,J+1,1),[-1-J 0 1+J],(J+1)^2,(J+1)^2);
            col2 = repmat([1/(h^2) -2/(h^2) 1/(h^2)],(J+1),1);
            Ak_2 = spdiags(repmat(col2,J+1,1),[-1-J 0 1+J],(J+1)^2,(J+1)^2);
            vcos = zeros((J+1)^2,1);
            vsin = zeros((J+1)^2,1);
            for j=1:J+1
                vcos(1+(j-1)*(J+1):j*(J+1)) = (cos(pi*((j-1)*h-1/2)));
                vsin(1+(j-1)*(J+1):j*(J+1)) = (sin(pi*((j-1)*h-1/2)));
            end
            Dcos = spdiags(vcos,0,(J+1)^2,(J+1)^2);
            Dsin = spdiags(vsin,0,(J+1)^2,(J+1)^2);
            Ak = 1/(2*pi)*(Dcos.^3)*(-2*Dsin*Ak_1+1/pi*Dcos*Ak_2);
            
            %% Matrix for 2-dim case, mixed derivative
            % TODO: check that border is recpested correctly like for col1, col2 in A1
            A1k = spdiags(repmat([-1 0 +1],(J+1)^2,1),[-1+J+1 0 1+J+1],(J+1)^2,(J+1)^2);
            A1k = A1k + spdiags(repmat([1 0 -1],(J+1)^2,1),[-1-J-1 0 1-J-1],(J+1)^2,(J+1)^2);
            A1k = 1/(4*h^2)/(pi^2)*(Dcos_1.^2)*((Dcos.^2)*A1k);
        end
        
        function [A1,Ak,Am] = matrices3Dim(J)
            h = 1/J;
            col1 = repmat([-1/(2*h) 0 +1/(2*h)],(J+1),1);
            col1(J+1,1) = 0; col1(1,3) = 0;
            col2 = repmat([1/(h^2) -2/(h^2) 1/(h^2)],(J+1),1);
            col2(J+1,1) = 0; col2(1,3) = 0;
            % The zeros look like they are in the wrong spots, but this is correct due
            % to how spdiags works
            
            % Matrix for 3-dim case, 1st direction
            A1_1 = spdiags(repmat(col1,(J+1)^2,1),[-1 0 1],(J+1)^3,(J+1)^3);
            A1_2 = spdiags(repmat(col2,(J+1)^2,1),[-1 0 1],(J+1)^3,(J+1)^3);
            vcos = zeros(J+1,1); vsin = zeros(J+1,1);
            for j=1:J+1
                ind = j;
                vcos(ind) = (cos(pi*((j-1)*h-1/2)));
                vsin(ind) = (sin(pi*((j-1)*h-1/2)));
            end
            vcos = repmat(vcos,(J+1)^2,1); vsin = repmat(vsin,(J+1)^2,1);
            Dcos = spdiags(vcos,0,(J+1)^3,(J+1)^3); Dsin = spdiags(vsin,0,(J+1)^3,(J+1)^3);
            A1 = 1/(2*pi)*(Dcos.^3)*(-2*Dsin*A1_1+1/pi*Dcos*A1_2);
            clear A1_1; clear A1_2;
            
            % Matrix for 3-dim case, 2nd direction k
            Ak_1 = spdiags(repmat(col1,(J+1)^2,1),[-(J+1) 0 (J+1)],(J+1)^3,(J+1)^3);
            Ak_2 = spdiags(repmat(col2,(J+1)^2,1),[-(J+1) 0 (J+1)],(J+1)^3,(J+1)^3);
            vcos = zeros((J+1)^2,1); vsin = zeros((J+1)^2,1);
            for j=1:J+1
                ind = (1:(J+1)) + (j-1)*(J+1);
                vcos(ind) = (cos(pi*((j-1)*h-1/2)));
                vsin(ind) = (sin(pi*((j-1)*h-1/2)));
            end
            vcos = repmat(vcos,(J+1),1); vsin = repmat(vsin,(J+1),1);
            Dcos = spdiags(vcos,0,(J+1)^3,(J+1)^3); Dsin = spdiags(vsin,0,(J+1)^3,(J+1)^3);
            Ak = 1/(2*pi)*(Dcos.^3)*(-2*Dsin*Ak_1+1/pi*Dcos*Ak_2);
            clear Ak_1; clear Ak_2;
            
            % Matrix for 3-dim case, 3rd direction m
            Am_1 = spdiags(repmat(col1,(J+1)^2,1),[-(J+1)^2 0 (J+1)^2],(J+1)^3,(J+1)^3);
            Am_2 = spdiags(repmat(col2,(J+1)^2,1),[-(J+1)^2 0 (J+1)^2],(J+1)^3,(J+1)^3);
            vcos = zeros((J+1)^3,1); vsin = zeros((J+1)^3,1);
            for j=1:J+1
                ind = (1:((J+1)^2)) + (j-1)*(J+1)^2;
                vcos(ind) = (cos(pi*((j-1)*h-1/2)));
                vsin(ind) = (sin(pi*((j-1)*h-1/2)));
            end
            Dcos = spdiags(vcos,0,(J+1)^3,(J+1)^3); Dsin = spdiags(vsin,0,(J+1)^3,(J+1)^3);
            Am = 1/(2*pi)*(Dcos.^3)*(-2*Dsin*Am_1+1/pi*Dcos*Am_2);
            clear Am_1; clear Am_2;
            
        end
    end
    
end

