% This is ADMM solver to obtain Z 
% input:
% X: data set with   N*d
% Y: label set with  n*l
% G: group information with  n*N 
% W: projection matrix with c*l
% output:
% Z: selection matrix with  N*c
% M: selection matrix with  N*c
% initialize
% Z0
% M0: relaxtion matrix with N*c
% min  ||GXX'ZW-Y||_{F}^{2} + lambda1 ||Z'||_{1,2}^{2} 
% s.t. Z>=0, 1^T Z = 1^T
%--------------------------------------------------------------------------
% Copyright @ Xing, 2018
%--------------------------------------------------------------------------

function [Z,M] = admm_solver(X, G, Z0, M0, W, Y, lambda1, mu, c)

n = size(Y,1);  %n samples
[N,d] = size(X); % N descriptors 

% X:N*d  % G:n*N  Z:N*c 
%% PARAMETER 
[N,d] = size(X);
[n,l] = size(Y);

%% initialize Z M
Z = Z0;
M = M0;
[N,c] = size(Z);


Lambda = zeros(N,c);


A = G*X*X'; % SIZE: n*N

verbose = true; % true/false: show/hide optimization steps
k = 1;
terminate = false;
maxIter = 40;  % setting
errThr = 10^-9; % setting
    
%% Optimize Z M W
while (~terminate)
        
        [Z, F] = reweighted_solver(M,Lambda,lambda1,mu,Z0); % updating Z
        [M, f] = Fista_solver(Lambda, mu, A, W, Z, M0, Y) ; % updating M
         
        Lambda = Lambda + mu .* (Z - M);
        
        err1 = errorCoef(Z0,Z); % convergence condition
        err2 = errorCoef(M0,M); % convergence condition
        err3 = errorCoef(Z,M); % convergence condition
        err4 = errorCoef(ones(1,N)*Z,ones(1,c)); % convergence condition
        
        Obj(k) = (norm(A*Z*W-Y,'fro'))^2+ lambda1*Regul_sparsity(Z);
        
        if ( k >= maxIter || ((err1 <= errThr) &(err2 <= errThr) & (err3 <= errThr) &(err4 <= errThr)) )
             terminate = true;
            if (verbose)
                fprintf('Terminating: \n');
                fprintf('||Z0-Z||= %1.2e,||M0-M||= %1.2e, ||Z-M||= %1.2e,||1Z-1||= %1.2e,iteration = %.0f \n\n',err1, err2, err3, err4, k);
            end
        else
           k = k + 1;
           if (verbose)
             if (mod(k,2)==0)
                fprintf('||Z0-Z||= %1.2e,||M0-M||= %1.2e, ||Z-M||= %1.2e,||1Z-1||= %1.2e,iteration = %.0f \n\n',err1, err2, err3, err4, k);
             end
           end
        end
        fun_Z{k} = F;
        fun_M{k} = f;
        Z0 = Z;
        M0 = M;
end
 



