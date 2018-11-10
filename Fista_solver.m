% updating M
% min_{M} ||GXX'MW-Y||_{F}^{2}+mu/2 ||Z-M||_{F}^{2}+ <Lambda, Z-M>
% Input: Lambda, mu, M, G(SIZE: n*N), X(SIZE: N*d), W, Z, Y
% Initiate: M0
% Output: M (size:N*c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function M = Fista_solver(Lambda, mu, G, X, W, Z, M0, Y)
function [M,f] = Fista_solver(Lambda, mu, A, W, Z, M0, Y)

%% Compute Lipschitz constant
% A = G*X*X'; % SIZE: n*N
B = 2* A'*A;
C = W*W';
[U1, Delta1] = eig(B); % B=U1*Delta1*V1';
[U2, Delta2] = eig(C); % C=U2*Delta2*V2'
c1 = (norm(inv(U1'),'fro'))^2;
c2 = (norm(inv(U2),'fro'))^2;
c3 = (norm(U1','fro'))^2;
c4 = (norm(U2,'fro'))^2;

delta1 = diag(Delta1);
delta2 = diag(Delta2);
Delta = delta1*delta2';
Constant = Delta+mu;
c5 = max(max(Constant))^2;

Lipschitz = c1*c2*c3*c4*c5; % Lipschitz constant

%% Initialize
M = M0;
S = M;
t = 1;

maxIter = 100;  % setting
errThr = 10^-10; % setting
[N,c] = size(Z); % size of Z and M is same
verbose = false;% true; % true/false: show/hide optimization steps
k = 1;
terminate = false;
    
%% Optimize M
while (~terminate)
        M_old = M;
        P = S-1/Lipschitz*(2*A'*(A*S*W-Y)*W'+mu*(S-Z-Lambda/mu));
        M  = solver_BCLS_closedForm(P);
        q = t-1;
        t = (1+sqrt(1+4*t^2))/2;
        S = M + q*(M-M_old)/t;
        
        err = errorCoef(M_old,M); % convergence condition
        
        f(k) = (norm(A*M*W-Y,'fro'))^2+ mu/2*(norm(Z-M+Lambda/mu,'fro'))^2;
        
        if ( k >= maxIter || (err <= errThr) )
             terminate = true;
            if (verbose)
                fprintf('Terminating: \n');
                fprintf('||M_old-M||= %1.2e,iteration = %.0f \n\n',err, k);
            end
        else
           k = k + 1;
           if (verbose)
             if (mod(k,2)==0)
                fprintf('||M_old-M||= %1.2e,iteration = %.0f \n\n',err, k);
             end
           end
        end
 end











