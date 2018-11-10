% updating Z
% min_{Z} lambda1 ||Z'||_{1,2}^{2}+mu/2 ||Z-M||_{F}^{2}+ <Lambda, Z-M>
% Input: lambda1, Lambda, mu, M
% Initiate: Z0, k
% Output: Z (size:N*c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z, F] = reweighted_solver(M,Lambda,lambda1,mu,Z0)

maxIter = 60;  % setting
errThr = 10^-10; % setting
elps = 10^-6;
[N,c] = size(M); % size of Z and M is same
q = zeros(1,c);
A = mu*M-Lambda;
verbose = false; %true; % true/false: show/hide optimization steps

for i = 1:1:N
    k = 1;
    terminate = false;
    m = M(i,:); % update i-th row
    a = A(i,:);
    z = Z0(i,:); % the i-th row of Z
    l = Lambda(i,:);
    
    while (~terminate)
        z_old = z;
        for j = 1:1:c
        q(j) = norm(z,1)/(abs(z(j))+elps)-1;
        end
        Q = diag(q);
        z = a * inv(mu*eye(c,c)+2*lambda1*Q);
        
        err = errorCoef_vector(z_old,z); % convergence condition
        
        f(k) = lambda1*z*Q*z'+mu/2*(norm(z-m+l/mu,2))^2;
        
        if ( k >= maxIter || (err <= errThr) )
             terminate = true;
            if (verbose)
                fprintf('Terminating: \n');
                fprintf('||z_old-z||= %1.2e,iteration = %.0f, row_index = %.0f \n\n',err, k, i);
            end
        else
           k = k + 1;
           if (verbose)
             if (mod(k,2)==0)
                fprintf('||z_old-z||= %1.2e,iteration = %.0f \n\n',err, k);
             end
           end
        end
    end
    Z(i,:) = z;
    F{i} = f;
end
        
        
        

        
        
    
