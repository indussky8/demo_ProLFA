% This is the main function to run the proposed algorithm
% input:
% X: data set with   N*d
% Y: label set with  n*l
% G: group information with  n*N 
% output:
% Z: selection matrix with  N*c
% W: projection matrix with c*l
% initialize
% Z0
% W0
% M0: relaxtion matrix with N*c
% min  ||GXX'ZW-Y||_{F}^{2} + lambda1 ||Z'||_{1,2}^{2} + lambda2 ||W||_{F}^2
% s.t. Z>=0, 1^T Z = 1^T
%--------------------------------------------------------------------------
% Copyright @ Xing, 2018
%--------------------------------------------------------------------------

% % clc, clear all
% % %% DATA CONSTRUCT
% % load image_0001_sift.mat
% % Data = features.data(1:20,:);
% % num(1) = size(features.data,1);
% % load image_0002_sift.mat
% % Data = [Data;features.data(1:20,:)];
% % num(2) = size(features.data,1);
% % load image_0003_sift.mat
% % Data = [Data;features.data(1:20,:)];
% % num(3) = size(features.data,1);
% % load image_0004_sift.mat
% % Data = [Data;features.data(1:20,:)];
% % num(4) = size(features.data,1);
% % 
% % num = [200,200,200,200];
% % 
% % n = length(num);
% % X = Data;
% % [N,d] = size(X);
% % G = zeros(n,N);
% % for i =1:1:n
% %     if i == 1
% %         G(i,1:num(i)) = 1/num(i);
% %     else
% %         G(i, sum(num(1:i-1))+1 : sum(num(1:i))) = 1/num(i);
% %     end
% % end
% % 
% % Y = [1,0;1,0;0,1;0,1];
% % 
% % %% PARAMETER 
% % [N,d] = size(X);
% % [n,l] = size(Y);
% % c = floor(N/2); % parameter setting

clc,clear all
% generate source and target sets
randn('state',111)
rand('state',111)
Ng = 35; % number of points in each group
X = [randn(2,Ng)+repmat([-5;0],1,Ng) randn(2,Ng)+repmat([5;0],1,Ng) randn(2,Ng)+repmat([0;5],1,Ng)]; % source set
X = X';
[N, d] = size(X);

n = 3;
num = [35, 35, 35];
G = zeros(n,N);
for i =1:1:n
    if i == 1
        G(i,1:num(i)) = 1/num(i);
    else
        G(i, sum(num(1:i-1))+1 : sum(num(1:i))) = 1/num(i);
    end
end

label = [1,2,3];
l = max(label);
Y = zeros(n,l);
for i =1:1:n
    Y(i,label(i))=1;
end
    
    
c = 3;%floor(N/2); % parameter setting

%% initialize Z M
PRE = ones(N,c);
index = randperm(N);
ind = index(1:c);
for i = 1:1:c
    PRE(ind(i),i) = 1;
end
Z0 = PRE;
M0 = Z0;
Z = Z0;
M = M0;
[N,c] = size(Z);
lambda1 = 10;
lambda2 = 10^-1;
mu = 10^-1;
Lambda = zeros(N,c);


A = G*X*X'; % SIZE: n*N

verbose = true; % true/false: show/hide optimization steps
k = 1;
terminate = false;
maxIter = 40;  % setting
errThr = 10^-10; % setting
    
%% Optimize Z M W
while (~terminate)
        
        [Z, F] = reweighted_solver(M,Lambda,lambda1,mu,Z0); % updating Z
         W   =  derivative_solver(lambda2, G,X,M, Y, c); % updating W
        [M, f] = Fista_solver(Lambda, mu, A, W, Z, M0, Y) ; % updating M
         
        Lambda = Lambda + mu .* (Z - M);
        
        err1 = errorCoef(Z0,Z); % convergence condition
        err2 = errorCoef(M0,M); % convergence condition
        err3 = errorCoef(Z,M); % convergence condition
        err4 = errorCoef(ones(1,N)*Z,ones(1,c)); % convergence condition
        
        Obj(k) = (norm(A*Z*W-Y,'fro'))^2+ lambda1*Regul_sparsity(Z) + lambda2* (norm(W,'fro'))^2;
        
        if ( k >= maxIter || ((err1 <= errThr) &(err2 <= errThr) & (err3 <= errThr) &(err4 <= errThr)) )
             terminate = true;
            if (verbose)
                fprintf('Terminating: \n');
                fprintf('||Z_k-Z||= %1.2e,||M_k-M||= %1.2e, ||Z-M||= %1.2e,||1Z-1||= %1.2e,iteration = %.0f \n\n',err1, err2, err3, err4, k);
            end
        else
           k = k + 1;
           if (verbose)
             if (mod(k,2)==0)
                fprintf('||Z_k-Z||= %1.2e,||M_k-M||= %1.2e, ||Z-M||= %1.2e,||1Z-1||= %1.2e,iteration = %.0f \n\n',err1, err2, err3, err4, k);
             end
           end
        end
        fun_Z{k} = F;
        fun_M{k} = f;
        Z0 = Z;
        M0 = M;
end
 



