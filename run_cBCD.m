clc,clear all
% generate source and target sets
% randn('state',111)
% rand('state',111)
% Ng = 35*1; % number of points in each group
% X = [randn(2,Ng)+repmat([-5;0],1,Ng) randn(2,Ng)+repmat([5;0],1,Ng) randn(2,Ng)+repmat([0;5],1,Ng)]; % source set
% X = X';
% [N, d] = size(X);
% 
% n = 3;
% num = [Ng, Ng, Ng];
% G = zeros(n,N);
% for i =1:1:n
%     if i == 1
%         G(i,1:num(i)) = 1/num(i);
%     else
%         G(i, sum(num(1:i-1))+1 : sum(num(1:i))) = 1/num(i);
%     end
% end
% 
% label = [1,2,3];
% l = max(label);
% Y = zeros(n,l);
% for i =1:1:n
%     Y(i,label(i))=1;
% end
%     
% c = 3;%floor(N/2); % parameter setting

load('moon.mat')
X = X';
[N, d] = size(X);
num = ones(1,10)*20;
n = length(num);
G = zeros(n,N);

%% Construct G
%First 
for i =1:1:n
    if i == 1
        G(i,1:num(i)) = 1/num(i);
    else
        G(i, sum(num(1:i-1))+1 : sum(num(1:i))) = 1/num(i);
    end
end
% %Second
% [idx,C] = kmeans(X,n);
% for i=1:1:n
%     G(i,find(idx==i)) = 1/length(find(idx==i));
% end

label = [ones(1,5),ones(1,5)*2];
l = max(label);
Y = zeros(n,l);
for i =1:1:n
    Y(i,label(i))=1;
end
c = 2;

%% initialize Z M
PRE = ones(N,c);
index = randperm(N);
ind = index(1:c);
for i = 1:1:c
    PRE(ind(i),i) = 1;
end
Z0 = PRE;
M0 = Z0;
% Z = Z0;
% M = M0;
[N,c] = size(Z0);
lambda1 = 2;
lambda2 = 0.5;
mu = 0.1;
Lambda = zeros(N,c);


A = G*X*X'; % SIZE: n*N

verbose = true; % true/false: show/hide optimization steps
k = 2;
terminate = false;
maxIter = 5;  % setting
errThr = 10^-1; % setting
Obj(1) = 10^1;

tic
%% Optimize Z M W
while (~terminate)
          W   =  derivative_solver(lambda2, G,X,M0, Y, c); % updating W
         [Z,M] = admm_solver(X, G, Z0, M0, W, Y, lambda1, mu, c);
        
        Obj(k) = (norm(A*Z*W-Y,'fro'))^2+ lambda1*Regul_sparsity(Z) + lambda2* (norm(W,'fro'))^2;
        err = abs(Obj(k)-Obj(k-1));
        
        if ( k >= maxIter || err <= errThr )
             terminate = true;
            if (verbose)
                fprintf('Terminating: \n');
                fprintf('||err_k-err_k-1||= %1.2e,iteration = %.0f \n\n',err, k);
            end
        else
           k = k + 1;
           if (verbose)
             if (mod(k,2)==0)
                fprintf('||err_k-err_k-1||= %1.2e,iteration = %.0f \n\n',err, k);
             end
           end
        end
       Z0 = Z;
       M0 = M;
end
toc
 plot (Obj)