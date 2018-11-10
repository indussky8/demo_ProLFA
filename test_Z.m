% test reweighted_solver function (updating Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear
N = 10;
c = 6;
M = eye(N,N);
M = M(:,1:c);
lambda1 = 10^-1;
mu = 10^-1;
Z0 = M;
Lambda = zeros(N,c);

[Z, F] = reweighted_solver(M,Lambda,lambda1,mu,Z0);