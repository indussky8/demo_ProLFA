% TEST Fista_solver function for updating M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear
N = 10;
c = 4;
l = 2;
n = 6;
Z = eye(N,N);
Z = Z(:,1:c);
lambda1 = 10^-1;
mu = 10^-1;
M0 = Z;
Lambda = zeros(N,c);
W = rand(c,l);
Y = [1,0;1,0;1,0;0,1;0,1;0,1];
A = rand(n,N);
[M,f] = Fista_solver(Lambda, mu,A, W, Z, M0, Y);
plot(f)
