% Sylvester equation
% updating W
% min_{W} lambda2 ||W||_{F}^{2}+ ||GXX'MW-Y||_{F}^{2}
% Input: lambda2, GXX'M Y
% Output: W (c*l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W =  Solver_W(lambda2, G,X,M, Y, c)
A = G*X*X'*M; % SIZE: n*c
[m, n] = size(A);  

W_A = A'*A+lambda2*eye(c);%d'*d'
W_B = Y'*Y; %c*c
W_Q = 2*A'*Y;
W = sylvester(W_A,W_B,W_Q);


