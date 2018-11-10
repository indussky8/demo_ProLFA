% updating W
% min_{W} lambda2 ||W||_{F}^{2}+ ||GXX'MW-Y||_{F}^{2}
% Input: lambda2, GXX'M Y
% Output: W (c*l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W =  derivative_solver(lambda2, G,X,M, Y, c)
A = G*X*X'*M; % SIZE: n*c
[m ,n] = size(A);  

W = inv(A'*A+lambda2*eye(c))*A'*Y;

end