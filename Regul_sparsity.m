%The L_{1,2}norm of Z
function regu_Z = Regul_sparsity(Z)
[N,c] = size(Z);
for i = 1:1:N
    z(i) = norm(Z(i,:),1);
end
regu_Z = (norm(z,2))^2;