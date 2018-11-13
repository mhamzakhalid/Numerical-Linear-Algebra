function [U,res] = jacobi(U,rhs,omega,N,nu1)

res = zeros(N+1);
ind = 2:N;
h=1/N;
for smoothing=1:nu1
    U(ind,ind) = 0.25*omega*( U(ind+1,ind) + U(ind-1,ind) +...
                               U(ind,ind+1) + U(ind,ind-1) +...
                               h^2.*rhs(ind,ind)) + (1 - omega)*U(ind,ind);                                        
    Au = 1/h^2.*matvec(U,N);
    error = norm(rhs(ind,ind) - Au(ind,ind));
end   

  res(ind,ind) = rhs(ind,ind) - Au(ind,ind);
end