function u = jacobi(u0,rhs,omega,N,nu1)

index = 2:N;
h=1/N;
tol=1e-8;
while error>tol 
    u1(index,index) = 0.25*(u0(index,index-1)+ u0(index-1,index) + ...
                            u0(index,index+1)+ u0(index+1,index)+...
                            h^2*rhs);
                   error=norm()
end