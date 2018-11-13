function [u]=mgv(u0,rhs,nu1,nu2,level,max_level)
n = length(rhs);
N = n -1;
m = (n+1)/2;
if (level == max_level)
[u,resvec,iter] = my_cg(u0,rhs,1.0e-8,1000);
else
u  = jacobi(u0,rhs,2/3,nu1);
rh = residual(u,rhs);
r2h = restriction(rh);
e2h = mgv(zeros(m,m),r2h,nu1,nu2,level+1,max_level);
%Interpolating on the finer grid
 eh = interpolation(e2h);
 u = u + eh;
 u = jacobi(u,rhs,2/3,nu2); 
end