function [u]=mgv(u0,rhs,N,nu1,nu2,level,max_level)
%
% function [u] = mgv(u0,rhs,N,nu1,nu2,level,max_level) performs
% one multigrid V-cycle for the 2D Poisson problem on the unit
% square [0,1]x[0,1] with initial guess u0 and righthand side rhs.
%
% input: u0 - initial guess
% rhs - righthand side
% N - u0 and rhs are (N+1)x(N+1) matrices
% nu1 - number of presmoothings
% nu2 - number of postsmoothings
% level - current level
% max_level - total number of levels
%
if (level == max_level)
[u,resvec,iter] = my_cg(u0,rhs,N,1.0e-8,1000);
else
[u , rh] = jacobi(u0,rhs,2/3,N,nu1);
r2h = restriction(rh,N);
e2h = mgv(zeros(N/2+1),r2h,N/2,nu1,nu2,level+1,max_level);
%Interpolating on finer grid
 eh = interpolation(e2h,N);
 u = u + eh;
 u = jacobi(u,rhs,2/3,N,nu2); 
end