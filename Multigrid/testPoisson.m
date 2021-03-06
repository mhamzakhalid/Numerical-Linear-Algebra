clear 


%% TASK 1 
% task = 'CG';
%% TASK 2
% task = 'MGV';
%% TASK 3
 task = 'PCG';

%% Initial Values
L = 5;
N = 2^L;
h = 1/N;
tol = 1e-12;
maxit = 100;
nu1 = 3;
nu2 = 3;
level = 1;
e = 1.;
it = 0;

%% Programs
switch task
    case 'CG'
        f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); %CG
        g = @(x,y) sin(2*pi*x).*cos(4*pi*y); %CG
    case {'MGV', 'PCG'}
        f = @(x,y) -1; 
        g = @(x,y) 4*y.*(1-y);
    otherwise 
        disp('Enter Correct case')
end
% Loading Geometry
[U rhs X Y] = geometryPoisson(N,f,g,task);

switch task
    case 'CG'
        [U,res,iter]=my_cg(U,rhs,N,tol,maxit);
        surf(X,Y,U)
        ts = sprintf('Conjugate Gradient Solution for N: %d',N);
        title(ts)
        figure
        semilogy(res)
        tc = sprintf('Convergence CG for N: %d',N);
        title(tc)
        xlabel('Iterations')
        ylabel('Residual')
        grid on
    case 'MGV' 
        max_level = L;
        Au0 = matvec(U,N);

        subplot(2,3,1) 
        surf(X,Y,U)
        title ('Initial guess')
  
        Au0 = 1/h^2.*matvec(U,N);
       while e>tol  
        U  = mgv(U,rhs,N,nu1,nu2,level,max_level);             
        Au = 1/h^2.*matvec(U,N); 
        res(it + 1) = norm(rhs(2:N,2:N) - Au(2:N,2:N));
        e = norm(rhs(2:N,2:N) - Au(2:N,2:N))/norm(rhs(2:N,2:N) - Au0(2:N,2:N));
        it = it +1;
            if it<6
                subplot(2,3,it+1)
                surf(X,Y,U)
                t = sprintf('V-cycle Iteration:%d',it);
                title (t)
            end
         
       end
       figure
       semilogy(res)
       tc = sprintf('Convergence Multigrid V-cycle for N: %d',N);
       title(tc)
       xlabel('Iterations')
       ylabel('Residual')
       grid on
       
    case 'PCG'
        max_level = L;
        Au0 = matvec(U,N);

        subplot(2,3,1) 
        surf(X,Y,U)
        title ('Initial guess')
        
       [U res it]  = precondCG(U,rhs,N,tol,maxit,nu1,nu2,X,Y);             

       figure
       semilogy(res)
       tc = sprintf('Convergence PCG for N: %d',N);
       title(tc)
       xlabel('Iterations')
       ylabel('Residual')
       grid on
end
        

