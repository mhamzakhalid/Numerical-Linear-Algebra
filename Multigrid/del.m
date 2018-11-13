clear
clc
%Initial values
N =40;
h = 1/N;
U = zeros(N+1);
tol =  1e-8;
maxit = 100;
%grid
x = 0:h:1;
y = 0:h:1;
[X,Y] = meshgrid(x,y);
%Defining Right Hand side and boundary
f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y);
g = @(x,y) sin(2*pi*x).*cos(4*pi*y);

%Boundary on U

U(:,1) = g(X(:,1),Y(:,1));
U(:,end) = g(X(:,end),Y(:,end));
U(1,:) = g(X(1,:)',Y(1,:)');
U(end,:) = g(X(end,:)',Y(end,:)');

%Right Hand side
rhs = h^2*f(X,Y);

%Implementing cg
[U,res,iter]=my_cg(U,rhs,N,tol,maxit);
mesh(U)
fprintf('Converged in %d iterations\nResidual:%e',iter,res)
    

   