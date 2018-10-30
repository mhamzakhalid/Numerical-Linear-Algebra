function x=fastpoissonPCG(b,n)
%T is (m by m) tridiag(-1,2,-1)
%x0 is the initial guess of size (m by m)
%b is the right hand side vector of size (m by m)

%Constructing T
%e = ones(n,1);
%T = spdiags([-e 2*e -e], -1:1, n, n);
T=full(gallery('tridiag',n,-1,2,-1));
%Reshaping vectors to matrix

bmat = reshape(b,[n,n]);
h=1/(n-1);

%Finding the eigen values
[Q,D]=eig(T);

%invD = diag(diag(1./D));

%Right Hand Side
bh = h^2*bmat;
G = Q'*bh*Q;
x1 = zeros(n);
lambda=diag(D);

for i=1:n
    for j=1:n
         x1(i,j)= G(i,j)/(lambda(i) + lambda(j));
    end
end

%End result 
xbar = Q*x1*Q';

% error(n) = norm(invD*xbar + xbar*invD - G);
%reshaping matrix to vector
 x = reshape(xbar,[n*n,1]);
%x = reshape(xbar',[],1);
end

