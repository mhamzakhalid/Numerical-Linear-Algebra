clear


%Initial values

%grids =[10,20,50];
%for gr = 1:length(grids)
    %N = grids(gr);
%for L=5:9    
    N = 2^5;
    h = 1/N;
    U = zeros(N+1);
    tol =  1e-12;
    maxit = 100;
    %grid
    x = 0:h:1;
    y = 0:h:1;
    [X,Y] = meshgrid(x,y);

    %% CG
    %%Defining Right Hand side and boundary
    % f = @(x,y) 20*pi^2*sin(2*pi*x).*cos(4*pi*y); %CG
    % g = @(x,y) sin(2*pi*x).*cos(4*pi*y); %CG
    % %Boundary on U
    % U(:,1) = g(X(:,1),Y(:,1));
    % U(:,end) = g(X(:,end),Y(:,end));
    % U(1,:) = g(X(1,:)',Y(1,:)');
    % U(end,:) = g(X(end,:)',Y(end,:)');
    %  
    % %Right Hand side
    % rhs = h^2*f(X,Y);
    %%Implementing cg
    %[U,res,iter]=my_cg(U,rhs,N,tol,maxit);
    %mesh(U)
    %fprintf('Converged in %d iterations\nResidual:%e',iter,res)


    %% MGV
    f = @(x,y) -1; 
    g = @(x,y) 4*y.*(1-y);
    %Boundary Conditions
    U(:,1) = g(X(:,1),Y(:,1));
    rnd = rand(N-1);
    U(2:N,2:N) = rnd(1:end,1:end);
    %%Right Hand side
    rhs = zeros(N+1);
    rhs(:,:) = h^2*f(X,Y);

    %% MGV
    nu1 = 3;
    nu2 = 3;
    max_level = 3;
    level = 1;
    e =1.;
    Au0 = matvec(U,N);
    %[U,res,iter]=my_cg(U,rhs,N,tol,maxit);
    %[U res] = jacobi(U,rhs,2/3,N,5);
     U1 = mgv2d(U,rhs,3);
     %mesh(U1)  
    U  = mgv(U,rhs,nu1,nu2,level,max_level);
    %figure
    %mesh(U)



