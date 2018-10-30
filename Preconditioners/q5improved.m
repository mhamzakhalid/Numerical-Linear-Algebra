clear

%loading the Stokes problem
load StokesProblem

%Setting up the right hand side
%A-1 b = g  Solve Ag = b
tol = 1e-6;
g0 = zeros(length(A),1);
[g,gi] = conj_grad(g0,A,b,tol);
rhs = B'*g - c;
tic
%Setting up the LHS
    y0 = zeros(size(B,2),1); 
    r0 = rhs;
    p0 = r0;
    error = 1; 
    i = 0;
    Abar = B'*mldivide(A,B); %For error evaluation
    while error>tol     
        d1 = B*p0;
        [d2,di] = conj_grad(zeros(length(d1),1),A,d1,tol);
        Ap0 = B'*d2;     
        alpha = dot(r0,r0)/dot(Ap0,p0);
        y = y0 + alpha*p0;
        r = r0 - alpha*Ap0;
        beta = dot(r,r)/dot(r0,r0);
        p = r + beta*p0;          
        error = norm(rhs - Abar*y);
        y0 = y;
        r0 = r;
        p0 = p;
        i = i + 1;
        err(i) = error;
    end
    
    %Plots
    semilogy(err)
    grid on
    title('Nested CG applied to Saddle Point Problem')
    xlabel('Iterations')
    ylabel('Error')
    %Computing x
    b_x = b-B*y;
    [x,xi] = conj_grad(zeros(length(A),1),A,b-B*y,tol);
    u=[x;y];
    t=toc;
O = zeros(size(B,2));
A_or = [A B;B' O];
sparse(A_or);   

Y = gmres(Abar,rhs,[],1e-12);
X = gmres(A,b-B*Y,[],1e-12);

fprintf('\nNESTED CG \nNo of Iterations: %d \nTime: %f\nError: %e \n \n',i,t,error)
