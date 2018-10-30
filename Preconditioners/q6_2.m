clear
clf
%N=44;
N=19;
dy=1/(N-1);
dx=1/(N-1); %step-length                                                   %Is dx step or some other factor as from question?
%k=sqrt(8/dx^2)+1; %Parameter in Helmotz Equation
k=1e2;
x=0:dx:1;
y=0:dy:1;

for i=1:N
    for j=1:N
        F(i,j)=exp(-50*((x(j) - 1/2)^2))*exp(-50*(y(i) - 1/2)^2);
    end 
end 
F = reshape(F,[N*N,1]);  % Right Hand side
%L-Matrix
e = ones(N*N,1);
f = ones(N*N,1);
f(N:N:N*N) = 0;
g = [1;f(1:end-1)];
Lap = spdiags([e,f, -4*e, g, e],[-N,-1,0,1,N],N*N,N*N);

%A-Matrix 
A = Lap + dx^2*k^2*eye(N*N);
A=sparse(A);

%Check Positive Definitness 
if eig(A)>0
    %Postive Definite Case
    disp('A is Positive Definite Matrix')
else
    %A not Positive Definite 
    disp('Not PDM')
end
tol = 1e-6;
disp(' ')
Upcg = pcg(A,F,1e-14);

%% Incomplete Cholesky factorization
tic
L = ichol(A);
Mcholesky = L*L';

U=zeros(N*N,1);
U=reshape(U,[N*N,1]);
U0=U;
[U,icholesky,err_chol] = preCholesky(U0,A,Mcholesky,F,tol);
 tpcg=toc;
fprintf('\nIncomplete Cholesky Factorization\nTime: %f\nIterations:%d\n\n',tpcg,icholesky);
A_nchol = Mcholesky\A;
estimatechol = (sqrt(cond(full(A_nchol))) - 1)/(sqrt(cond(full(A_nchol))+1));
for iter=1:6
       est_chol(iter) = 2*estimatechol^(iter-1)*sqrt(dot((A_nchol)*Upcg,Upcg));
end
semilogy(est_chol,'--')
hold on
semilogy(err_chol)
hold on
grid on

%% Fast Poisson Solvers
tic
U=zeros(N*N,1);
U=reshape(U,N*N,[]);
U0=U;
%PCG
r0 = F - A*U0;
err_fastp(1)=norm(r0);
z0=zeros(length(r0),1);
%Using Fast Poisson Solver
z0=fastpoissonPCG(r0,N); 
p0 = z0;
error=1.;
iter=0;
    while error>tol 
        Ap0 = A*p0;
        alpha = dot(r0,z0)/dot(Ap0,p0);
        U1 = U0 + alpha*p0;
        r1 = r0 - alpha*Ap0;
        %Applying Fast Poisson
                   
         z1=fastpoissonPCG(r1,N);   
         
        beta = dot(r1,z1)/dot(r0,z0);
        p1 = z1 + beta*p0;
        p0 = p1;
        r0 = r1;
        z0 = z1;
        U0 = U1;
        %error = norm(U1 - Upcg);
        error = norm(F - A*U0);
        iter = iter +1;
        err_fastp(iter+1) = error;       
    end
tfastP = toc;
fprintf('Fast Poisson Solver\nTime: %f\nIterations:%d\n\n',tfastP,iter)
Mfastp = -Lap;
A_fastp = Mfastp\A;
estimatefastp = (sqrt(cond(full(A_fastp))) - 1)/(sqrt(cond(full(A_fastp))+1));
for it=1:6
       est_fastp(it) = 2*estimatefastp^(it-1)*sqrt(dot((A_fastp)*Upcg,Upcg));
end
semilogy(est_fastp,'--')
hold on
semilogy(err_fastp)
hold on

 %% Jacobi
tic
U=zeros(N*N,1);
U=reshape(U,N*N,[]);
U0=U;
%PCG
r0 = F - A*U0;
err_jac(1)=norm(r0);
z0=zeros(length(r0),1);

[z0,er]=preJac(A,r0);

p0 = z0;
error=1.;
iter=0;
    while error>tol 
        Ap0 = A*p0;
        alpha = dot(r0,z0)/dot(Ap0,p0);
        U1 = U0 + alpha*p0;
        r1 = r0 - alpha*Ap0;
        %Applying Jacobi iterations for M=diagonal
        [z1,ejac]= preJac(A,r1);
        %z1 = M_jcb\r1;            
             
        beta = dot(r1,z1)/dot(r0,z0);
        p1 = z1 + beta*p0;
        p0 = p1;
        r0 = r1;
        z0 = z1;
        U0 = U1;
        %error = norm(x1 - Upcg);
        error = norm(F - A*U0);
        iter = iter +1;
        err_jac(iter+1) = error;
       % err_est(iter) = 2*estimate^(1*iter)*sqrt(Upcg'*(M\A)*Upcg);
        
    end
tjac = toc;
fprintf('Jacobi Method\nTime: %f\nIterations:%d\n\n',tjac,iter)
Mjac = diag(diag(A));
A_jac = Mjac\A;
estimatejac = (sqrt(cond(full(A_jac))) - 1)/(sqrt(cond(full(A_jac))+1));
for it=1:6
       est_jac(it) = 2*estimatejac^(it-1)*sqrt(dot((A_jac)*Upcg,Upcg));
end
semilogy(est_jac,'--')
hold on
semilogy(err_jac)
hold on


%% Symmetric Gauss-Seidal
 tic
U=zeros(N*N,1);
U=reshape(U,N*N,[]);
U0=U;
%PCG
r0 = F - A*U0;
err_symGS(1)=norm(r0);
z0=zeros(length(r0),1);
%Applying Symmetric Gauss Seidal
z0=symGS(A,r0,tol); 
p0 = z0;
error=1.;
iter=0;
    while error>tol 
        Ap0 = A*p0;
        alpha = dot(r0,z0)/dot(Ap0,p0);
        U1 = U0 + alpha*p0;
        r1 = r0 - alpha*Ap0;
        %Applying Symmetric Gauss Seidal iterations 
                    
         z1=symGS(A,r1,tol);   
        beta = dot(r1,z1)/dot(r0,z0);
        p1 = z1 + beta*p0;
        p0 = p1;
        r0 = r1;
        z0 = z1;
        U0 = U1;
       
        error = norm(F - A*U0);
        iter = iter +1;
        err_symGS(iter+1) = error;
       % err_est(iter) = 2*estimate^(1*iter)*sqrt(Upcg'*(M\A)*Upcg);
        
    end
    tsymGS = toc;
    fprintf('Symmetric GS Method\nTime: %f\nIterations:%d\n\n',tsymGS,iter)
    D = diag(diag(A));
    E = -triu(A,1);
    MsymGS= (D - E);
A_SymGS = MsymGS\A;
estimatesymGS = (sqrt(cond(full(A_SymGS))) - 1)/(sqrt(cond(full(A_SymGS))+1));
for it=1:6
       est_symGS(it) = 2*estimatesymGS^(it-1)*sqrt(dot((A_SymGS)*Upcg,Upcg));
end
semilogy(est_symGS,'--')
hold on
semilogy(err_symGS)

title('Performance of different preconditioners')
xlabel('Iteration Number')
ylabel('Error')
hold on
legend('Estimate ICF','Incomplete Cholesky Factorization','Estimate FPS','Fast Poisson Solver','Estimate Jacobi','Jacobi Method','Estimate SGS','Symmetric GS')

