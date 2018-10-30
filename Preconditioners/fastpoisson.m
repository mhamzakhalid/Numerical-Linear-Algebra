clear


N=4;
dy=1/(N-1);
dx=1/(N-1); %step-length                                                   %Is dx step or some other factor as from question?
%k=sqrt(8/dx^2)+1; %Parameter in Helmotz Equation
k=1e2;
F=zeros(N,N);
x=0:dx:1;
y=0:dy:1;
% x =linspace(0,1,N)';
% y =linspace(0,1,N)';
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

tol = 1e-14;
disp(' ')
%exact


Upcg = pcg(A,F,tol);
%Preconditioned CG

U=zeros(N*N,1);
U=reshape(U,N*N,[]);
U0=U;
%PCG
r0 = F - A*U0;
err_jac(1)=norm(r0);
z0=zeros(length(r0),1);

%[z0,er]=preJac(A,r0);
%z0 = Lap\r0;
z0=fastpoissonPCG(r0,N); 
p0 = z0;
error=1.;
iter=0;
    while error>tol 
        Ap0 = A*p0;
        alpha = dot(r0,z0)/dot(Ap0,p0);
        U1 = U0 + alpha*p0;
        r1 = r0 - alpha*Ap0;
        %Applying Jacobi iterations for M=diagonal
        %[z1,ejac]= preJac(A,r1);
        %z1 = Lap\r1;            
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
       % err_est(iter) = 2*estimate^(1*iter)*sqrt(Upcg'*(M\A)*Upcg);
        
    end
semilogy(err_fastp)