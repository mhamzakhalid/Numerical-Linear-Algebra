clear
clc

N=14;
dy=1/(N-1);
dx=1/(N-1); %step-length                                                   %Is dx step or some other factor as from question?
%k=sqrt(8/dx^2)+1; %Parameter in Helmotz Equation
k =1e2;
U=zeros(N*N,1);
U=reshape(U,[N*N,1]);
F=zeros(N,N);
U0=U;

x=0:dx:1;
y=0:dy:1;
for i=1:N
    for j=1:N
        F(i,j)=exp(-50*((x(j) - 1/2)^2 + (y(i) - 1/2)^2));
    end
end
F = reshape(F,[N*N,1]);  % Right Hand side

%L-Matrix
e = ones(N*N,1);
f = ones(N*N,1);
f(N:N:N*N) = 0;
g = [1;f(1:end-1)];
L = spdiags([e,f, -4*e, g, e],[-N,-1,0,1,N],N*N,N*N);

%A-Matrix 
A = L + dx^2*k^2*eye(N*N);
A=sparse(A);
%Check Positive Definitness 
if eig(A)>0
    %Postive Definite Case
    disp('A is Positive Definite Matrix')
    %Applying Conjugate Gradient Method
    x_ex = ones(length(A),1);
    b = A*x_ex;
    
    x0 = zeros(length(A),1);
    r0 = b - A*x0;
    p0 = r0;
    error = 1;
    tol = 1e-6;
    i = 0;
    estimate = (sqrt(condest(A)) - 1)/(sqrt(condest(A)+1));
    while error>tol     
        Ap0 = A*p0;
        alpha = dot(r0,r0)/dot(Ap0,p0);
        x = x0 + alpha*p0;
        r = r0 - alpha*Ap0;
        beta = dot(r,r)/dot(r0,r0);
        p = r + beta*p0;
        
        error = norm(b - A*x)/norm(b);
        
        err(i+1) = sqrt(dot(A*(x_ex - x),(x_ex - x)));
        
        res(i+1) = norm(b - A*x);
        x0 = x;
        r0 = r;
        p0 = p;
        err_est(i+1) = 2*estimate^i*sqrt(dot(A*x_ex,x_ex));
        i = i + 1;     
        
        
    end
    %Plots
  
    semilogy(err)
    hold on
    semilogy(res)
    hold on
    plot(err_est,'--')
    title('Error, Residual and Error estimate VS Iteration number')
    xlabel("Iteration No")
    ylabel("Y-axis")
    set(gca, 'FontName', 'Times New Roman')
    grid on 
    legend('Error','Residual','Error estimate')


else
    %A not Positive Definite 
    disp('Not PDM')
end

