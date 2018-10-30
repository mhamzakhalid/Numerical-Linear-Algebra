function [x1,iter,err]=preCholesky(x0,A,M,b,tol)
Upcg = pcg(A,b,tol,[]);
r0 = b - A*x0;
err(1)=norm(r0); 
z0 = M\r0;
p0 = z0;
error=1.;
iter=0;
    while error>tol 
        Ap0 = A*p0;
        alpha = dot(r0,z0)/dot(Ap0,p0);
        x1 = x0 + alpha*p0;
        r1 = r0 - alpha*Ap0;
        z1 = M\r1;
        beta = dot(r1,z1)/dot(r0,z0);
        p1 = z1 + beta*p0;
        p0 = p1;
        r0 = r1;
        z0 = z1;
        x0 = x1;
        error = sqrt(dot(A*(Upcg - x0),(Upcg - x0)));
        %error = norm(b - A*x0);
        iter = iter +1;
        err(iter+1) = error;
             
    end
       
end