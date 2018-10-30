function  [x,i]= conj_grad(x0,A,b,tol)
    r0 = b - A*x0;
    p0 = r0;
    error = 1; 
    i = 0;
    while error>tol   
        Ap0 = A*p0;
        alpha = dot(r0,r0)/dot(Ap0,p0);
        x = x0 + alpha*p0;
        r = r0 - alpha*Ap0;
        beta = dot(r,r)/dot(r0,r0);
        p = r + beta*p0;    
        error = norm(b - A*x)/norm(b);
        x0 = x;
        r0 = r;
        p0 = p;
        i = i + 1;
    end
end