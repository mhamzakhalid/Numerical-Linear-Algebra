function  [u,res,iter]=my_cg(u0,rhs,tol,maxit)
   %Initializations
    N = length(rhs)-1;
    b  = rhs;
    x0 = u0;
    x = x0;
    r0 = zeros(N+1);
    id = 2:N;
    %CG starts here
    Ax0 = matvec(x0,N);
    r0(id,id) = b(id,id) - Ax0(id,id);
    r = r0;
    p0 = r0;
    p = p0;
    error = 1; 
    iter = 0;
    while error>tol   || iter == maxit
        Ap0   = matvec(p0,N);
        alpha = sum(dot(r0,r0))/sum(dot(Ap0,p0));  %alpha=(r0,r0)/(Ap0,p0)
        x(id,id) = x0(id,id) + alpha*p0(id,id);    %x = x0 + alpha*p0
        r(id,id) = r0(id,id) - alpha*Ap0(id,id);   %r = r0 - alpha*Ap0
        beta = sum(dot(r,r))/sum(dot(r0,r0));      %beta = (r,r)/(r0,r0)
        p(id,id) = r(id,id) + beta*p0(id,id);      %p= r + beta*p0
        Ax = matvec(x,N);
        error = norm(b(id,id) - Ax(id,id))/norm(b(id,id) - Ax0(id,id));
        x0 = x;
        r0 = r;
        p0 = p;
        iter = iter + 1;
        res(iter) = norm(b(id,id) - Ax(id,id));
    end
    
    u= x;     
end


 