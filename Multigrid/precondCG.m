function [u res it] = precondCG(u0,rhs,N,tol,maxit,nu1,nu2,X,Y)
    h = 1/N;
    b  = rhs;
    x0 = u0;
    x = x0;
    r0 = zeros(N+1);
    id = 2:N;
    %PCG starts here
    Ax0 = 1/h^2.*matvec(x0,N);
    r0(id,id) = b(id,id) - Ax0(id,id);
    z0 = mgv(zeros(N+1,N+1),r0,N,nu1,nu2,1,4); 
    z = z0;
    r = r0;
    p0 = z0;
    p = p0;
    error = 1; 
    iter = 0;
    subplot(2,3,1)
    surf(u0)
    title('Intial Solution')
    while error>tol   || iter == maxit
        Ap0   = 1/h^2.*matvec(p0,N);
        alpha = sum(dot(r0,z0))/sum(dot(Ap0,p0));  %alpha=(r0,r0)/(Ap0,p0)
        x(id,id) = x0(id,id) + alpha*p0(id,id);    %x = x0 + alpha*p0
        r(id,id) = r0(id,id) - alpha*Ap0(id,id);   %r = r0 - alpha*Ap0
        z = mgv(zeros(N+1,N+1),r,N,nu1,nu2,1,4); 
        beta = sum(dot(r,z))/sum(dot(r0,z0));      %beta = (r,r)/(r0,r0)
        p(id,id) = z(id,id) + beta*p0(id,id);      %p= r + beta*p0
        Ax = 1/h^2.*matvec(x,N);
        error = norm(b(id,id) - Ax(id,id))/norm(b(id,id) - Ax0(id,id));
        x0 = x;
        r0 = r;
        p0 = p;
        z0 = z;
        iter = iter + 1;
        rd(iter) = norm(b(id,id) - Ax(id,id));
        
        if iter<6
                subplot(2,3,iter+1)
                surf(X,Y,x)
                t = sprintf('PCG Iteration:%d',iter);
                title (t)
        end
    end
    it = iter;
    res = rd;
    u= x;     
end