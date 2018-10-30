function  [x1,i] = step_desc(x,A,b,tol)
r0 = b - A*x;
p0=A*r0;
error = 1.;
i=0;
    while error>tol
           alpha = dot(p0,r0)/dot(p0,p0);
           x1 = x + alpha*r0;
           r = r0 - alpha*p0;
           p = A_rf*r;
           error = norm(b - A*x1,2);
           x = x1;
           p0 = p;
           r0 = r;
           i=i+1;
    end

end