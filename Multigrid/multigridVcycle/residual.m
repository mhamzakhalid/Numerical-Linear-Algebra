function rh =residual(U,rhs)
    n = length(rhs);
    ind = 2:n-1;
    rh = zeros(n,n);
    Au = matvec(U,rhs);
    rh(ind,ind) = rhs(ind,ind) -  Au(ind,ind);
end                      
                        