function U  = jacobi(U,rhs,omega,nu1)
        N = length(rhs) -1;
        ind = 2:N;
        h=1/N;
        for smoothing=1:nu1
            U(ind,ind) = 0.25*omega*( U(ind+1,ind) + U(ind-1,ind) +...
                                       U(ind,ind+1) + U(ind,ind-1) +...
                                       rhs(ind,ind)) + (1 - omega)*U(ind,ind);                                        
            %error(smoothing) = norm(residual(U,rhs,N));
        end   
end