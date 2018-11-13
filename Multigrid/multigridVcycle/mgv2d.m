function Z = mgv2d(U,B,gs_iter)
% 2d multigrid v-cycle for ch5ex10

[n n1] = size(B);

if ( n~=n1 )
    error('mgv2d: input not square');
end

N = n-1;
h = 2/N;

if ( n==4)
    % final solution
    Z = zeros(n,n);
    Z(2,2) = -h^2/4*B(2,2);
else
    % k GS iterations
    for iter=1:gs_iter
        for j = 2:n-1
            for i = 2:n-1
                U(i,j) = (U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1))/4 - h^2/4*B(i,j);
            end
        end
    end
 
    % compute residual
    cen = 2:n-1;
    dec = cen - 1;
    inc = cen + 1;
    R = zeros(n,n);
    R(cen,cen) = B(cen,cen) - 1/h^2*( ...
        -4*U(cen,cen) ...
        +U(dec,cen) ...
        +U(inc,cen) ...
        +U(cen,dec) ...
        +U(cen,inc) );

    % restrict residual
    cen_s = 3:2:n-2;
    dec_s = cen_s - 1;
    inc_s = cen_s + 1;

    m = (n+1)/2;
    Rhat = zeros(m,m);
    Rhat(2:m-1,2:m-1) = 1/16*( ...
        R(dec_s,dec_s) + ...
        R(dec_s,inc_s) + ...
        R(inc_s,dec_s) + ...
        R(inc_s,inc_s) + ...
        2*R(inc_s,cen_s) + ...
        2*R(dec_s,cen_s) + ...
        2*R(cen_s,inc_s) + ...
        2*R(cen_s,dec_s) + ...
        4*R(cen_s,cen_s));

    % recur
    Uhat = mgv2d(zeros(m,m),Rhat,gs_iter);

    % interpolate
    Ucor = zeros(n,n);
    Ucor(cen_s,cen_s) = Uhat(2:m-1,2:m-1);
    Ucor(2:2:n-1,cen_s) = .5*(Uhat(1:m-1,2:m-1)+Uhat(2:m,2:m-1));
    Ucor(cen_s,2:2:n-1) = .5*(Uhat(2:m-1,1:m-1)+Uhat(2:m-1,2:m));
    Ucor(2:2:n-1,2:2:n-1) = .25*( ...
        Uhat(1:m-1,1:m-1) + ...
        Uhat(2:m,1:m-1) + ...
        Uhat(1:m-1,2:m) + ...
        Uhat(2:m,2:m));
    
    Z = U + Ucor;
    
    % k steps of GS
    for iter=1:gs_iter
        for j = 2:n-1
            for i = 2:n-1
                Z(i,j) = (Z(i+1,j)+Z(i-1,j)+Z(i,j+1)+Z(i,j-1))/4 - h^2/4*B(i,j);
            end
        end
    end
end


