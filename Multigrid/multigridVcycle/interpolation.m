function eh = interpolation(e2h);
            n = length(e2h);
            m = (n+1)/2;
            ind = 3:2:n-2;
            eh = zeros(n,n);
            eh(ind ,ind ) = e2h(2:m-1,2:m-1);
            eh(2:2:n-1,ind ) = .5*(e2h(1:m-1,2:m-1)+e2h(2:m,2:m-1));
            eh(ind ,2:2:n-1) = .5*(e2h(2:m-1,1:m-1)+e2h(2:m-1,2:m));
            eh(2:2:n-1,2:2:n-1) = .25*( ...
                                        e2h(1:m-1,1:m-1) + ...
                                        e2h(2:m,1:m-1) + ...
                                        e2h(1:m-1,2:m) + ...
                                        e2h(2:m,2:m));
end