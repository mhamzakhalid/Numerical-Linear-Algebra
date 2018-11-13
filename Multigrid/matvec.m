function Ax = matvec(x,N)
    Ax = x; 
    index = 2:N;
    Ax(index,index) = 4*x(index,index)-x(index+1,index)-x(index,index+1)-...
                        x(index-1,index)-x(index,index-1);
end
        