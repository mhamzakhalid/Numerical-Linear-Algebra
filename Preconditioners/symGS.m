function [x]=symGS(A,b,tol)

D = diag(diag(A));
E = -triu(A,1);
F = -tril(A,-1);
x0 = zeros(length(A),1);
for k=1:5
    x12 = (D - E)\F*x0 + (D - E)\b;
    x1 = (D - F)\E*x12 + (D - F)\b;
    x0 = x1;
    error = norm(b - A*x1);
end
x = x1;
end
