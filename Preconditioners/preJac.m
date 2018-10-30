function [U1,error]=preJac(A,F)
A1=diag(diag(A));
A2=A1-A;
U0=zeros(length(A),1);
for k=1:5
    U1=mldivide(A1,(A2*U0 + F));
    U0=U1;
    error(k)=norm((A1-A2)*U1 - F);
end
end