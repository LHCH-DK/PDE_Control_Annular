function y=Precond_2D_Fou(x,CA,CB,C0,C1,C2,D,N,M)
%% Matrix-free inversion of the preconditioner
X1=x(1:(N-1)*M);
X2=x((N-1)*M+1:end);
X1=reshape(X1,N-1,M);
X2=reshape(X2,N-1,M);
for j=0:M-1
K=(j-M/2)^2;
for k=0:N-2
lambda=D(k+1,k+1);
sigma=CA+(CB*K+C0)*lambda;
kappa=(-sigma^2/(C2*lambda)-C1*lambda);
P(k+1,j+1)=kappa^(-1)*(X2(k+1,j+1)-(sigma*X1(k+1,j+1))/(C2*lambda));
Y(k+1,j+1)=(X1(k+1,j+1)-sigma*P(k+1,j+1))/(C2*lambda);
end
end
y=[Y(:);P(:)];
end