function P=JacobiPLegendre(x,N)

if N==0
P=ones(length(x),1);
elseif N==1
P=x;
else
L2=ones(length(x),1);
L1=x;

for k=2:N
P=(2*k-1)/k*x.*L1-(k-1)/k*L2;

L2=L1;
L1=P;
end

end
P=sqrt((2*N+1)/2)*P;


end