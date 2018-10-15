function [P,De]=JacobiPLegendreAndDerivative(x,N)

if N==0
P=ones(length(x),1);
De=zeros(length(x),1);
elseif N==1
P=x;
De=ones(length(x),1);
else
L2=ones(length(x),1);
L1=x;
LPrime2=zeros(length(x),1);
LPrime1=ones(length(x),1);
for k=2:N
P=(2*k-1)/k*x.*L1-(k-1)/k*L2;
De=LPrime2+(2*k-1).*L1;

L2=L1;
L1=P;

LPrime2=LPrime1;
LPrime1=De;
end

end
P=P;
De=De;
end