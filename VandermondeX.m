function Vx=VandermondeX(xj,N)

for k=0:N
[P,De]=JacobiPLegendreAndDerivative(xj,k);
Vx(:,k+1)=De;
end