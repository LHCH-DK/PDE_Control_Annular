function VFouX=FourierVandermondeX(xj,N)

for k=-N/2:N/2-1
VFouX(:,k+N/2+1)=k*exp(1i*k*xj);
end