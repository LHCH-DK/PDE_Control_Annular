function VFou=FourierVandermonde(xj,N)

for k=-N/2:N/2-1
VFou(:,k+N/2+1)=exp(1i*k*xj);
end