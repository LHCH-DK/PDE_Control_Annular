function [P,V] = JacobiVandermonde(x,a,b,n)

a1 = @(j) 2*(j+a)*(j+b)/((2*j+a+b+1)*(2*j+a+b));
a2 = @(j) (a^2-b^2)/((2*j+a+b+2)*(2*j+a+b));
a3 = @(j) 2*(j+1)*(j+a+b+1)/((2*j+a+b+2)*(2*j+a+b+1));

p0 = 0*x+1;
p1 = (a-b+(a+b+2)*x)/2;

V = zeros(n+1,length(x));
V(1,:) = p0;

if n >= 2
    P = ((a2(1)+x).*p1-a1(1)*1)/a3(1);
    Pm = p1;
    
    V(2,:) = p1;
    V(3,:) = P;
elseif n > 0
    P = p1;
    V(2,:) = p1;
else
    P = 0*x+1;
end
for k = 2:n-1
    tmp = P;
    P = ((a2(k)+x).*P-a1(k)*Pm)/a3(k);
    Pm = tmp;
    
    V(k+2,:) = P;
end
%gamma_n = 2^(a+b+1)*gamma(n+a+1)*gamma(n+b+1)/(factorial(n)*(2*n+a+b+1)*gamma(n+a+b+1));
%gamma_n = trapz(x,P.^2);
%P = 1/sqrt(gamma_n)*P;
V = V.';