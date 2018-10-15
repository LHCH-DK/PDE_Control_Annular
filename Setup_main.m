% Use Legendre polynomials
alpha=0; beta=0;
%Number of Legendre modes
Nvec=[51];
%Number of Fourier modes
NN=50;
%Number of trig modes
M=NN/2;
%Regularization paramter
epsi=10^(-2);
%Annulus
b=60;
a=30;
%Scaling parameter, kappa - Do not change.
BetaR=(b-a)^2/4;
%
%Box-constraints
ua=-40;
ub=40;
H=@(x)max(ua,min(x,ub));
HP=@(x)(x<=ub & x>=ua).*1;
%Desired state
zd=@(r,theta)(r>=40 & r<=60 & theta>=0 & theta<=pi/2 | r>=40 & r<=60 & theta>=pi & theta<=3/2*pi ).*4+0.*r.*theta;
% Rhs source-term
f=@(r,theta)r.*theta.*0;
%Max Number of SSN iterations
IterMax=10;
%Tolerance SSN solver
tolSSN=10^(-4);
%Tolerance KSP solver
tolKSP=10^(-9);