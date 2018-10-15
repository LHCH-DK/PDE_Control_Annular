clear all
close all
clc

run Setup_main.m

for i=1:length(Nvec)
N=Nvec(i);
% Fourier-collocation points
j=0:2*M-1;
thetaj=pi*j/M;
x=-1:0.01:1;
% LGL points
[xj,wj] = JacobiGL(alpha,beta,N);
cc=(b+a)/(b-a);
rj=(b-a)/2*(xj+cc);

%Mesh
[Rj,Thetaj]=meshgrid(rj,thetaj);

%Legendre Vandermonde matrix
[P,V] = JacobiVandermonde(xj,alpha,beta,N);
Vx=VandermondeX(xj,N);

%Fourier Vandermonde matrix 
VFou=FourierVandermonde(thetaj,NN);
VFouX=FourierVandermondeX(thetaj,NN);

% Fourier-like basis - Pass "D" (Dirichlet BC) or "N" (Neumann BC)
Msmall=GenerateMassMatrixLG(N,'D');
[VGa,VGaX]=GenerateVanderMondeMatricesLG(V,Vx,N,'D');
[Q,D]=eig(Msmall);
Ds=sparse(D);
Vm=VGa*Q;
VmX=VGaX*Q;


%% Calculate rhs 
%Rhs
auxFull=(b-a)^2/4*(xj+cc);
AuxFull=repmat(auxFull',NN,1);

FFou=(1/NN*VFou).'*(AuxFull.*f(Rj,Thetaj))*(repmat(wj,1,N-1).*Vm);
ZFou=(1/NN*VFou).'*(AuxFull.*zd(Rj,Thetaj))*(repmat(wj,1,N-1).*Vm);
FFou=FFou';
ZFou=ZFou';


%% Initialize SNN scheme
iter=0;
DiffY=1;
DiffP=1;
YApproxFou=zeros(NN,N+1);
PApproxFou=zeros(NN,N+1);
FouLeY=zeros((N-1)*NN,1);
FouLeP=zeros((N-1)*NN,1);

%Setup iterative method
Vinv=V\eye(N+1);
VFouinv=((VFou').')\eye(NN);
WJ=repmat(wj,1,N+1);
MatV=((WJ.*V)'*Vm);
MatVX=((WJ.*V)'*VmX);
MatF=((1/NN*VFou)'*(VFou));
MatFX=((1/NN*VFou)'*(VFouX));
  
while ((DiffY >tolSSN & DiffP >tolSSN) & iter<=IterMax);
%Inner SNN Loop

% Make physical functions
YP=YApproxFou.*PApproxFou;
YY=YApproxFou.^2;

Hn=H(1/epsi*PApproxFou);
HPn=1/epsi*HP(1/epsi*PApproxFou);
    
Raux=@(r,theta)2/(b-a).*r;
Rauxinv=@(r,theta)(b-a)/2*1./r;
RR=Raux(Rj,Thetaj);
RRinv=Rauxinv(Rj,Thetaj);

%Variable rhs
auxFull=(b-a)^2/4*(xj+cc);
AuxFull=repmat(auxFull',NN,1);

YYYn=[];
YYYn=(1/NN*VFou).'*(AuxFull.*(YApproxFou.^3))*(repmat(wj,1,N-1).*Vm);
YYYn=YYYn';
YYYn=YYYn(:);

YYPn=[];
YYPn=(1/NN*VFou).'*(AuxFull.*(YY.*PApproxFou))*(repmat(wj,1,N-1).*Vm);
YYPn=YYPn';
YYPn=YYPn(:);

HN=(1/NN*VFou).'*(AuxFull.*Hn)*(repmat(wj,1,N-1).*Vm);
HN=HN';
HN=HN(:);

HPN=(1/NN*VFou).'*(AuxFull.*(HPn.*PApproxFou))*(repmat(wj,1,N-1).*Vm);
HPN=HPN';
HPN=HPN(:);

%Define rhs
rhs1=ZFou(:)+6*YYPn;
rhs2=2*YYYn+FFou(:)+HN-HPN;
RhsFullFou=[rhs1;rhs2];

%% Solve subsystems using preconditioned KSP method
%tolKSP=10^(-9);
MaxAq=max(cc+xj);
MinAq=min(cc+xj);
wcAq=1/2*(MaxAq+MinAq);
wcCq=wcAq;
% 
MaxBq=max(1./(cc+xj));
MinBq=min(1./(cc+xj));
wcBq=1/2*(MaxBq+MinBq);
% 
MaxYP=max(max(Raux(Rj,Thetaj).*YP));
MinYP=min(min(Raux(Rj,Thetaj).*YP));
MaxYY=max(max(Raux(Rj,Thetaj).*YY));
MinYY=min(min(Raux(Rj,Thetaj).*YY));
wcYP=1/2*(MaxYP+MinYP);
wcYY=1/2*(MaxYY+MinYY);
%
CA=wcAq;
CB=wcBq;
C0=BetaR*wcYY;
C1=1/(epsi)*BetaR*wcCq;
C2=BetaR*(wcYP+wcCq);

%Solve by GMRES
FouLESol=gmres(@(x)KKTfun_2D_Fou(x,N,NN,VFou,VFouX,Vm,VmX,VFouinv,Vinv,RR,RRinv,MatF,MatFX,MatV,MatVX,BetaR,YP,YY,HPn),RhsFullFou,[],tolKSP,100,@(x)Precond_2D_Fou(x,CA,CB,C0,C1,C2,D,N,NN),[],[FouLeY;FouLeP]);

%Extract solution
FouLeY=FouLESol(1:(N-1)*NN,1);
FouLeP=FouLESol((N-1)*NN+1:end,1);

%Reshape solution for summation procedure
FouLEMatY=reshape(FouLeY,(N-1),NN);
FouLEMatP=reshape(FouLeP,(N-1),NN);

%
YApproxFou_old=YApproxFou;
PApproxFou_old=PApproxFou;
%Represent solution by summation procedure
YApproxFou=real((VFou').'*(FouLEMatY)'*Vm');
PApproxFou=real((VFou').'*(FouLEMatP)'*Vm');
%Compute control
UApproxFou=H(1/epsi*PApproxFou);

%Count
iter=iter+1;
DiffY=max(max(YApproxFou-YApproxFou_old));
DiffP=max(max(PApproxFou-PApproxFou_old));
end
end

%% Plot
jj=0:2*M;
thetajj=pi*jj/M;
[Rj,Thetajj]=meshgrid(rj,thetajj);
Xj=Rj.*cos(Thetajj);
Yj=Rj.*sin(Thetajj);

%Show Optimal control control
figure(1)
surf(Xj,Yj,[UApproxFou;UApproxFou(1,:)])
colorbar;
axis square
view(2)

%Compare state to desired state
figure(2);
h1=subplot(1,2,1);
surf(Xj,Yj,[YApproxFou;YApproxFou(1,:)])
colorbar;
set(h1,'zlim',[min(min([YApproxFou;YApproxFou(1,:)])),max(max([YApproxFou;YApproxFou(1,:)]))]);
axis square
view(2)

h2=subplot(1,2,2);
surf(Xj,Yj,zd(Rj,Thetajj))
colorbar;
set(h2,'zlim',[min(min(zd(Rj,Thetajj))),max(max(zd(Rj,Thetajj)))]);
axis square
view(2)
