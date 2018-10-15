function [VGa,VGaX]=GenerateVanderMondeMatricesLG(V,Vx,N,BC,varargin)

if BC=='D'
ScaleVec=1./sqrt(4*[0:N-2]+6);
ScaleMat=repmat(ScaleVec,N+1,1);
VGa=ScaleMat.*(V(1:N+1,1:N-1)-V(1:N+1,3:(N+1)));
VGaX=ScaleMat.*(Vx(1:N+1,1:N-1)-Vx(1:N+1,3:(N+1)));
elseif BC=='N'
K=0:N-2;
bk=-K.*(K+1)./((K+2).*(K+3));
ScaleVec=1./sqrt(-bk.*(4*[0:N-2]+6));
ScaleMat=repmat(ScaleVec,N+1,1);
ScaleMat(:,1)=1/2;
BK=repmat(bk,N+1,1);
VGa=ScaleMat.*(V(1:N+1,1:N-1)+BK.*V(1:N+1,3:(N+1)));
VGaX=ScaleMat.*(Vx(1:N+1,1:N-1)+BK.*Vx(1:N+1,3:(N+1)));
elseif BC=='R'
K=0:N-2;
Cr=varargin{1};
ak=zeros(1,N-1);
bk=((-K.^2-K)-2*Cr)./((K.^2+5*K+6)+2*Cr);
ScaleVec=1./sqrt(-bk.*(4*K+6));
ScaleMat=repmat(ScaleVec,N+1,1);
AK=repmat(ak,N+1,1);
BK=repmat(bk,N+1,1);
VGa=ScaleMat.*(V(1:N+1,1:N-1)+AK.*V(1:N+1,2:N)+BK.*V(1:N+1,3:(N+1)));
VGaX=ScaleMat.*(Vx(1:N+1,1:N-1)+AK.*Vx(1:N+1,2:N)+BK.*Vx(1:N+1,3:(N+1)));
end