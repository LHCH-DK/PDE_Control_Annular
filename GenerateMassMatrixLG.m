function M=GenerateMassMatrixLG(N,BC,varargin)

K=0:N-2;

if BC=='D'
ak=zeros(1,N-1);
bk=-ones(1,N-1);
ScaleVec=1./sqrt(-bk.*(4*K+6));
elseif BC=='N'
ak=zeros(1,N-1);
bk=-K.*(K+1)./((K+2).*(K+3));
ScaleVec=1./sqrt(-bk.*(4*K+6));
ScaleVec(1,1)=1/2;
elseif BC=='R'
Cr=varargin{1};
ak=zeros(1,N-1);
bk=((-K.^2-K)-2*Cr)./((K.^2+5*K+6)+2*Cr);
ScaleVec=1./sqrt(-bk.*(4*K+6));
end

M1=ScaleVec.^2.*(2./(2*K+1)+ak.^2.*2./(2*K+3)+bk.^2.*2./(2*K+5));
M2=2*(ScaleVec(1:end-1).*ScaleVec(2:end)).*(ak(1:end-1).*2./(2*K(1:end-1)+3)+ak(2:end).*bk(1:end-1).*2./(2*K(1:end-1)+5));
M3=2*(ScaleVec(1:end-2).*ScaleVec(3:end)).*bk(1:end-2).*2./(2*K(1:end-2)+5);

Mtemp=diag(M1,0)+diag(M2,1)+diag(M3,2);
M=1/2*(Mtemp+Mtemp');
