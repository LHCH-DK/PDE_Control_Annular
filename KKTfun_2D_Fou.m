function y=KKTfun_2D_Fou(x,N,M,VFou,VFouX,Vm,VmX,VFouinv,Vinv,Raux,RauxInv,MatF,MatFX,MatV,MatVX,BetaR,YP,YY,HPn)
  %% Matrix-free matrix-vector products
  x1=x(1:(N-1)*M);
  x2=x((N-1)*M+1:end); 
  
  Xhat=reshape(x1,N-1,M);
  X=(VFou').'*(Xhat)'*Vm';
  W=VFouinv*(Raux.*YP.*X)*Vinv';
  temp=MatF*W*MatV;
  temp=temp';
  MYPx1=temp(:);

  W=VFouinv*(Raux.*YY.*X)*Vinv';
  temp=MatF*W*MatV;
  temp=temp';
  MYYx1=temp(:);

  W=VFouinv*(Raux.*X)*Vinv';
  temp=MatF*W*MatV;
  temp=temp';
  CQFoux1=temp(:);

  X=(VFou').'*(Xhat)'*VmX';
  W=VFouinv*(Raux.*X)*Vinv';
  temp=MatF*W*MatVX;
  temp=temp';
  AQFoux1=temp(:);

  X=(VFouX').'*(Xhat)'*Vm';
  W=VFouinv*(RauxInv.*X)*Vinv';
  temp=MatFX*W*MatV;
  temp=temp';
  BQFoux1=temp(:);


  Xhat2=reshape(x2,N-1,M);
  X2=(VFou').'*(Xhat2)'*Vm';
  W=VFouinv*(Raux.*HPn.*X2)*Vinv';
  temp=MatF*W*MatV;
  temp=temp';
  MHPnx2=temp(:);

  W=VFouinv*(Raux.*YY.*X2)*Vinv';
  temp=MatF*W*MatV;
  temp=temp';
  MYYx2=temp(:);
 
  X2=(VFou').'*(Xhat2)'*VmX';
  W=VFouinv*(Raux.*X2)*Vinv';
  temp=MatF*W*MatVX;
  temp=temp';
  AQFoux2=temp(:);

  X2=(VFouX').'*(Xhat2)'*Vm';
  W=VFouinv*(RauxInv.*X2)*Vinv';
  temp=MatFX*W*MatV;
  temp=temp';
  BQFoux2=temp(:);

 y=[BetaR*(6*MYPx1+CQFoux1)+(AQFoux2+BQFoux2+3*BetaR*MYYx2);(AQFoux1+BQFoux1+3*BetaR*MYYx1)-BetaR*MHPnx2];
end