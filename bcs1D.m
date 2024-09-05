
% assign boundary conditions in the form:
%       \alpha*u + \beta*dudx = \gamma.

function [A,b]=bcs1D(A,b,gr,rhs,vali,vale,coeffA,coeffBr,cD)

nx=length(gr.xn);

% initial boundary
i=1;
ip=i+1;
cA=coeffA(gr.xp(i));
cB=coeffBr(gr.xp(i));
d=1/( gr.dxc(i)*gr.dxc(ip)*( gr.dxc(i) + gr.dxc(ip) ) );

den=(gr.dxc(1)-gr.dxn(1)/2)*gr.dxn(1)/2*gr.dxc(1);
dp=(gr.dxc(1)-gr.dxn(1)/2)^2/den;
dm=(gr.dxn(1)/2)^2/den;
dd=( (gr.dxn(1)/2)^2 - (gr.dxc(1)-gr.dxn(1)/2)^2 )/den;
c1=( vali(2)*dp + (2*gr.dxc(1)-gr.dxn(1))/2/gr.dxc(1)*( vali(1)+vali(2)*dd ) )/ ...
   ( vali(2)*dm + ( (2*gr.dxc(1)-gr.dxn(1))/2/gr.dxc(1)-1 )*( vali(1)+vali(2)*dd ) );
c2= - vali(3) / ( vali(2)*dm + ( (2*gr.dxc(1)-gr.dxn(1))/2/gr.dxc(1)-1 )*( vali(1)+vali(2)*dd ) );

A(i, i)     =  cA + cB*d*(gr.dxc(ip)^2-gr.dxc(i)^2) ...
               - cD/(gr.dxn(i)*gr.dxc(ip)) - cD/(gr.dxn(i)*gr.dxc(i)) ...
               - c1*cB*d*gr.dxc(ip)^2 + c1*cD/(gr.dxn(i)*gr.dxc(i));
A(i, ip)    =  cB*d*gr.dxc(i)^2 + cD/(gr.dxn(i)*gr.dxc(ip));
b(i)        =  rhs(i) + c2*cB*d*gr.dxc(ip)^2 - c2*cD/(gr.dxn(i)*gr.dxc(i));



% final boundary
i=nx-1;
ip=i+1;
cA=coeffA(gr.xp(i));
cB=coeffBr(gr.xp(i));
d=1/( gr.dxc(i)*gr.dxc(ip)*( gr.dxc(i) + gr.dxc(ip) ) );

den=(gr.dxc(nx)-gr.dxn(nx-1)/2)*gr.dxn(nx-1)/2*gr.dxc(nx);
dp=(gr.dxn(nx-1)/2)^2/den;
dm=(gr.dxc(nx)-gr.dxn(nx-1)/2)^2/den;
dd=( (gr.dxc(nx)-gr.dxn(nx-1)/2)^2 - (gr.dxn(nx-1)/2)^2 )/den;
c3=( vale(2)*dm + ( gr.dxn(nx-1)/2/gr.dxc(nx)-1 )*( vale(1)+vale(2)*dd ) )/ ...
( vale(2)*dp + gr.dxn(nx-1)/2/gr.dxc(nx)*( vale(1)+vale(2)*dd ) );
c4= vale(3) / ( vale(2)*dp + gr.dxn(nx-1)/2/gr.dxc(nx)*( vale(1)+vale(2)*dd ) );

A(i, i)     =  cA + cB*d*(gr.dxc(ip)^2-gr.dxc(i)^2) ...
               - cD/(gr.dxn(i)*gr.dxc(ip)) - cD/(gr.dxn(i)*gr.dxc(i)) ...
               + c3*cB*d*gr.dxc(i)^2 + c3*cD/(gr.dxn(i)*gr.dxc(ip));
A(i, i-1)   =  -cB*d*gr.dxc(ip)^2 + cD/(gr.dxn(i)*gr.dxc(i));
b(i)        =  rhs(i) - c4*cB*d*gr.dxc(i)^2 - c4*cD/(gr.dxn(i)*gr.dxc(ip));        


end
