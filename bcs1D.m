
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

c1=(2*vali(2)+vali(1)*gr.dxc(i))/(2*vali(2)-vali(1)*gr.dxc(i));
c2=-2*vali(3)*gr.dxc(i)/(2*vali(2)-vali(1)*gr.dxc(i));

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

c3=(-2*vale(2)+vale(1)*gr.dxc(i+1))/(-2*vale(2)-vale(1)*gr.dxc(i+1));
c4=-2*vale(3)*gr.dxc(i+1)/(-2*vale(2)-vale(1)*gr.dxc(i+1));

A(i, i)     =  cA + cB*d*(gr.dxc(ip)^2-gr.dxc(i)^2) ...
               - cD/(gr.dxn(i)*gr.dxc(ip)) - cD/(gr.dxn(i)*gr.dxc(i)) ...
               + c3*cB*d*gr.dxc(i)^2 + c3*cD/(gr.dxn(i)*gr.dxc(ip));
A(i, i-1)   =  -cB*d*gr.dxc(ip)^2 + cD/(gr.dxn(i)*gr.dxc(i));
b(i)        =  rhs(i) - c4*cB*d*gr.dxc(i)^2 - c4*cD/(gr.dxn(i)*gr.dxc(ip));        


end
