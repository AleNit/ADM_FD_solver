
% assemble coefficient matrix in sparse format

function [A,b]=getCoeffMat2D(gr,gt,rhs,ndof,coeffA,coeffBr,coeffBt,cD)

nr=length(gr.xn);
nt=length(gt.xn);

dt=gt.xn(2)-gt.xn(1);

indjp=2:nt;     indjp(end)=1;
indjm=0:nt-2;   indjm(1)=nt-1;

A=zeros(ndof,ndof);
b=zeros(ndof,1);

for j=1:nt-1
    jp=indjp(j);
    jm=indjm(j);

    for i=2:nr-2

        ip=i+1;
          
        c   = (nr-1)*(j-1)  + i;
        cip = (nr-1)*(j-1)  + i + 1;
        cim = (nr-1)*(j-1)  + i - 1;
        cjp = (nr-1)*(jp-1) + i;
        cjm = (nr-1)*(jm-1) + i;
                
        cA=coeffA(i,j);
        cBr=coeffBr(i,j);
        cBt=coeffBt(i,j);
        d=1/( gr.dxc(i)*gr.dxc(ip)*( gr.dxc(i) + gr.dxc(ip) ) );
        
        A(c,c)  =   cA + cBr*(gr.dxc(ip)^2 - gr.dxc(i)^2)*d + ...
                    - cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip)) + ...
                    - cD*gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i)) + ...
                    - cD*2/(gr.xc(ip)^2*dt^2);
        
        A(c,cip)=   cBr*gr.dxc(i)^2*d + ...
                    + cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));
        
        A(c,cim)=   -cBr*gr.dxc(ip)^2*d + ...
                    + cD*gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i));

        A(c,cjp)=   cBt/(2*gr.xp(i)*dt) + cD/(gr.xc(ip)^2*dt^2);
        
        A(c,cjm)=   -cBt/(2*gr.xp(i)*dt) + cD/(gr.xc(ip)^2*dt^2);

        b(c)    =   rhs(i,j);

    end
end


end
