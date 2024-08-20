
function [A,b]=bcs3Dc(gr,gt,gz,VA,I,J,h,b,rhs,vale,coeffA,coeffBr,coeffBt,coeffBz,cD)

nr=length(gr.xn);
nt=length(gt.xn);
nz=length(gz.xn);
ndof=(nr-1)*(nt-1)*(nz-1);


% grid spacing in periodic directions
dt=gt.xn(2)-gt.xn(1);
dz=gz.xn(2)-gz.xn(1);


% indices for periodicity
indjp=2:nt;     indjp(end)=1;
indjm=0:nt-2;   indjm(1)=nt-1;

indkp=2:nz;     indkp(end)=1;
indkm=0:nz-2;   indkm(1)=nz-1;


% initial boundary - mirror value
i=1;
ip=i+1;
d=1/( gr.dxc(i)*gr.dxc(ip)*( gr.dxc(i) + gr.dxc(ip) ) );


for k=1:nz-1
    kp=indkp(k);
    km=indkm(k);
    
    for j=1:nt-1
        jp=indjp(j);
        jm=indjm(j);

        cA=coeffA(i,j,k);
        cBr=coeffBr(i,j,k);
        cBt=coeffBt(i,j,k);
        cBz=coeffBz(i,j,k);

        c   = (nr-1)*(nt-1)*(k-1) + (nr-1)*(j-1) + i;
        jsym=mod( j+(nt-1)/2-1 , nt-1 );
        cim = (nr-1)*(nt-1)*(k-1) + (nr-1)*jsym + i;
        cip = (nr-1)*(nt-1)*(k-1) + (nr-1)*(j-1)  + i + 1;
        cjp = (nr-1)*(nt-1)*(k-1) + (nr-1)*(jp-1) + i;
        cjm = (nr-1)*(nt-1)*(k-1) + (nr-1)*(jm-1) + i;
        ckp = (nr-1)*(nt-1)*(kp-1) + (nr-1)*(j-1)  + i;
        ckm = (nr-1)*(nt-1)*(km-1) + (nr-1)*(j-1)  + i;        
        
        I(h)    =   c;
        J(h)    =   c;
        VA(h)   =   cA + cBr*(gr.dxc(ip)^2 - gr.dxc(i)^2)*d + ...
                    - cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip)) + ...
                    - cD*gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i)) + ...
                    - cD*2/(gr.xc(ip)^2*dt^2) + ...
                    - cD*2/dz^2;

        I(h+1)  =   c;
        J(h+1)  =   cip;        
        VA(h+1) =   cBr*gr.dxc(i)^2*d + ...
                    + cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));
        
        I(h+2)  =   c;
        J(h+2)  =   cim;        
        VA(h+2) =   -cBr*gr.dxc(ip)^2*d + ...
                    + cD*gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i));

        I(h+3)  =   c;
        J(h+3)  =   cjp;        
        VA(h+3) =   cBt/(2*gr.xp(i)*dt) + cD/(gr.xc(ip)^2*dt^2);
        
        I(h+4)  =   c;
        J(h+4)  =   cjm;        
        VA(h+4) =   -cBt/(2*gr.xp(i)*dt) + cD/(gr.xc(ip)^2*dt^2);

        I(h+5)  =   c;
        J(h+5)  =   ckp;        
        VA(h+5) =   cBz/(2*dz) + cD/dz^2;

        I(h+6)  =   c;
        J(h+6)  =   ckm;        
        VA(h+6) =   -cBz/(2*dz) + cD/dz^2;

        b(c)    =   rhs(i,j,k);   

        h=h+7;
    
    end
end



% final boundary
i=nr-1;
ip=i+1;
d=1/( gr.dxc(i)*gr.dxc(ip)*( gr.dxc(i) + gr.dxc(ip) ) );

for k=1:nz-1
    kp=indkp(k);
    km=indkm(k);

    for j=1:nt-1
        jp=indjp(j);
        jm=indjm(j);
    
        cA=coeffA(i,j,k);
        cBr=coeffBr(i,j,k);
        cBt=coeffBt(i,j,k);
        cBz=coeffBz(i,j,k);

        c   = (nr-1)*(nt-1)*(k-1) + (nr-1)*(j-1) + i;        
        cim = (nr-1)*(nt-1)*(k-1) + (nr-1)*(j-1) + i - 1;
        cjp = (nr-1)*(nt-1)*(k-1) + (nr-1)*(jp-1) + i;
        cjm = (nr-1)*(nt-1)*(k-1) + (nr-1)*(jm-1) + i;
        ckp = (nr-1)*(nt-1)*(kp-1) + (nr-1)*(j-1)  + i;
        ckm = (nr-1)*(nt-1)*(km-1) + (nr-1)*(j-1)  + i;

        c3=(-2*vale(j,k,2)+vale(j,k,1)*gr.dxc(ip))/(-2*vale(j,k,2)-vale(j,k,1)*gr.dxc(ip));
        c4=-2*vale(j,k,3)*gr.dxc(ip)/(-2*vale(j,k,2)-vale(j,k,1)*gr.dxc(ip));

        I(h)    =   c;
        J(h)    =   c;    
        VA(h)   =   cA + cBr*(gr.dxc(ip)^2 - gr.dxc(i)^2)*d + ...
                      - cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip)) + ...
                      - cD*gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i)) + ...
                      - cD*2/(gr.xc(ip)^2*dt^2) + ...
                      - cD*2/dz^2 + ...
                      + c3*cBr*d*gr.dxc(i)^2 + c3*cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));
    
        I(h+1)  =   c;
        J(h+1)  =   cim;        
        VA(h+1) =  -cBr*d*gr.dxc(ip)^2 + ...
                      + cD*gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i));
    
        I(h+2)  =   c;
        J(h+2)  =   cjp;        
        VA(h+2) =  cBt/(2*gr.xp(i)*dt) + cD/(gr.xc(ip)^2*dt^2);
        
        I(h+3)  =   c;
        J(h+3)  =   cjm;        
        VA(h+3) =  -cBt/(2*gr.xp(i)*dt) + cD/(gr.xc(ip)^2*dt^2); 

        I(h+4)  =   c;
        J(h+4)  =   ckp;        
        VA(h+4) =   cBz/(2*dz) + cD/dz^2;

        I(h+5)  =   c;
        J(h+5)  =   ckm;        
        VA(h+5) =   -cBz/(2*dz) + cD/dz^2;        
    
        b(c)    =  rhs(i,j,k) - c4*cBr*d*gr.dxc(i)^2 + ... 
                  -c4*cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));

        h=h+6;
    
    end
end


A=sparse(I,J,VA,ndof,ndof);

end
