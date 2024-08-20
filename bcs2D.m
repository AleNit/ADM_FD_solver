
function [A,b]=bcs2D(gr,gt,A,b,rhs,vale,coeffA,coeffBr,coeffBt,cD)

nr=length(gr.xn);
nt=length(gt.xn);

dt=gt.xn(2)-gt.xn(1);

indjp=2:nt;     indjp(end)=1;
indjm=0:nt-2;   indjm(1)=nt-1;


% initial boundary -------------------- > mirror value
i=1;
ip=i+1;
d=1/( gr.dxc(i)*gr.dxc(ip)*( gr.dxc(i) + gr.dxc(ip) ) );

for j=1:nt-1
    jp=indjp(j);
    jm=indjm(j);

    cA=coeffA(i,j);
    cBr=coeffBr(i,j);
    cBt=coeffBt(i,j);

    c   = (nr-1)*(j-1) + i;
    jsym=mod( j+(nt-1)/2-1 , nt-1 );
    cim = (nr-1)*jsym + i;
    cip = (nr-1)*(j-1)  + i + 1;
    cjp = (nr-1)*(jp-1) + i;
    cjm = (nr-1)*(jm-1) + i;

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



% final boundary
i=nr-1;
ip=i+1;
d=1/( gr.dxc(i)*gr.dxc(ip)*( gr.dxc(i) + gr.dxc(ip) ) );

for j=1:nt-1
    jp=indjp(j);
    jm=indjm(j);

    cA=coeffA(i,j);
    cBr=coeffBr(i,j);
    cBt=coeffBt(i,j);

    c   = (nr-1)*(j-1)  + i;
    cim = (nr-1)*(j-1)  + i - 1;
    cjp = (nr-1)*(jp-1) + i;
    cjm = (nr-1)*(jm-1) + i;

    c3=(-2*vale(j,1,2)+vale(j,1,1)*gr.dxc(ip))/(-2*vale(j,1,2)-vale(j,1,1)*gr.dxc(ip));
    c4=-2*vale(j,1,3)*gr.dxc(ip)/(-2*vale(j,1,2)-vale(j,1,1)*gr.dxc(ip));

    A(c,c)     =  cA + cBr*(gr.dxc(ip)^2 - gr.dxc(i)^2)*d + ...
                  - cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip)) + ...
                  - cD*gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i)) + ...
                  - cD*2/(gr.xc(ip)^2*dt^2) + ...
                  + c3*cBr*d*gr.dxc(i)^2 + c3*cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));

    A(c,cim)   =  -cBr*d*gr.dxc(ip)^2 + ...
                  + cD*gr.xn(i)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(i));

    A(c,cjp)   =  cBt/(2*gr.xp(i)*dt) + cD/(gr.xc(ip)^2*dt^2);
    
    A(c,cjm)   =  -cBt/(2*gr.xp(i)*dt) + cD/(gr.xc(ip)^2*dt^2);    

    b(c)       =  rhs(i,j) - c4*cBr*d*gr.dxc(i)^2 + ... 
                  -c4*cD*gr.xn(ip)/(gr.xc(ip)*gr.dxn(i)*gr.dxc(ip));

end
