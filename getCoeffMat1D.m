
% assemble coefficient matrix for the 1D problem

function A=getCoeffMat1D(gr,ndof,coeffA,coeffBr,cD)

nr=length(gr.xn);


%% full matrix
A=zeros(ndof,ndof);
for i=2:nr-2

    ip=i+1;
    im=i-1;

    cA=coeffA(gr.xp(i));
    cB=coeffBr(gr.xp(i));
    d=1/( gr.dxc(i)*gr.dxc(ip)*( gr.dxc(i) + gr.dxc(ip) ) );

    A(i,i)  =   cA + cB*(gr.dxc(ip)^2 - gr.dxc(i)^2)*d ...
                -cD/(gr.dxn(i)*gr.dxc(ip)) - cD/(gr.dxn(i)*gr.dxc(i));

    A(i,ip) =   cB*gr.dxc(i)^2*d + cD/(gr.dxn(i)*gr.dxc(ip));

    A(i,im) =   -cB*gr.dxc(ip)^2*d + cD/(gr.dxn(i)*gr.dxc(i));

end


end
