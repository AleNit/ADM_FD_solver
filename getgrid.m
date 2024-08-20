
function gr=getgrid(Lx,nx,rref,pbc)

% with periodic boundary conditions impose a uniform grid spacing
if (pbc)
    rref='lin';
end

% node grid
switch rref
    case 'lin'
        xn=linspace(0,Lx,nx);
    case 'tanh-i'               % refinement on the initial side
        c=2;
        xn=tanh(c.*(0:nx-1)./(nx-1));
        xn=xn.*Lx/xn(end);
        xn=Lx-xn; 
        xn=flip(xn);
    case 'tanh-e'               % refinement on the final side
        c=2;
        xn=tanh(c.*(0:nx-1)./(nx-1));
        xn=xn.*Lx/xn(end);           
    otherwise
        error('... grid law not implemented')
end

dxn=xn(2:end)-xn(1:end-1);

% cell center grid
xp=(xn(2:end)+xn(1:end-1)).*0.5;
dxp=xp(2:end)-xp(1:end-1);

% cell center grid with ghost cells
xc=[-xp(1),xp,2.*xn(end)-xp(end)];
dxc=xc(2:end)-xc(1:end-1);

% assign everything to the grd object
gr.xn=xn;
gr.xc=xc;
gr.xp=xp;
gr.dxn=dxn;
gr.dxc=dxc;
gr.dxp=dxp;

end
