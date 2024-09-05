
function gr=getgrid(Lx,nx,rref,pbc)

% with periodic boundary conditions impose a uniform grid spacing
if (pbc)
    rref='lin';
end

% node grid with ghost nodes
switch rref
    case 'lin'
        % xn=linspace(0,Lx,nx);
        dx=Lx/(nx-1);
        xni=-dx:dx:Lx+dx;
    case 'tanh-i'               % refinement on the initial side
        c=2;
        xni=tanh(c.*(-1:nx)./(nx-1));
        xni=xni.*Lx/xni(end-1);
        xni=Lx-xni; 
        xni=flip(xni);
    case 'tanh-e'               % refinement on the final side
        c=2;   
        xni=tanh(c.*(-1:nx)./(nx-1));
        xni=xni.*Lx/xni(end-1);           
    otherwise
        error('... grid law not implemented')
end

% node grid
xn=xni(2:end-1);
dxn=xn(2:end)-xn(1:end-1);

% cell center grid with ghost cells
xc=(xni(2:end)+xni(1:end-1)).*0.5;
dxc=xc(2:end)-xc(1:end-1);

% cell center grid
xp=xc(2:end-1);
dxp=xp(2:end)-xp(1:end-1);

% assign everything to the gr object
gr.xn=xn;
gr.xc=xc;
gr.xp=xp;
gr.dxn=dxn;
gr.dxc=dxc;
gr.dxp=dxp;

end
