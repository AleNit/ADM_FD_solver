
% solve steady advection-diffusion-migration equation in a cylindrical frame.
% The quantity q is transported by an incompressible flow with velocity u:
%       cD*Lapl(q) + grad(q)*cB - cA*q = rhs
% by centered finite differences on a cell-centered, non-uniform grid; 
% b.c.s are prescribed by a ghost cell method with generalized {\alpha,\beta,\gamma} coefficients; 
% Boundary conditions: r --> (inner-Dirichelet/Neumann), theta --> (periodic),  z --> (periodic)

% A. Nitti, Polytechnic University of Bari (2024)


clc
clear 
close all
clearAllMemoizedCaches
kf=1;


%% parameters for the convergence study
Lr=0.5;                                                         % domain length in the radial direction
Lz=2;                                                           % domain length in the axial direction
nrc=[4,8,16]+1;

% analytical solution and coefficient functions
qa      = @(r,t,z) r.^2.*cos(2.*t).*sin(pi.*z);
coeffA  = @(r,t,z)(sin(t.*2.0).*3.0)./r;
coeffBr = @(r,t,z)-sin(t.*2.0)+cos(r.*pi.*3.0).*sin(t);
coeffBt = @(r,t,z)cos(t.*2.0).*-2.0+cos(z.*pi);
coeffBz = @(r,t,z)-(z.*sin(t).*(cos(r.*pi.*3.0)-r.*pi.*sin(r.*pi.*3.0).*3.0))./r;
coeffD  = 0.8;
rhs     = @(r,t,z)r.*cos(t.*2.0).*sin(z.*pi).*(sin(t.*2.0)- ...
          cos(r.*pi.*3.0).*sin(t)).*-2.0+r.*sin(t.*2.0).*sin(z.*pi).* ...
          (cos(t.*2.0).*2.0-cos(z.*pi)).*2.0-r.^2.*pi.^2.*cos(t.*2.0).* ...
          sin(z.*pi).*(4.0./5.0)+r.*cos(t.*2.0).*sin(t.*2.0).*sin(z.*pi).* ...
          3.0-r.*z.*pi.*cos(t.*2.0).*cos(z.*pi).*sin(t).*(cos(r.*pi.*3.0)-r.* ...
          pi.*sin(r.*pi.*3.0).*3.0);




for C=1:length(nrc)

disp(['------------------ GRID ',num2str(C)])    

%% input parameters
nr=nrc(C);                                                          % number of nodes in the radial direction
nt=ceil(nr*4)+1;                                                          % number of nodes in the tangential direction
nz=ceil(Lz/Lr*nr)+1;                                                          % number of nodes in the axial direction
refr='lin';                                                  % nodes distribution int the radial direction

% boundary condition at r=Lr
valr1= @(t,z) cat(3,ones(nt-1,nz-1),zeros(nt-1,nz-1),qa(Lr,t,z));           % Dirichelet



%% pre-processing operations
% check if number of cell centers is even
if ( rem(nt-1,2)~=0 )
    error('... nt-1 must be even to fulfil the axis mirroring condition')
end

% create and plot grid
gr=getgrid(Lr,nr,refr,false(1));
gt=getgrid(2*pi,nt,'lin',false(1));
gz=getgrid(Lz,nz,'lin',false(1));
ndof=(nr-1)*(nt-1)*(nz-1);

% kf=plotgrid(gr.xn,gt.xn,gz.xn,Lr,Lz,kf);
disp(['ndof = ',num2str(ndof)])

% assemble coefficient matrix
tic
[Tp,Rp,Zp]=meshgrid(gt.xp,gr.xp,gz.xp);
RHS=rhs(Rp,Tp,Zp);
cA=coeffA(Rp,Tp,Zp);
cBr=coeffBr(Rp,Tp,Zp);
cBt=coeffBt(Rp,Tp,Zp);
cBz=coeffBz(Rp,Tp,Zp);
[VA,I,J,h,b]=getCoeffMat3Dc(gr,gt,gz,RHS,ndof,cA,cBr,cBt,cBz,coeffD);

% assign boundary conditions
[Zc,Tc]=meshgrid(gz.xp,gt.xp);
bcr1=valr1(Tc,Zc);
[A,b]=bcs3Dc(gr,gt,gz,VA,I,J,h,b,RHS,bcr1,cA,cBr,cBt,cBz,coeffD);
asst=toc;
disp(['matrix assembly time: ',num2str(asst)])



%% solve problem (BCGSTAB iterative solver)
maxit=10000;
tol=1.0e-6;
tic
qs = bicgstab(A,b,tol,maxit);
solt=toc;

q=reshape(qs,[nr-1,nt-1,nz-1]);
disp(['system solution time: ',num2str(solt)])



%% compute error norm and plot results
[Tp,Rp]=meshgrid(gt.xp,gr.xp);
k1=ceil( (nz-1)*1/3 );
q1=q(:,:,k1);
qa1=qa(Rp,Tp,gz.xp(k1));
k2=ceil( (nz-1)*2/3 );
q2=q(:,:,k2);
qa2=qa(Rp,Tp,gz.xp(k2));

figure(kf); kf=kf+1;
set(gcf,'Position',[100,100,1000,450])
tiledlayout(1,2);

nexttile
[X,Y,V1]=pol2cart(Tp,Rp,qa1);
surf(X,Y,V1)
shading interp
hold on
surf(X,Y,q1,'FaceColor','none')
axis equal
xlabel('x','interpreter','latex');    ylabel('y','interpreter','latex');    zlabel('q','interpreter','latex');
legend('analytical','numerical')
tit=strcat('r-$\theta$ plane, k=',num2str(k1),'/',num2str(nz-1));
title(tit,'interpreter','latex')

nexttile
[X,Y,V2]=pol2cart(Tp,Rp,qa2);
surf(X,Y,V2)
shading interp
hold on
surf(X,Y,q2,'FaceColor','none')
axis equal
xlabel('x','interpreter','latex');    ylabel('y','interpreter','latex');    zlabel('q','interpreter','latex');
legend('analytical','numerical')
tit=strcat('r-$\theta$ plane, k=',num2str(k2),'/',num2str(nz-1));
title(tit,'interpreter','latex')


% compute rmse
err1=( qa1-q1 )./(qa1+1);
disp(['relative rmse section 1: ',num2str(rms(err1,'all'))])

errC(C)=rms(err1,'all');


end



%% convergence plot
grids = Lr./nrc;
p =  abs((log(errC(2))-log(errC(1)))/(log(grids(2))-log(grids(1))));

figure(kf); kf=kf+1;
loglog(grids,errC,'-s','color','b');
hold on
y_slope = 0.1.*grids.^2;
loglog(grids,y_slope,'--','color','r');
xlabel('$\Delta r$','interpreter','latex','fontsize',14)
ylabel('RMSE','interpreter','latex','fontsize',14);
legend(strcat('p=',num2str(p)),'slope=2', 'location','northwest',...
    'interpreter','latex','fontsize',12)
title('r-$\theta$ plane','interpreter','latex','fontsize',14);


