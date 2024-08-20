
% solve the steady advection-diffusion-migration equation in a 2D polar frame.
% The quantity q is transported by an incompressible flow with velocity u:
%       cD*Lapl(q) + grad(q)*cB - cA*q = rhs
% Use centered finite differences on a cell-centered, non-uniform grid; 
% b.c.s are prescribed by a ghost cell method with generalized {\alpha,\beta,\gamma} coefficients; 
% Boundary conditions: r --> (inner-Dirichelet/Neumann), theta --> (periodic)

% A. Nitti, Polytechnic University of Bari (2024)


% clc
clear 
close all
clearAllMemoizedCaches
kf=1;


%% input parameters
Lr=0.5;
nr=21;                                                          % number of nodes
nt=33;                                                          % number of nodes
refr='tanh-e';                                                  % nodes distribution 

% analytical solution and coefficient functions
qa      = @(r,t) r.^2.*cos(2.*t);
coeffA  = @(r,t) (sin(t.*2.0).*3.0)./r;
coeffBr = @(r,t) -sin(t.*2.0)+cos(r.*pi.*3.0).*sin(t);
coeffBt = @(r,t) cos(t).^2.*-4.0+cos(t).*(cos(r.*pi.*3.0)- ...
          r.*pi.*sin(r.*pi.*3.0).*3.0)+2.0;
coeffD  = 0.1;
rhs     = @(r,t) r.*sin(t).*(cos(r.*pi.*3.0)+cos(t).*5.0- ...
          cos(t).^3.*1.0e+1-r.*pi.*sin(r.*pi.*3.0).*cos(t).^2.*6.0).*-2.0;

% boundary conditions at r=Lr
valr1= @(t) cat(3,ones(nt-1,1),zeros(nt-1,1),qa(Lr,t));         % Dirichelet
% valr1= @(t) cat(3,zeros(nt-1,1),ones(nt-1,1),2.*Lr.*cos(2.*t));         % Neumann



%% pre-processing operations
% check if number of cell centers is even
if ( rem(nt-1,2)~=0 )
    error('... nt-1 must be even to fulfil the axis mirroring condition')
end

% create and plot grid
gr=getgrid(Lr,nr,refr,false(1));
gt=getgrid(2*pi,nt,'lin',false(1));
ndof=(nr-1)*(nt-1);

% assemble coefficient matrix
[Tn,Rn]=meshgrid(gt.xp,gr.xp);
RHS=rhs(Rn,Tn);
cA=coeffA(Rn,Tn);
cBr=coeffBr(Rn,Tn);
cBt=coeffBt(Rn,Tn);
[A,b]=getCoeffMat2D(gr,gt,RHS,ndof,cA,cBr,cBt,coeffD);

% assign boundary conditions
bcr1=valr1(gt.xp');
[A,b]=bcs2D(gr,gt,A,b,RHS,bcr1,cA,cBr,cBt,coeffD);


%% solve problem
% compute matrix condition number
dA=decomposition(A);
tf = isIllConditioned(dA);
cnum=rcond(dA);
disp(['rciprocal of condition number of coeff. matrix: ',num2str(cnum)])
if (tf)
    error('ill-conditioned coefficient matrix:')    
end

% solve problem
us=A\b;

q=reshape(us,[nr-1,nt-1]);


%% compute error norm and plot results
figure(kf); kf=kf+1;
uan=qa(Rn,Tn);
[X,Y,Va]=pol2cart(Tn,Rn,uan);
surf(X,Y,Va)
shading interp
hold on
surf(X,Y,q,'FaceColor','none')
axis equal
xlabel('x');    ylabel('y');    zlabel('q');
legend('analytical','numerical')


% compute rmse
err1=( uan-q )./(uan+1);
disp(['relative rmse: ',num2str(rms(err1,'all'))])




