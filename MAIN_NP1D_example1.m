
% solve a steady 1D advection-diffusion-migration equation in one dimenson:
%       cD*d^2(q)/dx^2 + cB*d(q)/dx - cA*q = rhs
% Use centered finite differences on a cell-centered, non-uniform grid; 
% b.c.s are prescribed by a ghost cell method with generalized {\alpha,\beta,\gamma} coefficients; 

% A. Nitti, Polytechnic University of Bari (2024)


% clc
clear 
close all
kf=1;


%% input parameters
Lr=1;                                                           % domain length
nr=17;                                                          % number of nodes  
refr='tanh-i';                                                  % grid stretching law

% prescribe analitycal solution and coefficients
qa=@(r) tanh(3.*r);         % analytical solution
rhs=@(r) -1.8.*tanh(3.*r).*sech(3.*r).^2-r.*tanh(3.*r)+3.*r.*sech(3.*r).^2;
coeffD=0.1;
coeffA=@(r) -1.*r;
coeffB=@(r) 1.*r;

% boundary conditions
% valr0= [1,0,0];                                               % values at r=0 boundary;
valr0= [0,1,3];                                                 % values at r=0 boundary;
valr1= [1,0,qa(Lr)];                                            % values at r=Lr boundary


% create and plot grid
gr=getgrid(Lr,nr,refr,false(1));
ndof=nr-1;


%% solve problem
% assemble coefficient matrix
A=getCoeffMat1D(gr,ndof,coeffA,coeffB,coeffD);

% assign boundary conditions
b=rhs(gr.xp');
[A,b]=bcs1D(A,b,gr,b,valr0,valr1,coeffA,coeffB,coeffD);

% compute matrix condition number
dA=decomposition(A);
tf = isIllConditioned(dA);
if (tf)
    error('ill-conditioned coefficient matrix:')
    rcond(dA)
end

% solve the tridiagonal problem
q=A\b;


%% compute error norm and plot results
figure(kf); kf=kf+1;
plot(gr.xp,qa(gr.xp),'o')
hold on
plot(gr.xp,q,'-x')
axis equal
legend('analytical','numerical','Location','southeast')
xlabel('x')
ylabel('q')


% compute rmse
err1=(qa(gr.xp')-q)./(qa(gr.xp')+1);
disp(['relative rmse: ',num2str(rms(err1,'all'))])
