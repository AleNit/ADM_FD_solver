
function kf=plotgrid(rn,tn,xn,Lr,Lx,kf)

nt=length(tn);
nx=length(xn);

X1=repmat(Lr.*cos(tn),[nx,1]);
Y1=repmat(Lr.*sin(tn),[nx,1]);
Z1=repmat(xn',[1,nt]);

[T,R]=meshgrid(tn,rn);
[X2,Y2]=pol2cart(T,R);
Z2=zeros(size(Y2));

Z3=Lx.*ones(size(Y2));

figure(kf); kf=kf+1;
surf(X1,Y1,Z1,'FaceColor',[0.9,0.9,0.9]); hold on
surf(X2,Y2,Z2,'FaceColor',[0.9,0.9,0.9])
surf(X2,Y2,Z3,'FaceColor',[0.9,0.9,0.9])
xlabel('x'); ylabel('y'); zlabel('z');
axis equal
drawnow

end