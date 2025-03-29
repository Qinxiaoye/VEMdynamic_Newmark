% Copyright:
% Bingbing Xu
% Contact:
% xubingbingdut@foxmail.com, bingbing.xu@ikm.uni-hannover.de
% Institute of Continuum Mechanics, Leibniz Universität Hannover
% reference: 
% https://github.com/Qinxiaoye/VEM_mechanics2D
% https://github.com/Terenceyuyue/mVEM (another way for solid
%            mechanics)
% also see
% https://www.sciencedirect.com/science/article/pii/S0045782525001653 

clear;

alpha = 0.25250625;  % parameters for Newmark method
delta = 0.505;
% time step 
sumTime = 1E-3;
dt = 1E-5;

load mesh.mat;

k = 1;

E = 1E3; nu = 0.3; % plane stress
rho = 8E-9;
mat.E = E; mat.nu = nu; mat.rho = rho;

sumNode = size(node,1);
% global stiffness matrix

[GK,GM] = globalK(node,elem,mat);


% right hand vector


ff = zeros(sumNode*2,1);
nodeD = find(node(:,2)<0.001);
nodeT = find(node(:,2)>19.999);

fixNdf = [nodeT;nodeT+sumNode];
face = findFace(node,elem,nodeD);
press = [face,ones(length(face),1)*(-100)];

nodeForce = getForce(node,elem,press,'y');


[GK,GM] = boundaryCondition(GK,GM,fixNdf);
ff = ff+nodeForce;
ff(fixNdf) = 0;
F = sparse(ff);

 %displace for point 5 in y direction
dispT = [5+sumNode];  


x = NewMark(alpha,delta,GK,GM,GK*0,F,sumTime,dt,2,0,dispT);

figure;
plot(0:dt:sumTime,x(dispT,:),'-o')

xend = x(:,end);
ux = xend(1:sumNode);
uy = xend(sumNode+1:sumNode*2);
figure;
showsolution(node+[ux,uy],elem,uy);






