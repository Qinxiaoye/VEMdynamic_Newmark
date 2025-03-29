function x = NewMark(alpha,delta,K,M,C,F,sumTime,dt,ndim,ramp,dispT)
% New-Mark
% F is an external force and cannot change with time
% ramp == 0 is a step load, ramp == 1 is a ramp load, the initial value is 0 and the value is F at the end of the load step
% It can be improved so that the load application increases in a gradient manner
% dispT is the node number and degree of freedom number to be displayed

sumStep = fix(sumTime/dt);
sumNode = size(K,1)/ndim;

% initial condition
x = zeros(sumNode*ndim,sumStep+1);
xn = sparse(zeros(sumNode*ndim,1));
vn = xn;
% Acceleration
if ramp == 0
    an = M\(F-C*vn-K*xn);
else
    Fn = F*0;   
    an = M\(Fn-C*vn-K*xn);
end
% parameters
a0 = 1/(alpha*dt^2);
a1 = delta/(alpha*dt);
a2 = 1/(alpha*dt);
a3 = 1/(2*alpha)-1;
a4 = delta/alpha-1;
a5 = dt/2*(delta/alpha-2);
a6 = dt*(1-delta);
a7 = delta*dt;

% equivalent stiffness 
K = K + a0*M + a1*C;
% Cholesky (not necessary)
%[R,~] = chol(K);
% [L,U,P] = lu(K);

% ---------------------every time step------------------------------------
if size(dispT,1)>1
    dispT = dispT(1,:);
end
if size(dispT,1)>0
    figure;
    h = animatedline;
    axis auto;
    xlabel('Time');
    ylabel('Displacement');
    grid on;
    addpoints(h,0,0)
    drawnow
end

for n = 1:sumStep
    if ramp == 1
        Fn = F/sumStep*n;
        f = Fn + M*(a0*xn+a2*vn+a3*an)+C*(a1*xn+a4*vn+a5*an);  % 等效力
    else
        f = F + M*(a0*xn+a2*vn+a3*an)+C*(a1*xn+a4*vn+a5*an);  % 等效力
    end
    xnp = K\f;
    anp = a0*(xnp-xn)-a2*vn-a3*an;
    vnp = vn+a6*an+a7*anp;
    an = anp;
    vn = vnp;
    xn = xnp;
    x(:,n+1) = full(xn);
    if size(dispT,1)>0
        addpoints(h,n,x(dispT,n+1))
        drawnow
    end
end



