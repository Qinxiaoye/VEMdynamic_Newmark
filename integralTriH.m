function Int = integralTriH(n,nodeT,elemT,k,hK,c)
%integralTri approximates integrals in a polygonal domain with trianguation (nodeT,elemT).


[lambda,weight] = quadpts(n);
NT = size(elemT,1); 
if  k == 1
    Int = zeros(3,3);
elseif k == 2
    Int = zeros(6,6);
elseif k == 3
    Int = zeros(10,10);
end

for iel = 1:NT
    vT = nodeT(elemT(iel,:),:);
    area = 0.5*abs(det([[1;1;1],vT]));
    
    xy = lambda*vT;
    for p = 1:size(xy,1)
        f = myM([xy(p,1),xy(p,2)],k,hK,c);
        Int = Int + area*weight(p)*(f'*f);
    end
end