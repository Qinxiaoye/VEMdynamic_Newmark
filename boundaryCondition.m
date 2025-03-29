function [GK,GM] = boundaryCondition(GK,GM,fixNdf)

sumNdf = size(GK,1);

GK(fixNdf,:) = 0;
GK(:,fixNdf) = 0;
GK = GK+sparse(fixNdf,fixNdf,1,sumNdf,sumNdf);

GM(fixNdf,:) = 0;
GM(:,fixNdf) = 0;
GM = GM+sparse(fixNdf,fixNdf,1,sumNdf,sumNdf);