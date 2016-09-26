function [coef,f,y] = SourceWoundPDEs(nx,nt,u,DuDx)

global EC50 Qr ninitx Dp K km DrBar RhoRBar p0 r0 S_max decay

r=u(1);
p=u(2);

DrDx=DuDx(1);
DpDx=DuDx(2);

coef = [1; 1]; 

% PDGF-dependent scaling terms:
Beta=p./(km+p);
gamma = EC50/(km+EC50);
BarFracR=(Beta)./(gamma+Beta);

f = [DrBar*(1-((r)/K)).*BarFracR*DrDx; %-nchi.*r.*BarFracR.*Dp*(1-(r/K))
    Dp*DpDx];

y = [RhoRBar*r*BarFracR*(1-(r/K));
    -Beta.*(Qr*r)+S_max*exp(-decay*nt)*(r-r0)/K];

