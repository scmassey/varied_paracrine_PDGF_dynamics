% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = WoundBCs(xl,ul,xr,ur,t)

global EC50 Qr ninitx Dp K km DrBar RhoRBar p0 r0 S_max decay

pl = [0; 0];                               
ql = [1; 1];                                  
pr = pl;                            
qr = ql;                                  

