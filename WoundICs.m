% --------------------------------------------------------------------------

function u0 = WoundICs(x)

global EC50 Qr ninitx Dp K km DrBar RhoRBar p0 r0 S_max decay

u0(1)=r0;

% if within initial core:
if x < ninitx 
    u0(2)=p0;
else
    u0(2)=0;
end

u0=u0';
