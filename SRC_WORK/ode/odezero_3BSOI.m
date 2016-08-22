function [value,isterminal,direction] = odezero_3BSOI(t, yvar)
% ODEZERO_PLANE event routine in MATLAB ODE format.
%
% [VALUE,ISTERMINAL,DIRECTION] = ODEZERO_PLANE(T, YVAR, EVENT) 
% is an event routine in MATLAB ODE format (see event in MATLAB help). 
% The condition of the event triggering is a given value of one dimension 
% in YVAR. This dimension is given by EVENT.DIM. 
% The value of YVAR(EVENT.DIM) that actually triggers the event is given by
% the scalar EVENT.value.
%
% See also EVENT
%
% BLB 2015

Xmoon  = 0.987849418376566;
r3bsoi = 0.414136897651021;

%Event parameters
value      =  sqrt((yvar(1)-Xmoon)^2 + yvar(2)^2 +yvar(3)^2) - r3bsoi;
isterminal =  1;
direction  =  0;

end