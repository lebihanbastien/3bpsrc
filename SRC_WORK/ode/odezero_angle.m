%-------------------------------------------------------------------------%
% Event routine in ode format. Stops at a given angle wrt a given center 
% (usually a primary or a center of mass).
%       
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function [value,isterminal,direction] = odezero_angle(t, yvar, event)

%Coordinate relative to the center
xl = yvar(1) - event.center(1);
yl = yvar(2) - event.center(2);

%Current angle
phi = atan2(yl,xl);

%Event parameters
value      = phi - event.value;
isterminal = event.isterminal;
direction  = event.direction; 

end