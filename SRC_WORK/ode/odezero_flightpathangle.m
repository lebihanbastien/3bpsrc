%-------------------------------------------------------------------------%
% Event routine in ode format. Stops at a null terrestrial flight path angle
%       
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function [value,isterminal,direction] = odezero_flightpathangle(t,yvar, event)
% Position with respect to the center of the Earth
xrt = yvar(1:3) - event.center';

% Velocity wrt the Earth
vrt = yvar(4:6)';

%Event parameters
value = - vrt * xrt - event.value;
isterminal = event.isterminal;
direction  = event.direction;

end
