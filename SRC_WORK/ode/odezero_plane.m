%-------------------------------------------------------------------------%
% Event routine in ode format. Stops on a simple plane defined by:
%       yvar(event.dim)=event.value
%       
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function [value,isterminal,direction] = odezero_plane(t, yvar, event)
%Event parameters
value = yvar(event.dim)-event.value;
isterminal = event.isterminal;
direction =  event.direction;
end