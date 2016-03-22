%-------------------------------------------------------------------------%
% Event routine in ode format. Stops on a simple plane defined by:
%       yvar(1)=0
% 
%   Note that this routine is already included in odezero_plane
%
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function [value,isterminal,direction] = odezero_y(t,yvar)
%Event parameters
value = yvar(2);
isterminal = 1;
direction =  -1;
end