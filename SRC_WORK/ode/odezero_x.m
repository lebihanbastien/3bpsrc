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
function [value,isterminal,direction] = odezero_x(t,yvar)
%Event parameters
value = yvar(1);
isterminal = 1;
direction =  0;
end