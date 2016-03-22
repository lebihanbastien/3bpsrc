%-------------------------------------------------------------------------%
% Event routine for ode. Stops at a given angle wrt a given center (usually
% the Earth, or the Moon, or any primary).
%-------------------------------------------------------------------------%
function [value,isterminal,direction] = odezero_angle(t, yvar, event)

%Coordinate relative to the center
xl = yvar(1) - event.center(1);
yl = yvar(2) - event.center(2);

%phi = pi + 2 * atan(yl/ (xl + sqrt(xl^2 + yl^2)));
phi = atan2(yl,xl);
value = phi - event.value;  % when value = 0, an event is triggered
isterminal = event.isterminal; % terminate after the first event
direction =  event.direction;  % get all the increasing zeros
end