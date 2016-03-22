%Null terrestrial flight path angle
function [value,isterminal,direction] = odezero_flightpathangle(t,yvar, event)
% Position with respect to the center of the Earth
xrt = yvar(1:3) - event.center';
% Velocity wrt the Earth
vrt = yvar(4:6)';
value = - vrt * xrt - event.value;  % when value = event.value, an event is triggered
isterminal = 1;%event.isterminal; % terminate after the first event
direction = 1;%event.direction;  % get all the  zeros
