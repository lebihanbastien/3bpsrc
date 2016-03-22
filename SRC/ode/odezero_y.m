%Event routine @ y = 0 (1/2 orbit)
function [value,isterminal,direction] = odezero_y(t,yvar)
% when value is equal to zero, an event is triggered.
% set isterminal to 1 to stop the solver at the first event, or 0 to
% get all the events.
%  direction=0 if all zeros are to be computed (the default), +1 if
%  only zeros where the event function is increasing, and -1 if only
%  zeros where the event function is decreasing.
value = yvar(2);  % when value = 0, an event is triggered
isterminal = 1; % terminate after the first event
direction = -1;  % get all the decreasing zeros
end