%Event routine @ y = 0 (1/2 orbit)
function [value,isterminal,direction] = odezero_plane(t, yvar, event)
% when value is equal to zero, an event is triggered.
% set isterminal to 1 to stop the solver at the first event, or 0 to
% get all the events.
%  direction=0 if all zeros are to be computed (the default), +1 if
%  only zeros where the event function is increasing, and -1 if only
%  zeros where the event function is decreasing.
value = yvar(event.dim)-event.value;  % when value = 0, an event is triggered
isterminal = event.isterminal; % terminate after the first event
direction =  event.direction;  % get all the increasing zeros
end