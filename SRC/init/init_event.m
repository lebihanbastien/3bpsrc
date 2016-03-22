%-------------------------------------------------------------------------%
% Initializes an event structure
%-------------------------------------------------------------------------%
function [ event ] = init_event(type, value, isterminal, direction, center, cst)

%Event type during integration
event.type = type;
%Definition parameter for the event plane
event.value = value;
%Is the event terminal?
event.isterminal = isterminal;
%Dimension: dim = 1 for x, dim = 2 for y, dim = 0 if there is no event set
switch(event.type)
    case cst.manifold.event.type.X_SECTION
        event.dim = 1;
    case cst.manifold.event.type.Y_SECTION
        event.dim = 2;
    case cst.manifold.event.type.Z_SECTION
        event.dim = 3;
    case cst.manifold.event.type.FREE
        event.dim = 0;
end

%Decreasing (-1), increasing (1) or all (0) zeros;
event.direction = direction;
%Center (if the event is an angle event e.g. for lunar flyby detection).
event.center = center;

end

