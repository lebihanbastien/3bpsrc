%-------------------------------------------------------------------------%
% Initializes a manifold branch (with event) prior to real computation
%-------------------------------------------------------------------------%
function [ manifold_branch ] = init_manifold_branch_event(stability, way, varargin)

if(nargin == 3)
    %----------------------------------------------
    %if nargin == 3 we have the following inputs in varargin:
    % 1. event
    %----------------------------------------------
    %Stability (stable or unstable)
    manifold_branch.stability = stability;
    %Exterior or interior
    manifold_branch.way = way;
    %Event during integration
    manifold_branch.event = varargin{1};
    
elseif(nargin == 8)
    %----------------------------------------------
    % if nargin == 8 we have the following inputs in varargin, in this order:
    % 1. event_type
    % 2. event_value
    % 3. event_isterminal
    % 4. event_direction
    % 5. event_center
    % 6. cst
    %----------------------------------------------
    cst = varargin{6};
    %Stability (stable or unstable)
    manifold_branch.stability = stability;
    %Exterior or interior
    manifold_branch.way = way;
    %Event type during integration
    manifold_branch.event.type = varargin{1};
    %Definition parameter for the event plane
    manifold_branch.event.value = varargin{2};
    %Is the event terminal?
    manifold_branch.event.isterminal = varargin{3};
    %Dimension: dim = 1 for x, dim = 2 for y, dim = 0 if there is no event set
    switch(varargin{1})
        case cst.manifold.event.type.X_SECTION
            manifold_branch.event.dim = 1;
        case cst.manifold.event.type.Y_SECTION
            manifold_branch.event.dim = 2;
        case cst.manifold.event.type.Z_SECTION
            manifold_branch.event.dim = 3;
        case cst.manifold.event.type.FREE
            manifold_branch.event.dim = 0;
    end
   
    %Decreasing (-1), increasing (1) or all (0) zeros;
    manifold_branch.event.direction = varargin{4};
    
    %Center (if the event is an angle event e.g. for lunar flyby detection).
    manifold_branch.event.center = varargin{5};
    

else
    error('wrong number of arguments');
end

end




