function manifold_branch = init_manifold_branch(stability, way, varargin)
% INIT_MANIFOLD_BRANCH initializes a manifold leg.
%
% MANIFOLD_BRANCH = INIT_MANIFOLD_BRANCH(STABILITY, WAY) initializes a
% manifold leg with a given STABILITY (STABLE/UNSTABLE) and a given WAY
% (EXTERIOR/INTERIOR).
%
% MANIFOLD_BRANCH = INIT_MANIFOLD_BRANCH(STABILITY, WAY, EVENT) initializes
% a manifold leg with a given STABILITY (STABLE/UNSTABLE), a given WAY
% (EXTERIOR/INTERIOR), and a given EVENT structure (see INIT_EVENT.m).
%
% See also INIT_EVENT
%
% BLB 2016.
        
%--------------------------------------------------------------------------
switch(nargin)
    
    case 2
        %------------------------------------------------------------------
        % if nargin == 2 we do not have any additionnal inputs in varargin.
        % There is no event defined and the manifold is let free to be
        % integrated indefinitely.
        %------------------------------------------------------------------
        %Stability (stable or unstable)
        manifold_branch.stability = stability;
        %Exterior or interior
        manifold_branch.way = way;
        %Event type during integration
        manifold_branch.event.type = 'FREE';
        
    case 3
        %------------------------------------------------------------------
        % if nargin == 3 we have the following inputs in varargin:
        % 1. event (structure).
        %------------------------------------------------------------------
        %Stability (stable or unstable)
        manifold_branch.stability = stability;
        %Exterior or interior
        manifold_branch.way = way;
        %Event during integration
        manifold_branch.event = varargin{1};
        
    otherwise
        error('Wrong number of inputss');
end

end

