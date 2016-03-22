%-------------------------------------------------------------------------%
% init_manifold_branch( stability, way, cst)
% 
% Initializes a manifold branch (no event)
%       
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
% @return the structure manifold_branch
function [ manifold_branch ] = init_manifold_branch( stability, way, cst)
% Initializes the manifold_branch prior to real computation
%--------------------------------------------------------------------------

%Stability (stable or unstable)
manifold_branch.stability = stability;
%Exterior or interior
manifold_branch.way = way;
%Event type during integration
manifold_branch.event.type = cst.manifold.event.type.FREE;

end

