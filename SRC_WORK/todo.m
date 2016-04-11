%--------------------------------------------------------------------------
% TODO
%--------------------------------------------------------------------------
%
% 1. Careful with custom_odezero_2. It seems that there is a problem where
% the final time is reached and no event has been found (look at lfb for
% example).
%
% 2. Handle the issue with initplot2D. For now, only one libration point is
% dsiplayed when the 2D plot is initialized, which causes trouble when an
% orbit around another point is plotted after...