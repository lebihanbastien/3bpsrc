%-------------------------------------------------------------------------%
% flight_path_angle(yvc, center) 
%
% Flight path angle in radian. If the state in position/velocity is (x, v),
% the sinus of the flight path angle with respect to the center x0 
% is given by:
%
%   sin(fpa)  = - (x - xO).v / (|x - x0|*|v|)
%
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function fpa = flight_path_angle(yvc, center) 
    % Position wrt object's center
    erl = yvc(1:3) - center';
    
    %Calcul du sinus
    sinfpal = -  erl' * yvc(4:6)/ (norm(erl) * norm(yvc(4:6)));
    fpa = asin(sinfpal);  
end