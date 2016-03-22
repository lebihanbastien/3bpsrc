%Flight path angle in radian
function fpa = flight_path_angle(yvc, center) 
    % Position wrt object's center
    erl = yvc(1:3) - center';
    
    %Calcul du sinus
    sinfpal = -  erl' * yvc(4:6)/ (norm(erl) * norm(yvc(4:6)));
    fpa = asin(sinfpal);  
end