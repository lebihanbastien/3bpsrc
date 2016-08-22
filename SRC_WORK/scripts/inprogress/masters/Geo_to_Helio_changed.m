function yarc_bcp_B = Geo_to_Helio (yarc_bcp,tarc_bcp,theta_sun,cr3bp,cst)

yarc_bcp_B = zeros (size(yarc_bcp,1),6);
for i=1:size(yarc_bcp,1)
        %     theta_0=0.0;
        %     T = tarc_bcp(i) + theta_0;
        %     R_1 = [cos(T),-sin(T);sin(T),cos(T)];
        %     R_21 = [-sin(T),-cos(T);cos(T),-sin(T)];
        %     R_matrix = [R_1 , zeros(2); R_21 , R_1];
        %     yarc_bcp_In = (R_matrix*([yarc_bcp(i,1:2),yarc_bcp(i,4:5)]'))';%+[cr3bp.mu, 0,0,0]'))';
        %     L_AB = 1/(cst.sun.as);
        %     T_AB = 0.074801112292747;%30/360;%-cst.sun.omegaS;
        % 
        %     yarc_bcp_In_B = [yarc_bcp_In(1:2),yarc_bcp_In(4:5)/T_AB]*L_AB;
        %     T_B = T_AB*T+degtorad(theta_sun);
        %     R_1 = [cos(T_B),-sin(T_B);sin(T_B),cos(T_B)];
        %     R_21 = [-sin(T_B),-cos(T_B);cos(T_B),-sin(T_B)];
        %     R_matrix = [R_1 , zeros(2); R_21 , R_1];
        %     yarc_bcp_B(i,:) = (inv (R_matrix)*yarc_bcp_In_B')' -[0.999996959789866, 0, 0, 0];
        
            L_AB = 1/(cst.sun.as);
        	T_AB = 0.074801112292747;%30/360;%-cst.sun.omegaS;
            
            T = tarc_bcp(i);
            T_B = -0.925195985520347*T+degtorad(theta_sun);
            
%             as = 3.888111400000000e+02;
%             yarc_bcp(1) = 1-cr3bp.mu;%as*cos(T_B);
%             yarc_bcp(2) = 0;%as*sin(T_B);
            
            yarc_bcp_In_B = [yarc_bcp(i,1:3), yarc_bcp(i,4:6)];
            
            R_1 = [cos(T_B),-sin(T_B), 0 ;sin(T_B),cos(T_B), 0; 0 0 1];
            R_21 = [-sin(T_B),-cos(T_B), 0;cos(T_B),-sin(T_B), 0; 0, 0, 0];
            R_matrix = [R_1 , zeros(3); R_21 , R_1];
            yarc_bcp_B(i,:) = (inv(R_matrix)*yarc_bcp_In_B')';
end

end

