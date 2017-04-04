function  ysc  = selenographic2selenocentric(ysg, cr3bp, t)
%
% YSC  = SELENOGRAPHIC2SELENOCENTRIC(YSG, CR3BP) computes the coordinates 
% YSG into selenocentric coordinates.
%
% BLB 2016

% Selenographic 2 synodical
ysyn = selenographic2synodical(ysg, cr3bp);

% Synodical to selenocentric
ysc = synodical2selenocentric(ysyn, cr3bp, t);

end
