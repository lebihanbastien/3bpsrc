function ysg  = selenocentric2selenographic(ysc, cr3bp, t)
%
% YSG  = SELENOCENTRIC2SELENOGRAPHIC(YSC, CR3BP) computes the coordinates 
% YSC into selenographic coordinates.
%
% BLB 2016

% Selenocentric to synodical
ysyn = selenocentric2synodical(ysc, cr3bp, t);

% Synodical to selenographic
ysg = synodical2selenographic(ysyn, cr3bp);

end
