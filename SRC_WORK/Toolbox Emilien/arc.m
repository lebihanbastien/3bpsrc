function [arc,data,ind1,ind2] = arc(centre,r,vect1,vect2,depass)

axe = cross(vect1,vect2,2);
axe = axe/norm(axe);

ang = 180/pi*acos(dot(vect1,vect2)/(norm(vect1)*norm(vect2)));

t = unique([-depass:1:0 1:1:ang-0.00001 ang:1:depass+ang])'/180*pi;

XYZ = ones(size(t,1),1)*(r*vect1/norm(vect1));
for i=1:size(t,1)
    XYZ(i,:) = centre+(vrrotvec2mat([axe(1) axe(2) axe(3) t(i)])*XYZ(i,:)')';
end

arc = plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3));

ind1 = length(-depass:1:0)+1;
ind2 = length([-depass:1:0 1:1:ang-0.00001])+1;
data = XYZ;
