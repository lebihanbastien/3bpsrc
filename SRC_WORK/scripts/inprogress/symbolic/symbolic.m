%-------------------------------------------------------------------------%
% Symbolic derivatives of the BCFBP potential (sun-perturbed Earth-Moon
% problem
%-------------------------------------------------------------------------%
%init
x = sym('x');
y = sym('y');
z = sym('z');
theta = sym('theta');
ms = sym('ms');
as = sym('as');
mu = sym('mu');
%Radii
r1 = sqrt((x+mu)^2 + y^2 + z^2);
r2 = sqrt((x-1+mu)^2 + y^2 + z^2);
rs = sqrt((x-as*cos(theta))^2 + (y-as*sin(theta))^2 + z^2);
%Potential U4
u_4 = 1/2*(x^2+y^2) + (1-mu)/r1 + mu/r2 + ms/rs - ms/as^2*(x*cos(theta) + y*sin(theta));
%Order1
du_x = diff(u_4,x);
du_y = diff(u_4,y);
du_z = diff(u_4,z);
%Order 2
du_xx = diff(du_x,x)
du_xy = diff(du_x,y)
du_xz = diff(du_x,z)

du_yx = diff(du_y,x)
du_yy = diff(du_y,y)
du_yz = diff(du_y,z)

du_zx = diff(du_z,x)
du_zy = diff(du_z,y)
du_zz = diff(du_z,z)

