clear all
close all
clc

addpath Utilities/

f = figure('Color',[1 1 1]);
set(gcf,'units','centimeters','position',[10 10 10 8])

% view:
a = 0;
b = -30;
h = -30;

% ------- DEBUT DU GLOBE ------------------------------

r = .9;
Nlo = 8;
Nla = 11;
xc = [];
yc = [];
zc = [];
for lo = ((0:360/Nlo:359.999)-h)/180*pi;
    t=(0:1:360)'/180*pi;
    
    auxX = r*cos(t);
    auxZ = r*sin(t);
    auxY = r*zeros(size(t));
    
    auxXYZ = [auxX auxY auxZ];
    for ind = 1:size(auxXYZ,1)
        auxXYZ(ind,:) = (vrrotvec2mat([0 0 1 lo])*auxXYZ(ind,:)')';
    end
    xc = [xc; NaN; auxXYZ(:,1)];
    yc = [yc; NaN; auxXYZ(:,2)];
    zc = [zc; NaN; auxXYZ(:,3)];
end
for la = (-90:180/(Nla-1):90)/180*pi
    auxX = r*cos(la)*sin(t);
    auxY = r*cos(la)*cos(t);
    auxZ = r*sin(la)*ones(size(t));
    
    xc = [xc; NaN; auxX];
    yc = [yc; NaN; auxY];
    zc = [zc; NaN; auxZ];
end

XYZc = rotation([xc yc zc],vrrotvec2mat([0 0 1 -b/180*pi]));
cercles = plot3(XYZc(:,1),XYZc(:,2),XYZc(:,3),'k');
cercles.Color = 0*[1 1 1];
cercles.LineWidth = .3;
hold on
axis equal

tcont=(0:.2:360)'/180*pi;
auxX = r*sin(tcont);
auxY = 0*tcont+0.01;
auxZ = r*cos(tcont);
auxXYZ = [auxX auxY auxZ];
for ind = 1:size(auxXYZ,1)
    auxXYZ(ind,:) = (vrrotvec2mat([(vrrotvec2mat([0 0 1 a/180*pi])*[1 0 0]')' -b/180*pi])*vrrotvec2mat([0 0 1 a/180*pi])*auxXYZ(ind,:)')';
end
xc = auxXYZ(:,1);
yc = auxXYZ(:,2);
zc = auxXYZ(:,3);
cont = fill3(xc,yc,zc,'White');
cont.LineWidth = 1;
cont.EdgeColor = 'k';



auxX = r*sin(t);
auxY = r*cos(t);
auxZ = r*zeros(size(t));
auxXYZ = [auxX auxY auxZ];
xc = auxXYZ(:,1);
yc = auxXYZ(:,2);
zc = auxXYZ(:,3);
equat = plot3(xc,yc,zc);
equat.LineWidth = 1;
equat.Color = 'k';

load coast
lat = lat/180*pi;
long = -(long-h)/180*pi;
x = r*cos(lat).*sin(long);
y = r*cos(lat).*cos(long);
z = r*sin(lat);
cote = plot3(x,y,z);
cote.Color = 0*[1 1 1];
cote.LineWidth = .5;


view(a,b)


%--------------- Orbite -----------------------

t = (0:1:2.75*360)'/180*pi;
r0 = 1.2;
rho = r0.*(exp(0.5*t/(2*pi)));

dt = 0;
XYZorb = [rho.*cos(t-dt) rho.*sin(t-dt) 0*t];

i1 = 10;
i2 = -20;
XYZorb = rotation(XYZorb,vrrotvec2mat([1 0 0 -b/180*pi])*vrrotvec2mat([1 0 0 90/180*pi]));

orbite = plot3(XYZorb(:,1),XYZorb(:,2),XYZorb(:,3));
orbite.Color = 'k';



%--------------- Vecteurs -----------------------
npt = 441;
pt = XYZorb(npt,:);
pt0 = XYZorb(1,:);

xyz = [0 0 0;pt];
arrow3(xyz(1,:),xyz(2,:),'k',1.5,.2,.3);
ptT1 = pt;

xyz = [0 0 0;pt0];
arrow3(xyz(1,:),xyz(2,:),'k',1.5,.2,.3);

xyz = [pt0*1.05;1.4*pt0];
Vx1 = plot3(xyz(:,1),xyz(:,2),xyz(:,3));
Vx1.Color = 'k';
Vx1.LineStyle = '--';
Vx1.LineWidth = 1;
[arc0,data,~,ind2] = arc([0 0 0],1.6,pt0,pt,0);
arc0.Color = 'k';
arc0.LineStyle = '-';
ptaf = 1.6*pt/norm(pt);
arrow3(data(ind2-1,:),data(ind2,:),'k',1,.15,.3);

aux = 4*(XYZorb(npt+1,:)-XYZorb(npt-1,:))/norm(XYZorb(npt+1,:)-XYZorb(npt-1,:));
xyz = [pt;pt+aux];
arrow3(xyz(1,:),xyz(2,:),'k',1.5,.2,.3);
ptT2 = pt+aux;
xyz = [pt;pt+aux*0.4];
arrow3(xyz(1,:),xyz(2,:),'k',1.5,.2,.3);
ptT3 = pt+aux*0.4;

aux = cross(cross(pt,aux),pt);
aux = 3*aux/norm(aux);
xyz = [pt;pt+aux];
[ligne, pointe] = arrow3(xyz(1,:),xyz(2,:),'k',1,.2,.3);
ligne.LineStyle = '--';
ptT4 = pt+aux;

[arc,data,ind1,ind2] = arc(pt,2,ptT2-pt,ptT4-pt,8);
ind1=ind1-1;
ind2=ind2+1;
arc.Color = 'k';
arrow3(data(ind1-1,:),data(ind1,:),'k',1,.15,.3);
arrow3(data(ind2,:),data(ind2-1,:),'k',1,.15,.3);
ptT5 = data(ind1,:);


npt = 600;
arrow3(XYZorb(npt-1,:),XYZorb(npt,:),'k',1.5,.2,.3);
pto = XYZorb(npt,:);


text(ptT1(1)-.35,ptT1(2),ptT1(3)-0.15,'$\mathbf{r}$','Interpreter','latex','Fontsize',12);
text(pt0(1)-.1,pt0(2),pt0(3)-0.25,'$\mathbf{r}_0$','Interpreter','latex','Fontsize',12);
text(ptT2(1),ptT2(2),ptT2(3)+.25,'$\mathbf{v}$','Interpreter','latex','Fontsize',12);
text(ptT3(1)+.15,ptT3(2),ptT3(3)+.2,'$\mathbf{F}_{C_{th}}$','Interpreter','latex','Fontsize',12);
text(ptT4(1),ptT4(2),ptT4(3)-.35,'$\hat{\mathbf{x}}_O$','Interpreter','latex','Fontsize',12);
text(ptT5(1)+.05,ptT5(2),ptT5(3)-.45,'$\varphi$','Interpreter','latex','Fontsize',12);
text(ptaf(1)+.3,ptaf(2),ptaf(3)+.11,'$\nu$','Interpreter','latex','Fontsize',12);
text(pto(1)-1.2,pto(2),pto(3)+.1,'orbit','Interpreter','latex','Fontsize',12);

grid off
axis off

% 
zoom(1.5)
ax = gca;
ax.Position = ax.Position+[0 -.15 0 0];


% grid on
% axis on
% xlabel('x')
% ylabel('y')
% zlabel('z')


%%

exportFig(f,'Spiral_orbit')
