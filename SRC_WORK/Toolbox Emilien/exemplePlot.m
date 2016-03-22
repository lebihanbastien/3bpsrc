clear all
close all
clc



% données d'exemple .. -------------------------------------------------
mu = 3.986e14;
rmin = (6378+600)*1e3;
i=0;
for rmax = ([1125 2250 4500 9000 18000 36000]+6378)*1e3;
i=i+1;
a = (rmin+rmax)/2;
e = (rmax-rmin)/(rmax+rmin);

theta = (0:1:360)'/180*pi;

omega(:,i) = sqrt(mu/(a^3*(1-e^2)^3))*(1+e*cos(theta)).^2;
eta(:,i) = (mu/(a^3*(1-e^2)^3))*(1+e*cos(theta)).^3.*(-2*e*sin(theta));

end
% fin données d'exemple ------------------------------------------------




f1 = figure('Color',[1 1 1]);
set(gcf,'units','centimeters','position',[10 10 12 8]);
% les deux dernier nombre sont la taille en cm. Important, c'est ICI que
% l'on choisit la taille que la figure aura sur le papier...
% (les deux premier nombre sont la position sur l'écran, donc on s'en cogne.



% ma figure d'exemple
hold on
f1.CurrentAxes.ColorOrder = ((0.2:(1/0.741-0.2)/7:1/0.741))'*[0 0.4470 0.7410];
plot(theta,omega);
hold off
aux = f1.CurrentAxes.YAxis.Limits;
f1.CurrentAxes.XAxis.TickValues = (0:0.5:2)*pi;
f1.CurrentAxes.XAxis.TickLabels = {'$0$'; '$\pi/2$'; '$\pi$'; '$3\pi/2$'; '$2\pi$'};
axis tight
f1.CurrentAxes.YAxis.Limits = aux;
box off
grid on
xlabel('$\nu$ (rad)') % ecrire les légende en latex. Elles seront mises a la bonne polie et taille par la suite.
ylabel('$\omega$ (rad s$^{-1}$)')



f1 = prepareFig(f1); % cette fonction permet de définir les tailles des textes de la fig. (en pt de taille latex)




set(f1, 'WindowButtonMotionFcn', @mouseMove); % pour pouvoir positionner correctement des annotations ; ballade ta souris surla figure.


annotation(f1,'textarrow',[0.8 0.66],[0.66 0.664],...
    'String','$e = 0.036$','Interpreter','latex','Fontsize',12);
annotation(f1,'textarrow',[0.8 0.66],0.588*[1 1],...
    'String','$e = 0.106$','Interpreter','latex','Fontsize',12);
annotation(f1,'textarrow',[0.8 0.66],0.495*[1 1],...
    'String','$e = 0.218$','Interpreter','latex','Fontsize',12);
annotation(f1,'textarrow',[0.8 0.66],0.40*[1 1],...
    'String','$e = 0.376$','Interpreter','latex','Fontsize',12);
annotation(f1,'textarrow',[0.8 0.66],0.32*[1 1],...
    'String','$e = 0.555$','Interpreter','latex','Fontsize',12);
annotation(f1,'textarrow',[0.8 0.66],0.27*[1 1],...
    'String','$e = 0.717$','Interpreter','latex','Fontsize',12);


% export en pdf. La taille est celle définie en créant la figure.
exportFig(f1,'omegaNu')


