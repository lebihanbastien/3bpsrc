function [] = figtoprint(fig, filename, cview)
% figtoprint(fig, filename, cview) prints the figure in the handle fig, in
% the following format: filename.fig, .eps, .pdf, .png.

%% Size
fig.Position = [2046 55 981 826];

%% Margins
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

%% Change the orientation of the plot
view(cview);

%% Print
% .fig
savefig(fig, [filename, '.fig'])

% .pdf direct
%print(fig,filename,'-dpdf')

% .eps then .pdf
print(fig,filename,'-depsc')
eval(sprintf('!epstopdf %s.eps', filename)); 

% .png
print(fig,filename,'-dpng')

end
