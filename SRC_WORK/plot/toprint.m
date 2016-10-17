% To print the current (last) figure
% Once you think that the figure looks good, just change the name of the
% filename.
%
% BLBL 2016

filename = 'plot/J_earth_moon';
fig = gcf;

%% Size
%fig.Position = [2046 55 981 826];

%% Margins
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

%% Change the orientation of the 3D plot, if it exists
%view([70 20]);
%view([0 90]);

%% Background
set(gcf, 'Color', 'none');

%% Print
% .fig
savefig(fig, [filename, '.fig'])

% .pdf direct
print(fig,filename,'-dpdf')

% .eps then .pdf
%print(fig,filename,'-depsc')
%seval(sprintf('!epstopdf %s.eps', filename)); 

%% Print .png
% Native matlab
% print(fig,filename,'-dpng')

% -png with no background
%set(gca, 'Color', 'w');
eval(sprintf('export_fig %s.png', filename));

% Transparency (not working)
% [A, bitmapData, transparency] = imread(sprintf('%s.png', filename));
% imwrite(A, sprintf('%s.png', filename), 'png', 'Alpha', 0);