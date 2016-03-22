function exportFig(f,nom,type)

% Fonction pour sauvegarder vos figure en PDF.

% INPUTS
% f est le handle de la figure
% nom est un string, le nom a donner. Ne pas inclure d'extension !!

% Have fun !

% Emilien

if (nargin<3)
    type = 2;
end

if (abs(type)==1)
    nom = [nom '.eps'];
    snam='figureStyle';
    s=hgexport('readstyle',snam);
    s.Width = num2str(f.Position(3));
    s.Height = num2str(f.Position(4));
    
    hgexport(f,nom,s);
    
    code = ['!epstopdf ', nom];
    
    eval(code);
    if (type>0)
        delete(nom)
    end
else
    
    nom = [nom '.pdf'];
    snam='figureStyle';
    s=hgexport('readstyle',snam);
    s.Width = num2str(f.Position(3));
    s.Height = num2str(f.Position(4));
    
    
    f.PaperSize = [f.Position(3:4)];
    f.PaperPosition = [0 0 f.Position(3:4)];
    
    saveas(f,nom)
end