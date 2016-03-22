function f = prepareFig(f)
axe = f.CurrentAxes;
axe.Title.Interpreter = 'latex';
axe.Title.FontSize = 12;
axe.XLabel.Interpreter = 'latex';
axe.YLabel.Interpreter = 'latex';
axe.Box = 'off';
axe.TickLabelInterpreter = 'latex';
axe.FontSize = 12;
axe.XLabel.FontSize = 12;
axe.YLabel.FontSize = 12;
grid on
axe.GridLineStyle = ':';
axe.GridAlpha = .5;