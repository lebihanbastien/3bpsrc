function mouseMove (event,data)
f = gcf;
Width = f.Position(3);
Height = f.Position(4);

C = f.CurrentPoint;

disp(['(X,Y) = (', num2str(C(1,1)/Width), ', ',num2str(C(1,2)/Height), ')']);
