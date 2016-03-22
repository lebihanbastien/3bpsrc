function [X,Y,Z] = rotateSurf(x,y,z,P)

siz = size(x);
x = reshape(x,[prod(siz) 1]);
y = reshape(y,[prod(siz) 1]);
z = reshape(z,[prod(siz) 1]);

aux = rotation([x y z],P);
X = aux(:,1);
Y = aux(:,2);
Z = aux(:,3);

X = reshape(X,siz);
Y = reshape(Y,siz);
Z = reshape(Z,siz);