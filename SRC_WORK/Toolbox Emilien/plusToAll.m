function B = plusToAll(A,a)

B = zeros(size(A));
for i = 1:size(A,2)
    B(:,i) = A(:,i)+a(i);
end