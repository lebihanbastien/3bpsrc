%-------------------------------------------------------------------------%
%   Set the the nrxnc matrix M in the vector y(shift+1:end) 
%-------------------------------------------------------------------------%
function y = matrixToVector(y, M, nr, nc, shift)

for i = 1 : nr
    for j = 1 : nc
        m = nc*(i-1) + j;
        y(shift+m) = M(i,j);
    end
end

end


