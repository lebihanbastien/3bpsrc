%-------------------------------------------------------------------------%
%   Set the vector y(shift+1:end) in the nrxnc matrix M
%
% Author: BLB
% Version: 1.0
% Year: 2015
%-------------------------------------------------------------------------%
function M = vectorToMatrix(y, nr, nc, shift)

M = eye(nr, nc);
for i = 1 : nr
    for j = 1 : nc
        m = nc*(i-1) + j;
        M(i,j) = y(shift+m);
    end
end

end

