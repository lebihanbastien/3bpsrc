function M = vectorToMatrix(y, nr, nc, shift)
% VECTORTOMATRIX custom routine to set a vector into a matrix with a given
% shift.
%
% VECTORTOMATRIX(Y, M, NR, NC, SHIFT) puts the vector Y(SHIFT+1:end) in the
% NR x NC matrix M.
%
% See also MATRIXTOVECTOR
%
% BLB 2015

M = vectorToMatrix_rc(y, nr, nc, shift, 'ROWWISE');

end


function M = vectorToMatrix_rc(y, nr, nc, shift, type)
% VECTORTOMATRIX custom routine to set a vector into a matrix with a given
% shift.
%
% VECTORTOMATRIX(Y, M, NR, NC, SHIFT) puts the vector Y(SHIFT+1:end) in the
% NR x NC matrix M.
%
% See also MATRIXTOVECTOR
%
% BLB 2015

switch(type)
    case('COLUMNWISE')
        M = reshape(y(shift+1:nr*nc+shift), [nr, nc]);
    case('ROWWISE')
        M = eye(nr, nc);
        for i = 1 : nr
            for j = 1 : nc
                m = nc*(i-1) + j;
                M(i,j) = y(shift+m);
            end
        end
end




end