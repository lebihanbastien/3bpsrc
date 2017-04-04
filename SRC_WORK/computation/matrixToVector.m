function y = matrixToVector(y, M, nr, nc, shift)
% MATRIXTOVECTOR custom routine to set a matrix into a vector with a given
% shift.
%
% MATRIXTOVECTOR(Y, M, NR, NC, SHIFT) puts the NR x NC matrix M in the
% vector Y(SHIFT+1:end)
%
% See also VECTORTOMATRIX
%
% BLB 2015

y = matrixToVector_rc(y, M, nr, nc, shift, 'ROWWISE');

end


function y = matrixToVector_rc(y, M, nr, nc, shift, type)
% MATRIXTOVECTOR custom routine to set a matrix into a vector with a given
% shift.
%
% MATRIXTOVECTOR(Y, M, NR, NC, SHIFT) puts the NR x NC matrix M in the
% vector Y(SHIFT+1:end)
%
% See also VECTORTOMATRIX
%
% BLB 2015
switch(type)
    case('COLUMNWISE')
        y(shift+1:nr*nc+shift) = reshape(M, [nr*nc,1]);
    case('ROWWISE')
        for i = 1 : nr
            for j = 1 : nc
                m = nc*(i-1) + j;
                y(shift+m) = M(i,j);
            end
        end
end

end



