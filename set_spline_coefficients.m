function S = set_spline_coefficients(S,c)
%SET_SPLINE_COEFFICIENTS Update interior velocities in the fixed tangent bases.
% Positions, bases, K, LU factors, and both endpoint velocities stay fixed.
n = numel(S.t)-1;
d = S.M.dim;
if numel(c) ~= d*(n-1) 
    error('Incorrect coefficient-vector length.'); 
end
S.c = c(:);
% Preserve the initializer diagnostic separately, while keeping this field
% accurate after a correction. It is not the actual acceleration residual.
S.relative_residual = norm(S.matrix*S.c-S.rhs,inf)/(1+norm(S.rhs,inf));
for j = 2:n
    row = (j-2)*d+(1:d);
    S.V{j} = zeros(size(S.X{j}));
    for l = 1:d
        S.V{j} = S.V{j}+S.c(row(l))*S.basis{j}{l};
    end
end
end
