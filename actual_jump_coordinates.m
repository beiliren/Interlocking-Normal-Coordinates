function [r,maximum,jumps] = actual_jump_coordinates(S)
%ACTUAL_JUMP_COORDINATES Coordinates of a(t_j^+) - a(t_j^-).
% E_j is orthonormal, so the maximum block 2-norm equals the intrinsic jump.
[jumps,R] = acceleration_jumps(S);
n = numel(S.t)-1;
d = S.M.dim;
r = zeros(d*(n-1),1);
for j = 2:n
    row = (j-2)*d+(1:d);
    for l = 1:d
        r(row(l)) = S.M.inner(S.X{j},S.basis{j}{l},R{j-1});
    end
end
maximum = max([0;jumps]);
end
