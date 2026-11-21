function [jumps,R] = acceleration_jumps(S)
%ACCELERATION_JUMPS Actual covariant acceleration limits, not model residuals.
% Each limit is evaluated in its own interval. Both resulting vectors belong
% to T_{X_j}M and may be subtracted there without any extra transport.
M = S.M; 
n = numel(S.t)-1; 
jumps = zeros(n-1,1); 
R = cell(n-1,1);
for j = 2:n
    [~,~,al] = manifold_hermite(M,S.X{j-1},S.X{j},S.V{j-1},S.V{j}, ...
        S.t(j)-S.t(j-1),1);
    [~,~,ar] = manifold_hermite(M,S.X{j},S.X{j+1},S.V{j},S.V{j+1}, ...
        S.t(j+1)-S.t(j),0);
    R{j-1} = ar-al;
    jumps(j-1) = sqrt(max(0,M.inner(S.X{j},R{j-1},R{j-1})));
end
end
