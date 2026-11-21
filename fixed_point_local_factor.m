function [qbound,A] = fixed_point_local_factor(S)
%FIXED_POINT_LOCAL_FACTOR Measured LOCAL derivative bound; solver does not use it.
% Phi'=-K^{-1}N'. The full N has a linear previous-velocity block and a
% quadratic current-velocity block, so only TWO blocks per row are needed.
% Polarization evaluates B[w,z] exactly in exact arithmetic. This avoids a
% global finite-difference Jacobian and division by tiny increments.
% qbound is a measured bound at this point, not a certificate over a ball.
M = S.M;
n = numel(S.t)-1;
d = M.dim;
m = d*(n-1);
DN = zeros(m);
for j = 2:n
    i = j-1;
    p = S.X{i};
    q = S.X{j};
    dt = S.t(j)-S.t(i);
    u = M.log(p,q);
    w = M.dlog(p,q,S.V{j});
    row = (j-2)*d+(1:d);
    for l = 1:d
        z = M.dlog(p,q,S.basis{j}{l});
        [~,~,bp] = M.jet(p,u,w+z,0*u);
        [~,~,bm] = M.jet(p,u,w-z,0*u);
        value = (bp-bm)/4; % derivative of B[w,w]/2 is B[w,z].
        for k = 1:d
            DN(row(k),row(l)) = M.inner(q,S.basis{j}{k},value);
        end
        if i > 1
            e = S.basis{i}{l};
            [~,Le,~] = M.jet(p,u,e,0*u);
            value = (Le-M.dlog(q,p,e))/dt;
            col = (i-2)*d+l;
            for k = 1:d
                DN(row(k),col) = M.inner(q,S.basis{j}{k},value);
            end
        end
    end
end
A = -solve_velocity_system(S,DN);
qbound = 0;
for j = 1:n-1
    value = 0;
    for l = 1:n-1
        value = value+norm(A((j-1)*d+(1:d),(l-1)*d+(1:d)),2);
    end
    qbound = max(qbound,value);
end
end
