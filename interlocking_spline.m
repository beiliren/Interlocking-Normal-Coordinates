function S = interlocking_spline(M,t,X,v0,vn)
%INTERLOCKING_SPLINE Assemble and factor K; compute the fixed-point initializer.
% t: strictly increasing vector; X: cell array of n+1 manifold points.
% v0, vn: tangent endpoint velocities. S.V are the computed knot velocities.
% Unknowns are coefficients in intrinsic tangent bases. Consequently every
% computed velocity is tangent without a projection or optimization step.
t = t(:); 
n = numel(t)-1; 
d = M.dim;
if numel(X) ~= n+1 || n < 1 || any(diff(t)<=0)
    error('Inconsistent interpolation data.');
end
E = cell(n+1,1); 
V = cell(n+1,1);
for j = 1:n+1 
    E{j} = M.basis(X{j}); 
    V{j} = zeros(size(X{j})); 
end
V{1} = v0; 
V{n+1} = vn;
for j = [1,n+1]
    c = coordinates(M,X{j},E{j},V{j}); 
    reconstructed = 0*V{j};
    for k = 1:d 
        reconstructed = reconstructed+c(k)*E{j}{k}; 
    end
    if norm(V{j}-reconstructed,'fro')>1e-9*(1+norm(V{j},'fro'))
        error('An endpoint velocity is not tangent.');
    end
end
K = spalloc(d*(n-1),d*(n-1),3*d*d*(n-1)); 
rhs = zeros(d*(n-1),1);
% Only three blocks per interior row are assembled. The off-diagonal blocks
% use the exact logarithm differential, not an identity or parallel transport.
dt = diff(t);
for j = 2:n
    row = (j-2)*d+(1:d);
    r = 3*M.log(X{j},X{j+1})/dt(j)^2-3*M.log(X{j},X{j-1})/dt(j-1)^2;
    rhs(row) = coordinates(M,X{j},E{j},r);
    K(row,row) = 2*(1/dt(j-1)+1/dt(j))*eye(d);
    for side = 1:2
        if side == 1 
            k = j-1; 
            step = dt(j-1); 
        else
            k = j+1; 
            step = dt(j); 
        end
        D = zeros(d);
        for l = 1:d
            D(:,l) = coordinates(M,X{j},E{j},M.dlog(X{j},X{k},E{k}{l}))/step;
        end
        if k == 1 || k == n+1
            % Endpoint velocities are prescribed, so move these terms to RHS.
            rhs(row) = rhs(row)-D*coordinates(M,X{k},E{k},V{k});
        else
            col = (k-2)*d+(1:d); 
            K(row,col) = D;
        end
    end
end
factor = struct;
if n > 1
    % Sparse LU: P*K*Q=L*U. Reuse these factors in every fixed-point step.
    [factor.L,factor.U,factor.P,factor.Q] = lu(K);
    c = factor.Q*(factor.U\(factor.L\(factor.P*rhs)));
    for j = 2:n
        row = (j-2)*d+(1:d);
        for l = 1:d 
            V{j} = V{j}+c(row(l))*E{j}{l}; 
        end
    end
    residual = norm(K*c-rhs,inf)/(1+norm(rhs,inf));
    % This model residual is NOT the actual covariant-acceleration jump.
else
    residual = 0;
    c = zeros(0,1);
end
S = struct('M',M,'t',t,'X',{X},'V',{V},'basis',{E}, ...
    'matrix',K,'rhs',rhs,'relative_residual',residual, ...
    'c',c,'factor',factor,'initializer_relative_residual',residual);
end

function c = coordinates(M,p,E,v)
c = zeros(numel(E),1);
for l = 1:numel(E) 
    c(l) = M.inner(p,E{l},v); 
end
end
