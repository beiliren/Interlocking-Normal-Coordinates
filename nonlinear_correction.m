function N = nonlinear_correction(S)
%NONLINEAR_CORRECTION Full correction for the ORIGINAL block matrix K.
% N_j = E_j^{-1}[(L_{j-1}-D_{j,j-1})v_{j-1}/Delta_{j-1}
%                + B_{j-1}[D_{j-1,j}v_j,D_{j-1,j}v_j]/2].
% The small linear geometric term is essential. Omitting it while retaining
% the original K would solve a different acceleration-matching problem.
% M.jet with ddz=0 evaluates L*dz and B[dz,dz] without numerical derivatives.
M = S.M;
n = numel(S.t)-1;
d = M.dim;
N = zeros(d*(n-1),1);
for j = 2:n
    p = S.X{j-1}; q = S.X{j};
    dt = S.t(j)-S.t(j-1);
    u = M.log(p,q);
    w = M.dlog(p,q,S.V{j});
    [~,Lv,~] = M.jet(p,u,S.V{j-1},0*u);
    [~,~,Bww] = M.jet(p,u,w,0*u);
    value = (Lv-M.dlog(q,p,S.V{j-1}))/dt+Bww/2;
    row = (j-2)*d+(1:d);
    for l = 1:d 
        N(row(l)) = M.inner(q,S.basis{j}{l},value); 
    end
end
end
