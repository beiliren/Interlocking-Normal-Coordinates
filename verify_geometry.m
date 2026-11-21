function diagnostics = verify_geometry()
%VERIFY_GEOMETRY Checks inverse maps, jets, metrics and a cubic reproduction.
rng(1729); 
diagnostics = struct;
for name = {'sphere','spd2'}
    M = manifold_geometry(name{1}); 
    width = numel(M.initial)/4;
    p = reshape(M.initial(1:width),M.shape); 
    E = M.basis(p);
    u = 0*p; 
    v = 0*p; 
    a = 0*p;
    for k = 1:M.dim
        u = u+.17*randn*E{k}; 
        v = v+.2*randn*E{k}; 
        a = a+.1*randn*E{k};
    end
    q = M.exp(p,u); 
    [~,dv,~] = M.jet(p,u,v,0*u);
    inv_error = tnorm(M,p,M.dlog(p,q,dv)-v);
    assert(inv_error<3e-11,'Inverse differential test failed.');
    assert(tnorm(M,p,M.log(p,q)-u)<3e-11,'Log/exp inverse test failed.');
    G = zeros(M.dim);
    for k = 1:M.dim
        for l = 1:M.dim 
            G(k,l) = M.inner(p,E{k},E{l}); 
        end
    end
    assert(norm(G-eye(M.dim),'fro')<3e-11,'Basis is not orthonormal.');
    [x,vel,acc] = M.jet(p,u,v,a); 
    step = 2e-5;
    [~,vm,~] = M.jet(p,u-step*v+.5*step^2*a,v-step*a,a);
    [~,vp,~] = M.jet(p,u+step*v+.5*step^2*a,v+step*a,a);
    fd = (vp-vm)/(2*step);
    if strcmp(M.name,'sphere') 
        fd = fd+(vel'*vel)*x; 
    else
        fd = fd-vel*(x\vel); 
    end
    jet_error = tnorm(M,x,fd-acc);
    assert(jet_error<2e-8,'Analytic acceleration test failed.');
    transported = M.pt(p,q,v);
    assert(abs(M.inner(p,v,v)-M.inner(q,transported,transported))<3e-11);
    % A scalar cubic along one geodesic is an exact Riemannian cubic.
    t = [0;.13;.4;.9;1.4;2]; 
    X = cell(size(t)); 
    truth = X; 
    direction = .4*E{1}+.15*E{end};
    for j = 1:numel(t)
        f = .2*t(j)+.07*t(j)^2-.015*t(j)^3;
        df = .2+.14*t(j)-.045*t(j)^2;
        ddf = .14-.09*t(j);
        [X{j},truth{j},~] = M.jet(p,f*direction,df*direction,ddf*direction);
    end
    S = interlocking_spline(M,t,X,truth{1},truth{end}); 
    reproduction = 0;
    for j = 1:numel(t) 
        reproduction = max(reproduction,tnorm(M,X{j},S.V{j}-truth{j})); 
    end
    assert(reproduction<2e-10,'Cubic/geodesic reproduction failed.');
    assert(max(acceleration_jumps(S))<2e-9);
    % A perturbed velocity vector tests the FULL nonlinear identity away
    % from the special reference solution, where cancellations are small.
    P = set_spline_coefficients(S,S.c+.15*randn(size(S.c)));
    r = actual_jump_coordinates(P);
    N = nonlinear_correction(P);
    splitting_error = norm(r+2*(P.matrix*P.c-P.rhs+N),inf);
    assert(splitting_error<2e-8,'Full nonlinear residual identity failed.');
    direct = solve_velocity_system(P,P.rhs-N);
    residual = P.c+solve_velocity_system(P,r/2);
    update_error = norm(direct-residual,inf);
    assert(update_error<2e-10,'Incorrect fixed-point update sign or factor.');
    diagnostics.(M.name) = struct('inverse_error',inv_error,'jet_error',jet_error, ...
        'reproduction_error',reproduction,'splitting_error',splitting_error, ...
        'update_error',update_error);
end
disp(diagnostics);
end

function n = tnorm(M,p,v)
n = sqrt(max(0,M.inner(p,v,v)));
end
