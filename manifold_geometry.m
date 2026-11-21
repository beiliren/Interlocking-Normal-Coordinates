function M = manifold_geometry(kind)
%MANIFOLD_GEOMETRY Exact maps and analytic two-jets for S^2 and SPD(2).
% All tangent bases are orthonormal for the stated Riemannian metric.
% No Optimization Toolbox, Symbolic Toolbox or external manifold library.
% M.pt(q,p,v) transports FROM q TO p. M.dlog(p,q,v) maps a vector at q
% into the tangent space at p; these two operations are different.
% M.jet(p,z,dz,ddz) returns position, physical velocity and COVARIANT
% acceleration of exp_p(z(t)), with p fixed. M.initial is a geometry-test
% fixture; benchmark_config.m supplies the actual reference initial states.
switch lower(kind)
    case 'sphere'
        M.name = 'sphere'; 
        M.dim = 2; 
        M.shape = [3,1];
        M.inner = @(p,u,v) u(:)'*v(:);
        M.basis = @sphere_basis;
        M.exp = @sphere_exp; 
        M.log = @sphere_log;
        M.pt = @sphere_pt; 
        M.dlog = @sphere_dlog;
        M.jet = @sphere_jet; 
        M.ode = @sphere_ode;
        M.initial = [1;0;0; 0;.7;.2; 0;.0075;-.0125; 0;-.005;.015];
    case 'spd2'
        M.name = 'spd2'; 
        M.dim = 3; 
        M.shape = [2,2];
        M.inner = @(p,u,v) trace((p\u)*(p\v));
        M.basis = @spd_basis;
        M.exp = @spd_exp; 
        M.log = @spd_log;
        M.pt = @spd_pt; 
        M.dlog = @spd_dlog;
        M.jet = @spd_jet; 
        M.ode = @spd_ode;
        X = eye(2); 
        V = [.35,.2;.2,-.1];
        A = [.1,-.08;-.08,.15]; 
        B = [-.06,.04;.04,.03];
        M.initial = [X(:);V(:);A(:);B(:)];
    otherwise
        error('Unknown manifold. Use sphere or spd2.');
end
end

function E = sphere_basis(p)
Q = null(p(:)'); 
E = cell(1,size(Q,2));
for k = 1:numel(E) 
    E{k} = Q(:,k); 
end
end

function q = sphere_exp(p,u)
r = norm(u);
if r > 1e-10 
    s = sin(r)/r; 
else
    s = 1-r^2/6+r^4/120; 
end
q = cos(r)*p+s*u;
end

function u = sphere_log(p,q)
c = max(-1,min(1,p'*q)); 
u = q-c*p; 
sn = norm(u); 
r = atan2(sn,c);
if r > pi-1e-7 
    error('Logarithm is too close to the sphere cut locus.'); 
end
if sn > 1e-14 
    u = (r/sn)*u; 
end
end

function w = sphere_pt(q,p,v)
% Parallel transport FROM q TO p along the short great-circle arc.
if 1+q'*p < 1e-12 
    error('Parallel transport at antipodes is ambiguous.'); 
end
w = v-((v'*p)/(1+q'*p))*(q+p);
end

function w = sphere_dlog(p,q,v)
% Undo parallel transport, then rescale only the transverse component.
% The radial component has unit scale; the transverse scale is r/sin(r).
u = sphere_log(p,q); 
r = norm(u); 
w = sphere_pt(q,p,v);
if r < 1e-8 
    return; 
end
e = u/r; 
wpar = (e'*w)*e;
w = wpar+(r/sin(r))*(w-wpar);
end

function [x,v,a] = sphere_jet(p,z,dz,ddz)
% Analytic derivatives of exp_p(z(t)); no differencing across a junction.
% Differentiate sinc(sqrt(q)) as a function of q = z'*z. The series removes
% cancellation near q=0 and makes acceleration at the chart center stable.
q = z'*z;
if q < 1e-3
    s = 1-q/6+q^2/120-q^3/5040+q^4/362880;
    sq = -1/6+q/60-q^2/1680+q^3/90720;
    sqq = 1/60-q/840+q^2/30240-q^3/1995840;
else
    r = sqrt(q); 
    s = sin(r)/r;
    sq = (r*cos(r)-sin(r))/(2*r^3);
    sqq = ((3-r^2)*sin(r)-3*r*cos(r))/(4*r^5);
end
c = cos(sqrt(q)); 
cq = -s/2; 
cqq = -sq/2;
dq = 2*(z'*dz); 
ddq = 2*(dz'*dz+z'*ddz);
x = c*p+s*z;
v = cq*dq*p+sq*dq*z+s*dz;
ddx = (cqq*dq^2+cq*ddq)*p+(sqq*dq^2+sq*ddq)*z ...
    +2*sq*dq*dz+s*ddz;
a = ddx+(v'*v)*x;
% On the unit sphere, covariant acceleration is ddx + ||v||^2*x.
end

function dy = sphere_ode(~,y)
Y = reshape(y,3,4); 
x = Y(:,1); 
v = Y(:,2); 
a = Y(:,3); 
b = Y(:,4);
Ydot = [v, a-(v'*v)*x, b-(v'*a)*x, ...
    -(v'*v)*a+(a'*v)*v-(v'*b)*x];
dy = Ydot(:);
end

function X = sm(X)
X = (X+X')/2;
end

function S = spow(X,power)
[Q,D] = eig(sm(X)); 
d = diag(D);
if min(d) <= 0 
    error('A matrix left the positive definite cone.'); 
end
S = sm(Q*diag(d.^power)*Q');
end

function E = spd_basis(X)
% Congruence by X^(1/2) converts a Frobenius-orthonormal symmetric basis
% into an affine-invariant orthonormal basis at X.
S = spow(X,.5);
B = { [1,0;0,0], [0,0;0,1], [0,1;1,0]/sqrt(2) };
E = cell(1,3);
for k = 1:3 
    E{k} = sm(S*B{k}*S); 
end
end

function Y = spd_exp(X,U)
S = spow(X,.5); 
Si = spow(X,-.5);
Y = sm(S*expm(sm(Si*U*Si))*S);
end

function U = spd_log(X,Y)
S = spow(X,.5); 
Si = spow(X,-.5); 
[Q,D] = eig(sm(Si*Y*Si));
d = diag(D);
if min(d) <= 0 
    error('Nonpositive logarithm argument.'); 
end
U = sm(S*Q*diag(log(d))*Q'*S);
end

function W = spd_pt(Q,P,V)
% E*Q*E' = P; congruence by E is transport along the affine-invariant geodesic.
S = spow(Q,.5); 
Si = spow(Q,-.5);
E = S*spow(Si*P*Si,.5)*Si;
W = sm(E*V*E');
end

function W = spd_dlog(X,Y,V)
% Frechet derivative of the symmetric matrix logarithm, via divided
% differences in an eigenbasis. No assumption of commuting matrices is made.
S = spow(X,.5); 
Si = spow(X,-.5); 
[Q,D] = eig(sm(Si*Y*Si));
d = diag(D); 
F = zeros(2);
if min(d) <= 0 
    error('Nonpositive logarithm argument.'); 
end
for i = 1:2
    for j = 1:2
        gap = d(i)-d(j);
        if abs(gap)<1e-7*max(d(i),d(j))
            % Stable even series for (log(d_i)-log(d_j))/(d_i-d_j).
            mid = (d(i)+d(j))/2; 
            t = gap/(2*mid);
            F(i,j) = (1+t^2/3+t^4/5)/mid;
        else
            F(i,j) = (log(d(i))-log(d(j)))/gap;
        end
    end
end
W = sm(S*Q*(F.*(Q'*Si*V*Si*Q))*Q'*S);
end

function [X,V,A] = spd_jet(P,Z,dZ,ddZ)
S = spow(P,.5); 
Si = spow(P,-.5);
L = sm(Si*Z*Si); 
E = sm(Si*dZ*Si); 
F = sm(Si*ddZ*Si);
% A 3-by-3 block exponential includes the full noncommuting second
% derivative. The (1,3) block encodes one half of the second time derivative.
T = [L,E,F/2; zeros(2),L,E; zeros(2),zeros(2),L];
B = expm(T);
X = sm(S*B(1:2,1:2)*S); 
V = sm(S*B(1:2,3:4)*S);
ddX = sm(2*S*B(1:2,5:6)*S);
A = sm(ddX-V*(X\V));
% The affine-invariant connection gives A = X'' - X'*inv(X)*X'.
end

function R = spd_curvature(X,U,V,W)
Xi = inv(X);
R = sm(-.25*(U*Xi*V*Xi*W-V*Xi*U*Xi*W ...
    -W*Xi*U*Xi*V+W*Xi*V*Xi*U));
end

function dy = spd_ode(~,y)
% Reference cubic: nabla_t v=a, nabla_t a=b, nabla_t b=-R(a,v)v.
% Matrices are stacked columnwise; symmetrization removes ODE roundoff.
Y = reshape(y,4,4);
X = sm(reshape(Y(:,1),2,2)); 
V = sm(reshape(Y(:,2),2,2));
A = sm(reshape(Y(:,3),2,2)); 
B = sm(reshape(Y(:,4),2,2));
dX = V; 
dV = A+V*(X\V);
dA = B+.5*(V*(X\A)+A*(X\V));
dB = -spd_curvature(X,A,V,V)+.5*(V*(X\B)+B*(X\V));
dy = [dX(:);dV(:);dA(:);dB(:)];
end
