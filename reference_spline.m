function R = reference_spline(M,c,n)
%REFERENCE_SPLINE Build a genuine piecewise Riemannian cubic by IVPs.
% Start at c.initial_position, c.initial_velocity, c.initial_acceleration,
% c.initial_jerk. Integrate each IVP for c.segment_duration. Carry the END
% position, velocity and COVARIANT acceleration into the next IVP, and reset
% only its covariant jerk. This gives a C2 reference, generally with jumps
% in the third covariant derivative. It is not one long cubic IVP.
%
% The benchmarks use PRESCRIBED jerk tables c.jerk_coefficients. Thus every
% jerk table is fixed input shared by every compared interpolant.
% If a jerk table is absent, an optional feedback rule is
%   b_i = -damping*a_i/H + tracking*(w_i-v_i)/H^2.
% The target w_i is specified by c.target_coefficients in an orthonormal
% frame. These are benchmark inputs, independent of every interpolant.
% Frame coefficients are converted to actual tangent vectors at each knot.
H = c.segment_duration;
if n < 1 || n ~= floor(n) || H <= 0
    error('Invalid reference mesh.');
end
if isfield(c,'jerk_coefficients')
    if size(c.jerk_coefficients,1) < n
        error('Missing prescribed jerks.');
    end
elseif ~isfield(c,'target_coefficients') || size(c.target_coefficients,1) < n
    error('Supply a jerk table or target coefficients.');
end
R = struct('M',M,'H',H,'n',n,'T',n*H,'mode','integrate');
R.t = (0:n)'*H;
R.X = cell(n+1,1);
R.V = R.X;
R.A = R.X;
R.jerk_start = cell(n,1);
R.jerk_end = cell(n,1);
R.segments = cell(n,1);
R.target = cell(n,1);
R.frame = cell(n,1);
x = c.initial_position;
v = c.initial_velocity;
a = c.initial_acceleration;
E = c.initial_frame;
opts = odeset('RelTol',2e-13,'AbsTol',1e-14,'MaxStep',H/32);
for i = 1:n
    if strcmp(M.name,'spd2')
        % The positive square root gives a canonical frame. Its first two
        % vectors span the curved determinant-one SPD(2) submanifold.
        B = M.basis(x);
        E = {(B{1}-B{2})/sqrt(2),B{3},(B{1}+B{2})/sqrt(2)};
    end
    w = zeros(size(x));
    if isfield(c,'target_coefficients')
        for j = 1:numel(E)
            w = w+c.target_coefficients(i,j)*E{j};
        end
    end
    if i == 1
        b = c.initial_jerk;
    elseif isfield(c,'jerk_coefficients')
        b = zeros(size(x));
        for j = 1:numel(E)
            b = b+c.jerk_coefficients(i,j)*E{j};
        end
    else
        b = -c.damping*a/H+c.tracking*(w-v)/H^2;
    end
    R.X{i} = x;
    R.V{i} = v;
    R.A{i} = a;
    R.jerk_start{i} = b;
    R.target{i} = w;
    R.frame{i} = E;
    initial = [x(:);v(:);a(:);b(:)];
    sol = ode45(M.ode,[0,H],initial,opts);
    [xn,vn,an,bn] = reference_state(M,deval(sol,H));
    R.segments{i} = sol;
    R.jerk_end{i} = bn;
    if strcmp(M.name,'sphere')
        % Transport this frame between consecutive knots. Gram--Schmidt
        % removes floating-point drift, without changing the curve data.
        E = {M.pt(x,xn,E{1}),M.pt(x,xn,E{2})};
        E{1} = E{1}/norm(E{1});
        E{2} = E{2}-E{1}*(E{1}'*E{2});
        E{2} = E{2}/norm(E{2});
    end
    x = xn;
    v = vn;
    a = an;  % Carry all three end jets; replace ONLY jerk next time.
end
R.X{n+1} = x;
R.V{n+1} = v;
R.A{n+1} = a;
end
