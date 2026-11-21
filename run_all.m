function run_all(outdir,reference_mode)
%RUN_ALL Piecewise-IVP benchmarks with 4 NEW knots after the identity.
% addpath('matlab'); run_all
% run_all('my_results','saved') replays supplied high-accuracy reference jets.
% All curves share positions and clamped PHYSICAL endpoint velocities.
% Energy uses sampled positions and scaled logarithmic differences only.
root = fileparts(fileparts(mfilename('fullpath')));
cfg = benchmark_config();
if nargin < 1
    outdir = fullfile(root,'computed');
end
if nargin < 2
    reference_mode = 'integrate';
end
for folder = {outdir,fullfile(outdir,'data'),fullfile(outdir,'figures')}
    if ~exist(folder{1},'dir')
        mkdir(folder{1});
    end
end
datadir = fullfile(outdir,'data');
figdir = fullfile(outdir,'figures');
diary(fullfile(outdir,'run_log.txt'));
cleanup = onCleanup(@()diary('off'));
fprintf('MATLAB %s; piecewise-IVP reference; mode %s\n',version,reference_mode);
checks = verify_geometry();
save(fullfile(datadir,'geometry_checks.mat'),'checks');
for names = {'sphere','spd2'}
    M = manifold_geometry(names{1});
    c = cfg.(M.name);
    if strcmp(reference_mode,'saved')
        D = load(fullfile(root,'data',[M.name,'_reference_dense.mat']));
        R = struct('M',M,'H',D.H,'n',numel(D.knot_t)-1, ...
            'T',D.t(end),'mode','saved','dense',D,'t',D.knot_t, ...
            'X',{D.knot_X},'V',{D.knot_V},'A',{D.knot_A}, ...
            'jerk_start',{D.jerk_start},'jerk_end',{D.jerk_end});
    elseif strcmp(reference_mode,'integrate')
        R = reference_spline(M,c,max(cfg.knot_counts));
        D = sample_reference(R,cfg.reference_samples_per_segment);
    else
        error('Reference mode must be integrate or saved.');
    end
    save(fullfile(datadir,[M.name,'_reference_dense.mat']),'-struct','D');
    record_reference_junctions(R,datadir);
    jerk = c.jerk_coefficients;
    rows = cell(size(jerk,1),2+size(jerk,2));
    for j = 1:size(jerk,1)
        rows(j,:) = [{j,(j-1)*R.H},num2cell(jerk(j,:))];
    end
    jerknames = {'segment','time'};
    for j = 1:size(jerk,2)
        jerknames{end+1} = sprintf('coefficient%d',j);
    end
    write_csv(fullfile(datadir,[M.name,'_prescribed_jerks.csv']),jerknames,rows);
    case_rows = {};
    energy_rows = {};
    separation_rows = {};
    histories = struct;
    for N = cfg.knot_counts
        C = make_curves(R,N,1,cfg);
        [errors,knot] = maximum_errors(C,cfg);
        [data,E] = case_output(C,cfg,datadir);
        histories.(sprintf('N%d',N)) = C.fp.history;
        [~,J] = actual_jump_coordinates(C.S);
        case_rows(end+1,:) = {N,C.T,C.T/C.n,C.fp.iterations, ... 
            errors(1,1),errors(1,2),errors(1,3),knot,J};
        for i = 1:height(E)
            energy_rows(end+1,:) = table2cell(E(i,:));
        end
        for k = 3:numel(C.methods)
            method = C.methods{k};
            distance = 0;
            % This is unscaled intrinsic separation, not an enlarged overlay.
            for j = 1:numel(data.t)
                x = data.interlocking_points{j};
                y = data.([method,'_points']);
                distance = max(distance,tnorm(M,x,M.log(x,y{j})));
            end
            separation_rows(end+1,:) = {N,method,distance};
        end
        fprintf('%s: %d new knots (%d total); %d corrections; actual jump %.3e\n', ...
            M.name,N,N+1,C.fp.iterations,J);
    end
    case_names = {'N','T','mesh_step','iterations','position','velocity','acceleration','knot_velocity','actual_jump'};
    cases = cell2table(case_rows,'VariableNames',case_names);
    write_csv(fullfile(datadir,[M.name,'_cases.csv']),case_names,case_rows);
    energy_names = {'N','M','delta','method','energy'};
    energies = cell2table(energy_rows,'VariableNames',energy_names);
    write_csv(fullfile(datadir,[M.name,'_energy_resolution.csv']),energy_names,energy_rows);
    write_csv(fullfile(datadir,[M.name,'_separation.csv']), ...
        {'N','method','maximum_distance'},separation_rows);

    % The number-of-knots tests change duration. A SEPARATE mesh refinement
    % of the fixed ten-segment reference is needed to test Theorem 3.
    rows = {};
    fprows = {};
    refhist = struct;
    for r = cfg.refinement_subdivisions
        C = make_curves(R,cfg.refinement_reference_knots,r,cfg);
        [errors,knot] = maximum_errors(C,cfg);
        [~,J] = actual_jump_coordinates(C.S);
        [~,J0] = actual_jump_coordinates(C.S0);
        [~,J1] = actual_jump_coordinates(C.Sfirst);
        q = NaN;
        if any(cfg.local_factor_subdivisions==r)
            q = fixed_point_local_factor(C.S0);
        end
        identity = actual_jump_coordinates(C.S0)+ ...
            2*(C.S0.matrix*C.S0.c-C.S0.rhs+nonlinear_correction(C.S0));
        assert(norm(identity,inf)<2e-8,'Incorrect nonlinear splitting.');
        for k = 2:numel(C.methods)
            key = C.methods{k};
            kv = NaN;
            jump = 0;
            if strcmp(key,'interlocking')
                kv = knot;
                jump = J;
            end
            rows(end+1,:) = {C.n,R.H/r,key,errors(k-1,1),errors(k-1,2), ... 
                errors(k-1,3),kv,jump};
        end
        fprows(end+1,:) = {r,C.n,R.H/r,C.fp.iterations,J0,J1,J,q,q/(R.H/r)^2}; 
        refhist.(sprintf('r%d',r)) = C.fp.history;
        save_knots(C,fullfile(datadir,sprintf('%s_refinement_r%d.mat',M.name,r)));
    end
    refinement = cell2table(rows,'VariableNames', ...
        {'N','h','method','position','velocity','acceleration','knot_velocity','actual_jump'});
    refinement = add_orders(refinement);
    write_csv(fullfile(datadir,[M.name,'_refinement.csv']), ...
        refinement.Properties.VariableNames,table2cell(refinement));
    fpnames = {'subdivisions','N','h','iterations','initial_jump','first_jump','final_jump','local_factor','local_factor_over_h2'};
    fixedpoint = cell2table(fprows,'VariableNames',fpnames);
    write_csv(fullfile(datadir,[M.name,'_fixed_point.csv']),fpnames,fprows);
    save(fullfile(datadir,[M.name,'_histories.mat']),'histories','refhist');
    plot_benchmark(M.name,datadir,cases,energies,refinement,fixedpoint,cfg,figdir);
end
fprintf('Finished: exact C2 concerns the limit; finite outputs meet actual-jump tolerance.\n');
end

function D = sample_reference(R,L)
% The configured grids divide this uniform, piecewise-IVP dense reference.
t = linspace(0,R.T,R.n*L+1)';
X = cell(size(t));
V = X;
A = X;
B = X;
for k = 1:numel(t)
    [X{k},V{k},A{k},B{k}] = reference_evaluate(R,t(k));
end
D = struct('t',t,'X',{X},'V',{V},'A',{A},'B',{B},'H',R.H, ...
    'knot_t',R.t,'knot_X',{R.X},'knot_V',{R.V},'knot_A',{R.A}, ...
    'jerk_start',{R.jerk_start},'jerk_end',{R.jerk_end});
end

function record_reference_junctions(R,datadir)
M = R.M;
rows = {};
for j = 2:R.n
    [xl,vl,al,bl] = reference_evaluate(R,R.t(j),'left');
    [xr,vr,ar,br] = reference_evaluate(R,R.t(j),'right');
    rows(end+1,:) = {j-1,R.t(j),tnorm(M,xl,M.log(xl,xr)), ... 
        tnorm(M,xr,M.pt(xl,xr,vl)-vr),tnorm(M,xr,M.pt(xl,xr,al)-ar), ...
        tnorm(M,xr,M.pt(xl,xr,bl)-br)};
end
write_csv(fullfile(datadir,[M.name,'_reference_junctions.csv']), ...
    {'knot','time','position_jump','velocity_jump','acceleration_jump','jerk_jump'},rows);
end

function C = make_curves(R,N,r,cfg)
M = R.M;
n = N*r;
t = linspace(0,N*R.H,n+1)';
X = cell(n+1,1);
V = X;
A = X;
for j = 1:n+1
    if mod(j-1,r) == 0
        k = (j-1)/r+1;
        X{j} = R.X{k};
        V{j} = R.V{k};
        A{j} = R.A{k};
    else
        [X{j},V{j},A{j}] = reference_evaluate(R,t(j));
    end
end
S0 = interlocking_spline(M,t,X,V{1},V{end});
[S,fp,first] = refine_fixed_point(S0,cfg.fixed_point_tolerance,cfg.fixed_point_max_iterations);
if ~fp.converged
    error('Fixed-point failure: %s N=%d r=%d (%s).',M.name,N,r,fp.reason);
end
Y = cell(size(X));
for j = 1:n+1
    Y{j} = M.log(X{1},X{j});
end
chart = euclidean_clamped(t,Y,V{1},M.dlog(X{1},X{end},V{end}));
methods = {'reference','interlocking','single_chart'};euclid=[];
if strcmp(M.name,'sphere')
    euclid = euclidean_clamped(t,X,V{1},V{end});
    methods{end+1} = 'radial_projection';
end
C = struct('M',M,'R',R,'Nreference',N,'subdivisions',r,'n',n,'T',N*R.H, ...
    't',t,'X',{X},'truth',{V},'reference_A',{A},'S0',S0,'S',S,'Sfirst',first, ...
    'fp',fp,'chart',chart,'euclid',euclid,'methods',{methods}, ...
    'P',polynomial_cache(S),'P0',polynomial_cache(S0));
end

function P = polynomial_cache(S)
M = S.M;
n = numel(S.t)-1;
U = cell(n,1);
W = U;
Q = U;
Qi = U;
for i = 1:n
    U{i} = M.log(S.X{i},S.X{i+1});
    W{i} = M.dlog(S.X{i},S.X{i+1},S.V{i+1});
    if strcmp(M.name,'spd2')
        [E,D] = eig(S.X{i});
        Q{i} = E*diag(sqrt(diag(D)))*E';
        Qi{i} = E*diag(1./sqrt(diag(D)))*E';
    end
end
P = struct('U',{U},'W',{W},'V',{S.V},'root',{Q},'inverse_root',{Qi});
end

function [z,dz,ddz] = polynomial(C,key,i,s)
P = C.P;
if strcmp(key,'initializer')
    P = C.P0;
end
dt = C.t(i+1)-C.t(i);
u = P.U{i};
v = P.V{i};
w = P.W{i};
z = (3*s^2-2*s^3)*u+dt*(s-2*s^2+s^3)*v+dt*(s^3-s^2)*w;
dz = (6*s-6*s^2)*u/dt+(1-4*s+3*s^2)*v+(3*s^2-2*s)*w;
ddz = (6-12*s)*u/dt^2+(-4+6*s)*v/dt+(6*s-2)*w/dt;
end

function [x,v,a] = curve_jet(C,key,t,interval)
if nargin < 4
    interval = [];
end
switch key
    case 'reference'
        [x,v,a] = reference_evaluate(C.R,t);
    case {'interlocking','initializer'}
        i = interval;
        if isempty(i)
            i = min(C.n,find(C.t<=t,1,'last'));
        end
        s = (t-C.t(i))/(C.t(i+1)-C.t(i));
        [z,dz,ddz] = polynomial(C,key,i,s);
        [x,v,a] = C.M.jet(C.X{i},z,dz,ddz);
    case 'single_chart'
        [z,dz,ddz] = euclidean_eval(C.chart,t);[x,v,a]=C.M.jet(C.X{1},z,dz,ddz);
    case 'radial_projection'
        [y,dy,ddy] = euclidean_eval(C.euclid,t);rho=norm(y);
        if rho < 1e-8
            error('Ambient curve approaches projection singularity.');
        end
        dr = y'*dy/rho;
        ddr = (dy'*dy+y'*ddy-dr^2)/rho;
        x = y/rho;
        v = dy/rho-y*dr/rho^2;
        ddx = ddy/rho-2*dy*dr/rho^2-y*ddr/rho^2+2*y*dr^2/rho^3;
        a = ddx+(v'*v)*x;
end
end

function x = curve_point(C,key,t)
% Position-only evaluation keeps energy independent of acceleration jets.
switch key
    case 'reference'
        x = reference_evaluate(C.R,t);
    case 'interlocking'
        i = min(C.n,find(C.t<=t,1,'last'));
        s = (t-C.t(i))/(C.t(i+1)-C.t(i));
        z = polynomial(C,key,i,s);
        if strcmp(C.M.name,'spd2')
            P = C.P;
            Z = P.inverse_root{i}*z*P.inverse_root{i};
            x = P.root{i}*expm((Z+Z')/2)*P.root{i};
            x = (x+x')/2;
        else
            x = C.M.exp(C.X{i},z);
        end
    case 'single_chart'
        z = euclidean_eval(C.chart,t);
        x = C.M.exp(C.X{1},z);
    case 'radial_projection'
        y = euclidean_eval(C.euclid,t);
        if norm(y) < 1e-8
            error('Projection singularity.');
        end
        x = y/norm(y);
end
end

function [err,knot] = maximum_errors(C,cfg)
err = zeros(numel(C.methods)-1,3);
knot = 0;
M = C.M;
for j = 1:C.n+1
    knot = max(knot,tnorm(M,C.X{j},C.S.V{j}-C.truth{j}));
end
for i = 1:C.n
    for s = linspace(0,1,cfg.error_points_per_interval)
        t = C.t(i)+s*(C.t(i+1)-C.t(i));
        [xr,vr,ar] = reference_evaluate(C.R,t);
        for k = 2:numel(C.methods)
            key = C.methods{k};
            [x,v,a] = curve_jet(C,key,t,i);
            value = [tnorm(M,xr,M.log(xr,x)),tnorm(M,xr,M.pt(x,xr,v)-vr), ...
                tnorm(M,xr,M.pt(x,xr,a)-ar)];
            err(k-1,:) = max(err(k-1,:),value);
        end
    end
end
end

function [data,E] = case_output(C,cfg,datadir)
M = C.M;
L = max(cfg.energy_subintervals);
tt = linspace(0,C.T,L+1)';
samples = struct('t',tt);
erows = {};
for k = 1:numel(C.methods)
    key = C.methods{k};
    points = cell(size(tt));
    for j = 1:numel(tt)
        points{j} = curve_point(C,key,tt(j));
    end
    samples.(key) = points;
    for count = cfg.energy_subintervals
        assert(mod(L,count)==0);
        selected = points(1:L/count:end);
        erows(end+1,:) = {C.Nreference,count,C.T/count,key,discrete_cubic_energy(M,selected,C.T)};
    end
end
E = cell2table(erows,'VariableNames',{'N','M','delta','method','energy'});
% Report the measured ordering honestly; the solver never optimizes energy.
for count = cfg.energy_subintervals
    subset = E(E.M==count,:);
    own = subset.energy(strcmp(subset.method,'interlocking'));
    other = subset.energy(~strcmp(subset.method,'reference') & ~strcmp(subset.method,'interlocking'));
    assert(all(own<other),'Measured example energy ordering changed.');
end
save(fullfile(datadir,sprintf('%s_N%d_energy_samples.mat',M.name,C.Nreference)),'-struct','samples');
tp = linspace(0,C.T,cfg.plot_subintervals+1)';
data = struct('t',tp,'N',C.Nreference,'knot_points',{C.X},'knot_times',C.t);
for k = 1:numel(C.methods)
    key = C.methods{k};
    X = cell(size(tp));
    V = X;
    A = X;
    speed = zeros(size(tp));
    for j = 1:numel(tp)
        [X{j},V{j},A{j}] = curve_jet(C,key,tp(j));
        speed(j)=tnorm(M,X{j},V{j});
    end
    data.([key,'_points']) = X;
    data.([key,'_velocity']) = V;
    data.([key,'_acceleration']) = A;
    data.([key,'_speed']) = speed;
end
save(fullfile(datadir,sprintf('%s_N%d_curves.mat',M.name,C.Nreference)),'-struct','data');
save_knots(C,fullfile(datadir,sprintf('%s_N%d_data.mat',M.name,C.Nreference)));
write_csv(fullfile(datadir,sprintf('%s_N%d_iterations.csv',M.name,C.Nreference)), ...
    {'iteration','actual_jump','coefficient_increment','ratio'},num2cell(C.fp.history));
if strcmp(M.name,'spd2')
    entries = tp;
    names = {'time'};
    for k = 1:numel(C.methods)
        key = C.methods{k};
        points = data.([key,'_points']);
        values = zeros(numel(tp),3);
        for j = 1:numel(tp)
            x = points{j};
            values(j,:) = [x(1,1),x(1,2),x(2,2)];
        end
        entries = [entries,values];
        names = [names,{[key,'_x11'],[key,'_x12'],[key,'_x22']}];
    end
    write_csv(fullfile(datadir,sprintf('spd2_N%d_components.csv',C.Nreference)),names,num2cell(entries));
end
for k = 1:numel(C.methods)
    for j = 1:C.n+1
        [x,v,~] = curve_jet(C,C.methods{k},C.t(j));
        assert(tnorm(M,C.X{j},M.log(C.X{j},x))<2e-9);
        if j == 1 || j == C.n+1
            assert(tnorm(M,C.X{j},M.pt(x,C.X{j},v)-C.truth{j})<2e-8);
        end
    end
end
end

function save_knots(C,path)
t = C.t;
X = C.X;
truth = C.truth;
reference_A = C.reference_A;
V = C.S.V;
V_linear = C.S0.V;
V_first = C.Sfirst.V;
jump = acceleration_jumps(C.S);
jump_linear = acceleration_jumps(C.S0);
jump_first = acceleration_jumps(C.Sfirst);
history = C.fp.history;
save(path,'t','X','truth','reference_A','V','V_linear','V_first','jump','jump_linear','jump_first','history');
end

function tab = add_orders(tab)
for field = {'position','velocity','acceleration'}
    value = tab.(field{1});
    order = NaN(height(tab),1);
    for keys = unique(tab.method,'stable')'
        indices = find(strcmp(tab.method,keys{1}));
        for j = 2:numel(indices)
            a = indices(j-1);
            b = indices(j);
            order(b) = log(value(a)/value(b))/log(tab.h(a)/tab.h(b));
        end
    end
    tab.([field{1},'_order']) = order;
end
end

function value = tnorm(M,p,v),value = sqrt(max(0,M.inner(p,v,v)));end
function write_csv(path,names,rows)
f = fopen(path,'w');
if f < 0
    error('Cannot write output.');
end
c = onCleanup(@()fclose(f));
fprintf(f,'%s',names{1});
for j = 2:numel(names)
    fprintf(f,',%s',names{j});
end
fprintf(f,'\n');
for i = 1:size(rows,1)
    for j = 1:size(rows,2)
        if j > 1
            fprintf(f,',');
        end
        v = rows{i,j};
        if ischar(v)
            fprintf(f,'%s',v);
        else
            fprintf(f,'%.17g',v);
        end
    end
    fprintf(f,'\n');
end
end
