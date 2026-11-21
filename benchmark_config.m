function cfg = benchmark_config()
%BENCHMARK_CONFIG Fixed, reproducible piecewise-IVP spline data.
% N is the number of NEW knots after the starting point: N+1 in total.
cfg.knot_counts = 5;%[4,6,8,10];
cfg.refinement_reference_knots = 5;%10;
cfg.refinement_subdivisions = [1,2,4,8,16];
cfg.local_factor_subdivisions = [1,2,4,8,16];
cfg.reference_samples_per_segment = 2048;%4096;
cfg.error_points_per_interval = 33;
cfg.plot_subintervals = 512;
cfg.energy_subintervals = 2048;%[2048,4096,8192];
cfg.figure_visible = 'on';
cfg.sphere_show_single_chart_speed = false;
cfg.fixed_point_tolerance = 1e-9;
cfg.fixed_point_max_iterations = 200;

% Prescribed jerks: each row is a fresh covariant jerk at a reference knot.
cfg.sphere.segment_duration = 1;
cfg.sphere.initial_position = [1;0;0];
cfg.sphere.initial_velocity = [0;1.8;0.12];
cfg.sphere.initial_acceleration = [0;0;0.3];
cfg.sphere.initial_frame = {[0;1;0],[0;0;1]};
% The spherical frame is transported along each successive short knot geodesic.
cfg.sphere.jerk_coefficients = [
    0,-0.5700;
    -0.6936,1.9304;
    -0.6459,-0.6209;
    -1.2812,0.4240;
    1.0791,-1.0745
];
cfg.sphere.initial_jerk = zeros(size(cfg.sphere.initial_position));
for j = 1:numel(cfg.sphere.initial_frame)
    cfg.sphere.initial_jerk = cfg.sphere.initial_jerk + ...
        cfg.sphere.jerk_coefficients(1,j)*cfg.sphere.initial_frame{j};
end

% Prescribed jerks: each row is a fresh covariant jerk at a reference knot.
cfg.spd2.segment_duration = 1;
cfg.spd2.initial_position = eye(2);
cfg.spd2.initial_velocity = [1.8/sqrt(2),0; 0,-1.8/sqrt(2)];
cfg.spd2.initial_acceleration = [0,0.4/sqrt(2); 0.4/sqrt(2),0];
A = [1,0;0,-1]/sqrt(2);
B = [0,1;1,0]/sqrt(2);
C = eye(2)/sqrt(2);
cfg.spd2.initial_frame = {A,B,C};
% At X, the SPD frame is X^(1/2)*{A,B,C}*X^(1/2), recomputed at each knot.
cfg.spd2.jerk_coefficients = [
    -2.2416,2.0884,0;
    1.8608,-8.1305,0;
    8.6887,-0.6722,0;
    8.2040,2.4290,0;
    4.5980,13.6021,0
];
cfg.spd2.initial_jerk = zeros(size(cfg.spd2.initial_position));
for j = 1:numel(cfg.spd2.initial_frame)
    cfg.spd2.initial_jerk = cfg.spd2.initial_jerk + ...
        cfg.spd2.jerk_coefficients(1,j)*cfg.spd2.initial_frame{j};
end

% The jerk table is part of the reference input. Interpolation never changes
% it or selects a jerk by comparing energies. New knots are generated only
% by IVPs; position, velocity and covariant acceleration are carried across.
end
