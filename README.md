This code is to reproduce the numerical results in the paper,

"Interlocking Normal Coordinates for Approximating Riemannian Cubic Splines"

by Lyle Noakes and Erchuan Zhang.

## Run the MATLAB code

From the extracted package root:

```matlab
addpath(fullfile(pwd,'matlab'))
replot_saved_results                     % redraw supplied arrays
report = check_saved_data;               % recompute saved velocities, jumps and energies
run_saved_benchmarks                     % complete replay using supplied reference jets
run_all                                 % new piecewise-IVP integrations and full experiment
```

Base MATLAB R2019b or newer is sufficient, without extra toolboxes. `run_all` writes `computed/`. The default maximum actual acceleration-jump tolerance is 1e-9, with at most 200 fixed-point corrections. Failure is reported explicitly. Use `run_all` for new initial conditions, jerk tables or sampling grids.
