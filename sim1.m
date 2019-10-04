% Master thesis F2019
% Numerical integration of analytical solution
% Svend E. Andersen
clc; clearvars; %close all

% Initiate simulation
sim = Simulator;
sim.fric_mod = 'ambrosio';

% Enable magnet
load('magnet-force/rot_sim_exp315.mat', 'exp')
sim.set_magnet(11.9, 240, exp)

% Initial conditions
sim.y_0(6) = -17*2*pi;

% Solve
sim.solve([0 12.1])

% Calculate forces etc.
sim.postprocess();

%sim.export('mat')

% Plot stuff
pt = Plottools();

%pt.orbit(sim.r_OD(1, :), sim.r_OD(2, :), sim.clearance)
%pt.lateral(sim.time, sim.r_OD(1, :), sim.r_OD(2, :))
%pt.indent(sim.time, sim.d, sim.event_times)
%pt.indentf(sim.time, sim.d, sim.fn, sim.event_times)
%pt.states(sim.time, sim.solution(:, 1:2:end), {'\Gamma', '\beta', '\theta', ...
          %'x_{ih}', 'y_{ih}', 'x_{oh}', 'y_{ih}'})
