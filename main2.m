% Master thesis F2019
% Numerical integration of analytical solution
% Svend E. Andersen
%
% Substitution definition:
%   y = [gamma, gamma_d, beta, beta_d, theta, theta_d, x_ih, x_ih_d, y_ih,...
%        y_ih_d, x_mh, x_mh_d, y_mh, y_mh_d]^T
clc; clearvars; close all

% Initiate simulation
sim = Simulator;

sim.fric_mod = 'ambrosio';

% Enable magnet
load('exc_17hz_315deg.mat', 'exc')
sim.set_magnet(9.9, 245, exc)

% Initial conditions
sim.y_0(6) = -17*2*pi;

% Solve
%sim.o15_reltol = 1e-8; % Relative tolerance for the ode15s solver
%sim.o15_abstol = 1e-8; % Absolute tolerance for the ode15s solver

sim.solve([0 10.05])

% Calculate forces etc.
sim.postprocess();

%sim.export('basic')

% Plot stuff
pt = Plottools();

%pt.orbit(sim.r_OD(1, :), sim.r_OD(2, :), sim.clearance)

pt.lateral(sim.time, sim.r_OD(1, :), sim.r_OD(2, :))

pt.orbit(sim.r_OD(1, round((3/8)*end):end), sim.r_OD(2, round((3/8)*end):end), ...
         sim.clearance)

%pt.states(sim.time, sim.solution(:, 1:2:end), {'\Gamma', '\beta', '\theta', ...
          %'x_{ih}', 'y_{ih}', 'x_{oh}', 'y_{ih}'})
