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

sim.xi       = 0;
sim.m0       = 30e-3;
sim.e        = 10e-3;
sim.fric_mod = 'ambrosio';

% Initial conditions
sim.y_0(1)  = 0.002796;
sim.y_0(6)  = 18.9*2*pi;
sim.y_0(9)  = -3.72e-6;
sim.y_0(13) = -3.626e-6;

% Solve
sim.o45_reltol = 1e-9;
sim.o45_abstol = 1e-9;
sim.o15_reltol = 1e-9;
sim.o15_abstol = 1e-9;

sim.solve([0 0.1])

% Calculate forces etc.
sim.postprocess();

sim.export('basic')

% Plot stuff
pt = Plottools();

pt.debugplot('./exp.txt')
