% Master thesis F2019
% Numerical integration of analytical solution
% Svend E. Andersen
%
% Substitution definition:
%   y = [gamma, gamma_d, beta, beta_d, theta, theta_d, x_ih, x_ih_d, y_ih,...
%        y_ih_d, x_mh, x_mh_d, y_mh, y_mh_d]^T
clc
clearvars
close all
format long

% Initiate system class with parameters 'xi', 'm0' and 'e'
s = Rotorsystem(0.006, 30e-3, 10e-3);

% Define the contact model
cmod = Nikravesh(s.r_s, s.r_r);

% Motor angular velocity [Hz]
Omega = 18.9*2*pi;

% % Integration % %
% Integration time span [s]
tspan = [0, 0.5];

% Initial conditions
ih_offset = -3.72e-6;

y_0       = zeros(14,1);

y_0(1)    = 0.002796;
y_0(9)    = ih_offset;
y_0(13)   = -3.626e-6;
y_0(6)    = Omega;

% Solver options
options_ode45 = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, ...
                       'Events', @(t,y) impact_detect(t, y, s, 1));
options_ode15 = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, ...
                       'Events', @(t,y) impact_detect(t, y, s, -1));

% Solution containers
t_total  = 0;   % time
s_total  = 0;   % contact state
y_total  = [];  % y
te_total = [];  % event times

loc_tst = tspan(1); % Integration time starting point

tic
while t_total(end) ~= tspan(2)

  gap = s.calc_gap(y_0);

  if gap < 0
    contact_state = 0;
    [t,y,te,ye,ie] =  ode45(@(t,y) dydt(t,y,s,cmod,contact_state), ...
                      [loc_tst,tspan(2)], y_0, options_ode45);
  else
    contact_state = 1;
    [t,y,te,ye,ie] = ode15s(@(t,y) dydt(t,y,s,cmod,contact_state), ...
                     [loc_tst,tspan(2)], y_0, options_ode15);
  end

  % Collect results
  t_total  = [ t_total(1:end-1)  ; t ];
  y_total  = [ y_total(1:end-1,:); y ];
  s_total  = [ s_total(1:end-1)  ; contact_state*ones(length(t),1) ];
  te_total = [te_total; te];

  % Assign new initial conditions
  loc_tst = t(end);
  y_0 = y(end,:);

  t_total(end)
end
toc

fprintf('%i perimeter crossings detected\n', length(te_total))

return

% Write results to file
fileID = fopen('test.txt','w');
fprintf(fileID,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',...
  [t_total, y_total]');
fclose(fileID);

% Get contact forces
for i = 1:length(t_total)
  y_i = y_total(i,:)';

  if s_total(i) == 0, state = 0; else, state = 1;  end

  [F_cxs(i), F_cys(i), deltas(i) ] = contact_force(y_i, s, cmod, state);
end

return
%animate(t_total,y_total,s_total,te_total,s,F_cxs,F_cys,deltas)

% Get rotor centre in the contact plane
pos1 = zeros(3, length(t_total));
for i = 1:length(t_total)
  T_gamma = s.T_gam(y_total(i,1));
  T_beta  = s.T_bet(y_total(i,3));
  T_theta = s.T_the(y_total(i,5));
  pos1(:,i) = T_gamma.' * (T_beta.'* [0;0;s.l_OC] );
end


% Plots %

% Force-identation plot
%figure(); hold on
%plot(deltas, sqrt(F_cxs.^2+F_cys.^2), 'k'); grid on; hold on
%xlabel('Indentation [m]')
%ylabel('Contact force [N]')

% Orbit plot
figure('units','normalized','outerposition',[0.5 0 0.6 1])
subplot(1,2,1)
plot(t_total, pos1(1,:)*1e3,'k'); grid on; hold on
plot(t_total, pos1(2,:)*1e3,'r');
xlabel('Time [s]')
ylabel('Deflection [mm]')
legend('x','y')

subplot(1,2,2)
plot(pos1(1,:)*1e3, pos1(2,:)*1e3,'b'); grid on; hold on
plot( s.cl*cos(linspace(0,2*pi))*1e3,...
     (s.cl*sin(linspace(0,2*pi))+ih_offset)*1e3,'r')
xlabel('Orbit [mm]')
axis equal


% State vector plot
figure('units','normalized','outerposition',[0.5 0 0.6 1])
j = 1;
for i = 1:2:size(y_total,2)
  subplot(size(y_total,2)/2, 1, j)
  plot(t_total,y_total(:,i));

  posi = get(gca, 'Position');
  posi(1) = 0.055;
  posi(3) = 0.9;
  set(gca, 'Position', posi)

  j = j + 1;
end

figure()
subplot(2,1,1)
plot(t_total, sign(delta_ds))
subplot(2,1,2)
plot(t_total, s_total)

figure()
plot(t_total, sign(delta_ds)); hold on
plot(t_total, s_total)

return

% Force plot
figure('units','normalized','outerposition',[0.5 0 0.6 1])
subplot(5,1,1)
plot(t_total, sqrt(F_cxs.^2 + F_cys.^2)); grid on
xlabel('Time [s]')
ylabel('Resulting force [N]')

posi = get(gca, 'Position');
posi(1) = 0.055;
posi(3) = 0.9;
set(gca, 'Position', posi)

j = 2;
for i = 7:2:size(y_total,2)
  subplot(5, 1, j)
  plot(t_total,y_total(:,i));

  posi = get(gca, 'Position');
  posi(1) = 0.055;
  posi(3) = 0.9;
  set(gca, 'Position', posi)

  j = j + 1;
end


return

% FFT %
% Length of signal (round up to nearest thousand)
L = 1000*(ceil(length(t_total)/1000));

% Equidistant time and data
ts = linspace(tspan(1), tspan(2), L);
timesig = interp1(t_total, pos1(2,:), ts);

% Sampling frequency
Fs = L/(tspan(2)-tspan(1));

% Frequency vector
f = Fs*(0:(L/2))/L;

Y = fft(timesig);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure()
plot(f, P1); grid on
title('Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 50])

