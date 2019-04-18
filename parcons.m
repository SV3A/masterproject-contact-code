function parcons
% 'parcons' sets all physical parameters and constants.

  assignin('base', 'g', 0);        % Acceleration of gravity [m/s^2]

  % Shaft ----------------------------------------------------------------------
  % Position vectors [m]
  assignin('base', 'l_OM', 171.7e-3); % From pivot point to the PMB
  assignin('base', 'l_OG', 195.9e-3); % From pivot point to the centre of
                                      % gravity of the shaft
  assignin('base', 'l_OD', 259.7e-3); % From pivot point to the disc

  assignin('base', 'm_tot', 1.486);   % Shaft mass (including discs) [kg]

  % Mass moment of inertia components [kg*m^2]
  assignin('base', 'I_xx', 73399360 * 1e-9);
  assignin('base', 'I_yy', 73399360 * 1e-9);
  assignin('base', 'I_zz',   981365 * 1e-9);

  % Stiffness of the magnetic bearing [N/m]
  assignin('base', 'K_mx', 3.4675*10^4);
  assignin('base', 'K_my', 3.4675*10^4);

  % Damping ratio [-]
  assignin('base', 'xi', 0.000);

  % Houses ---------------------------------------------------------------------
  % Position vectors [m]
  assignin('base', 'l_OIH', 0e-3); % From pivot point to inner house
  assignin('base', 'l_OMH', 0e-3); % From pivot point to middle

  % House masses [kg]
  assignin('base', 'm_ih', 1.7);
  assignin('base', 'm_mh', 7.17);

  % Stiffness of the force transducer [N/m]
  assignin('base', 'K_ft', 8.88e7);

  % Stiffness of the beams [N/m]
  assignin('base', 'K_vg', 1.23e7);
  assignin('base', 'K_hg', 2.40e7);

  % Damping coefficient of the beams [Ns/m]
  assignin('base', 'D_vg', 238);
  assignin('base', 'D_hg', 1576);

end
