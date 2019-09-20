function y_dot = dydt(~, y, s, cmod, state)
% 'dydt' serves the equation of motion as a system of six first order equations.
%
% INPUT:
%   y     : State vector
%   s     : Rotor system object
%   cmod  : The contact model object
%   state : Contact indicator (0|1)
% OUTPUT:
%   y_dot: Accelerations and velocities
%

  l_OD  = s.l_OD;
  l_OM  = s.l_OM;
  l_OC  = s.l_OC;
  I_xx  = s.I_xx;
  I_yy  = s.I_yy;
  I_zz  = s.I_zz;
  K_mx  = s.K_mx;
  K_my  = s.K_my;
  K_ft  = s.K_ft;
  K_vg  = s.K_vg;
  K_hg  = s.K_hg;
  D_x   = s.D_x;
  D_y   = s.D_y;
  D_vg  = s.D_vg;
  D_hg  = s.D_hg;
  m0    = s.m0;
  e     = s.e;
  m_ih  = s.m_ih;
  m_mh  = s.m_mh;


  % Calculate contact forces
  [F_cx, F_cy] = contactForce(y, s, cmod, state);


% Don't edit anything after this as it may be overwritten by Maple!
% * * * * *
  y_dot(1)  = y(2);
  y_dot(2)  = -(cos(y(5)) * e * l_OD * m0 * y(6) ^ 2 + sin(y(5)) * e * l_OD * m0 * y(6) ^ 2 + l_OM ^ 2 * cos(y(1)) * D_y * y(2) + l_OM ^ 2 * cos(y(1)) * K_my * y(1) - I_xx * sin(y(3)) * y(4) * y(2) - sin(y(3)) * y(2) * I_yy * y(4) + I_zz * sin(y(3)) * y(4) * y(2) - l_OC * cos(y(1)) * F_cy + I_zz * y(4) * y(6)) / I_xx / cos(y(3));
  y_dot(3)  = y(4);
  y_dot(4)  = -(-D_y * sin(y(1)) * sin(y(3)) * l_OM ^ 2 * y(2) - y(1) * K_my * sin(y(1)) * sin(y(3)) * l_OM ^ 2 - cos(y(5)) * e * l_OD * m0 * y(6) ^ 2 + sin(y(5)) * e * l_OD * m0 * y(6) ^ 2 + D_x * cos(y(3)) * y(4) * l_OM ^ 2 + sin(y(3)) * y(2) ^ 2 * I_xx * cos(y(3)) - I_zz * cos(y(3)) * sin(y(3)) * y(2) ^ 2 + K_mx * cos(y(3)) * y(3) * l_OM ^ 2 + F_cy * sin(y(1)) * sin(y(3)) * l_OC - I_zz * cos(y(3)) * y(2) * y(6) + F_cx * cos(y(3)) * l_OC) / I_yy;
  y_dot(5)  = y(6);
  y_dot(6)  = (I_zz * cos(y(5)) * sin(y(3)) * e * l_OD * m0 * y(6) ^ 2 + I_zz * sin(y(5)) * sin(y(3)) * e * l_OD * m0 * y(6) ^ 2 + D_y * I_zz * cos(y(1)) * sin(y(3)) * l_OM ^ 2 * y(2) + y(1) * I_zz * K_my * cos(y(1)) * sin(y(3)) * l_OM ^ 2 + cos(y(3)) ^ 2 * y(2) * y(4) * I_xx ^ 2 - cos(y(3)) ^ 2 * y(2) * y(4) * I_xx * I_yy - cos(y(3)) ^ 2 * y(2) * y(4) * I_xx * I_zz - I_xx * I_zz * sin(y(3)) ^ 2 * y(4) * y(2) - I_yy * I_zz * sin(y(3)) ^ 2 * y(4) * y(2) + I_zz ^ 2 * sin(y(3)) ^ 2 * y(4) * y(2) - F_cy * I_zz * cos(y(1)) * sin(y(3)) * l_OC + I_zz ^ 2 * sin(y(3)) * y(4) * y(6)) / I_xx / cos(y(3)) / I_zz;
  y_dot(7)  = y(8);
  y_dot(8)  = -0.1e1 / m_ih * (D_vg * y(8) - D_vg * y(12) + K_vg * y(7) - K_vg * y(11) - F_cx);
  y_dot(9)  = y(10);
  y_dot(10) = (-2 * K_ft * y(9) + 2 * K_ft * y(13) + F_cy) / m_ih;
  y_dot(11) = y(12);
  y_dot(12) = (D_vg * y(8) - D_vg * y(12) - 2 * K_ft * y(11) + K_vg * y(7) - K_vg * y(11)) / m_mh;
  y_dot(13) = y(14);
  y_dot(14) = -1 / m_mh * (D_hg * y(14) - 2 * K_ft * y(9) + 2 * K_ft * y(13) + K_hg * y(13));

  y_dot = y_dot';

end
% Maple part last updated: 26-08-2019 17:38
