classdef Rotorsystem < handle
  % 'Rotorsystem' is a class containing all properties relevant to a particular
  % rotor-dynamical system.  The resulting object is pass-by-reference.
  % The state vector substitution is given as:
  %   y = [gamma, gamma_d, beta, beta_d, theta, theta_d,
  %        x_ih, x_ih_d, y_ih, y_ih_d, x_mh, x_mh_d, y_mh, y_mh_d]^T

  properties (Constant)
    r_r = 25e-3/2;      % Rotor radius [m]
    r_s = 29e-3/2;      % Stator radius [m]

    e_x =  0.027;       % Unbalance parameter [m]
    e_y = -0.072;       % Unbalance parameter [m]
    m0  = 1e-3;         % Unbalance mass [kg]

    l_OM = 171.7e-3;    % Position vector, pivot point to the PMB [m]
    l_OG = 195.9e-3;    % Position vector, pivot point to centre of gravity [m]
    l_OD = 259.7e-3;    % Position vector, pivot point to the disc [m]
    l_OC = 428.7e-3;    % Position vector, pivot point to contact point [m]
    l_OE = 72.5e-3;     % Position vector, pivot point to excitation point [m]

    I_xx = 7.085371e-2; % Mass moment of inertia component [kg*m^2]
    I_yy = 7.085371e-2; % Mass moment of inertia component [kg*m^2]
    I_zz = 0.106132e-2; % Mass moment of inertia component [kg*m^2]

    k_xx = 3.09e4;      % Stiffness of the magnetic bearing in x [N/m]
    k_yy = 3.09e4;      % Stiffness of the magnetic bearing in y [N/m]
    k_xy = -2.456e3;    % Cross stiffness term [N/m]
    k_yx = -292;        % Cross stiffness term [N/m]

    d0_xx = 40.9;       % Damping coefficient in x at theta = 0 [N*s/m]
    d_xx  = 8.48;       % Damping coefficient in x [N*s/m]
    d_yy  = 7.7;        % Damping coefficient in y [N*s/m]

    l_OIH = 428.7e-3;   % Position vector, pivot point to inner house [m]
    l_OMH = 428.7e-3;   % Position vector, pivot point to middle [m]

    m_ih = 1.79;        % Inner house mass [kg]
    m_mh = 8.44;        % Middle house mass [kg]

    k_ft1 = 8.42e6;     % Stiffness of the force transducer [N/m]
    k_ft2 = 1.26e7;     % Stiffness of the force transducer [N/m]

    k_vg = 9.709e8;     % Stiffness of vertical beams [N/m]
    k_hg = 1.41e9;      % Stiffness of horizontal beams [N/m]

    d_vg = 2168;        % Damping coefficient of vertical beams [Ns/m]
    d_hg = 12983;       % Damping coefficient of horizontal beams [Ns/m]

  end


  properties
    cl  % Clearance [m]

    % External magnet

    mag_enabled; % Boolean for knowing if the magnet has been enabled
    mag_app_t; % Absolute time, rel. to sim. start, of when to apply the magnet
    mag_app_angle; % Angle to sync the magnet application with
    mag_forcedata; % Force data nx2 vector containing time [s] force [N]
    mag_flag;  % Boolean used to switch on the magnet wrt. angle
  end


  methods
    function obj = Rotorsystem(mag_app_t, mag_app_angle, mag_forcedata)
      % Constructor function.
      % INPUT:

      % Calculate clearance
      obj.cl = obj.r_s - obj.r_r;

      % Handle external excitation
      if nargin > 0

        % Set trigger time and magnet flag
        obj.mag_enabled   = true;
        obj.mag_flag      = false;
        obj.mag_app_t     = mag_app_t;
        obj.mag_app_angle = mag_app_angle;
        obj.mag_forcedata = mag_forcedata;
      else
        obj.mag_enabled   = false;
      end
    end


    function T_gamma = T_gam(~, gamma)
      % Transformation matrix from I to B1

      T_gamma = [1       0          0
                 0  cos(gamma) sin(gamma)
                 0 -sin(gamma) cos(gamma)];
    end


    function T_beta = T_bet(~, beta)
      % Transformation matrix from B1 to B2

      T_beta = [cos(beta) 0 -sin(beta)
                    0     1      0
                sin(beta) 0  cos(beta)];
    end


    function T_theta = T_the(~, theta)
      % Transformation matrix from B2 to B3

      T_theta = [ cos(theta) sin(theta) 0
                 -sin(theta) cos(theta) 0
                       0          0     1];
    end


    function r_OR = rot_centrepos(obj, y)
      % Gets position of the rotor centre in the contact plane in the inertial
      % coordinate system.
      %
      % INPUT:
      %   y: State vector.

      r_OR = obj.T_gam(y(1))' * (obj.T_bet(y(3))' * [0; 0; obj.l_OC]);
    end


    function indent = calc_indent(obj, y)
      % Calculates the size of the indent between the rotor and stator.
      %
      % INPUT:
      %   y: State vector.
      %
      % OUTPUT:
      %   indent: The size of the 'indent', where the indent is represented as
      %   the rotor centre position minus the clearance, thus the output is
      %   negative when the rotor is not in contact with the stator, and
      %   positive otherwise, during which the output signifies the penetration
      %   of the rotor into the stator material.

      % Get rotor centre position
      r_OR = obj.rot_centrepos(y);

      % Calculate the effective clearance (or penetration)
      indent = sqrt( (r_OR(1)-y(7))^2 + (r_OR(2)-y(9))^2 ) - obj.cl;
    end


    function alpha = contact_ang(obj, y)
      % Calculates the contact angle between houses and rotor.
      % INPUT:
      %   y: State vector.

      % Get rotor centre position
      r_OR = obj.rot_centrepos(y);

      alpha = atan2(r_OR(2)-y(9), r_OR(1)-y(7));
    end


    function v_rc = rotor_linvel(obj, y)
      % Calculates the linear velocity of the rotor centre in
      % the inertial coordinate system.
      % INPUT:
      %   y: State vector.

      % Absolute angular velocity of the reference frame B2
      Omega_B2 = obj.T_bet(y(3)) * obj.T_gam(y(1)) * [y(2); 0; 0] + ...
                 obj.T_bet(y(3)) * [0; y(4); 0];

      % Absolute linear velocity of the rotor centre in the contact plane
      v_rc = obj.T_gam(y(1))' * (obj.T_bet(y(3))' * ...
             cross( Omega_B2, [0; 0; obj.l_OC] ));
    end


    function v_rel = abs_rel_velocity(obj, y)
      % Calculates the relative velocity between the rotor and
      % the stator in the inertial coordinate system.
      % INPUT:
      %   y: State vector.

      v_rel = obj.rotor_linvel(y) - [y(7); y(9); 0];
    end


    function vt_rel = tan_rel_velocity(obj, y)
      % Calculates the tangential compoent of the relative
      % velocity between the rotor and stator.
      % INPUT:
      %   y: State vector.

      v_rel = obj.abs_rel_velocity(y);

      % Get contact angle between rotor and inner house
      alpha = obj.contact_ang(y);

      vt_rel = obj.r_r*y(6) - v_rel(1)*sin(alpha) + v_rel(2)*cos(alpha);
    end


    function [delta_d, v_rel_r] = pen_rate(obj, y)
      % Calculates the penetration rate.
      % INPUT:
      %   y: State vector.

      % Get contact angle between rotor and inner house
      alpha = obj.contact_ang(y);

      % Unit vector in the radial direction
      e_r = [cos(alpha); sin(alpha); 0];

      % Relative velocity between the rotor and the inner house (stator)
      v_rel = obj.abs_rel_velocity(y);

      % Projection onto the radial direction
      v_rel_r = dot(v_rel, e_r)*e_r;

      % Determine direction of the identation rate
      if sign(v_rel_r(1)) == sign(e_r(1)) && sign(v_rel_r(2)) == sign(e_r(2))
        dir =  1;
      else
        dir = -1;
      end

      % Penetration rate, i.e. the magnitude of the radial velocity
      delta_d = norm( v_rel_r )*dir;
    end


    function [F_excx, F_excy] = magnetForce(obj, t, y)
      % Delivers the force from the external electro magnet.

      % If the time 't' is between the application time and the current time +
      % the length of the external magnet force data, apply the magnet (if
      % theta = 0) else if rotating, then check the angle.
      % When the angle crosses the threashhold switch the magnet on and
      % interpolate the force from the external data with respect to the 'local'
      % time, i.e. 't - t_magapply'.
      if t >= obj.mag_app_t && t < obj.mag_app_t + obj.mag_forcedata(end, 1)

        if ~obj.mag_flag
          if abs(y(6)) > 0
            % Calculate the instantanious angle
            angle = 360 * (abs(y(5))/(2*pi) - floor(abs(y(5))/(2*pi)));

            % Apply magnet when the angle crosses the specified angle
            if angle > obj.mag_app_angle && angle < obj.mag_app_angle+2 && ...
               ~obj.mag_flag

              obj.mag_flag = true;

              % Print feedback to console
              fprintf('Applied the magnet at %.2f˚\n', angle)
            end
          else
            obj.mag_flag = true;
            fprintf('Applied the magnet at t = %.2f\n', t)
          end
        end

        % When magnet applied interpolate force
        if obj.mag_flag
          F_excx = interp1(obj.mag_forcedata(:, 1), obj.mag_forcedata(:, 2), ...
                           t - obj.mag_app_t);
          F_excy = 0;
          return
        end
      end
      F_excx = 0;
      F_excy = 0;

    end
  end % methods
end % class
