classdef Nikravesh < Contactmodels
% 'Nikravesh' is a subclass of 'Contactmodels' and defines the contact model set
% forth in the article:
%   "A Contact Force Model With Hysteresis Damping for Impact Analysis of
%   Multibody Systems"
%   In: Journal of Mechanical Design, Year 1990, Volume 112, Issue 3, pp. 369.
%   By: Lankarani, H. M. and Nikravesh, P. E.

  properties (Constant)
    E_r  = 70e9;  % Young's modulus sphere (rotor material part) [Pa]
    E_s  = 100e9; % Young's modulus sphere (stator material part) [Pa]
    nu_r = 0.35;  % Poisson's ratio for the rotor material [-]
    nu_s = 0.35;  % Poisson's ratio for the stator material [-]
  end


  properties
    n             % Penetration exponent
    ce            % Coefficient of restitution
    K             % Generalized material parameter
    delta_d_init  % Relative velocity at impact instant
  end


  methods

    function obj = Nikravesh(friction_model, r_s, r_r)
      % Constructor function.
      % INPUT:
      %   r_s: stator inner radius
      %   r_r: rotor radius

      % Send friction model to super class constructor
      obj = obj@Contactmodels(friction_model);

      obj.name = "Nikravesh";
      obj.print_name;
      obj.n = 3/2;
      obj.ce = 0.5;
      %obj.K = (4/( 3*(((1-obj.nu_r^2)/obj.E_r) + ((1-obj.nu_s^2)/obj.E_s))))*...
      %(r_r*r_s / (r_s - r_r))^(0.5);

      % Material parameters
      h_r = (1 - obj.nu_r^2)/(pi*obj.E_r);
      h_s = (1 - obj.nu_s^2)/(pi*obj.E_s);

      obj.K = ( 4/(3*pi*(h_r + h_s)) ) * ((r_r*r_s)/(r_r+r_s))^(0.5);
      obj.delta_d_init = 0;
    end


    function Fn = calc_fn(obj, d, d_dot)
      % Calculates the magnitude of the radial normal force
      % INPUT:
      %   d    : penetration
      %   d_dot: relative velocity between rotor and stator

      Fn = obj.K * d^obj.n * ( 1 + (3*(1-obj.ce^2)/4)*d_dot/obj.delta_d_init );

    end

  end
end % class
