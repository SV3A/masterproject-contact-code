% scraps
return

% Simulation time [s]
s = 6;

% Time step and number of integration points
dt = 1e-7;
n = s/dt;

% Initial conditions
gamma_nm1 = 0;
beta_nm1  = 0;
theta_nm1 = 0;

gamma_d = 0;
beta_d = 0;
theta_d = omega;


T_gamma = T_gf(gamma_nm1);
T_beta  = T_bf(beta_nm1);
T_theta = T_tf(theta_nm1);


n_exp = 0;  % Export vector counter

gamma = gamma_nm1 + gamma_d*dt + (1/2)*0*dt^2;
beta  = beta_nm1  + beta_d*dt  + (1/2)*0*dt^2;
theta = theta_nm1 + theta_d*dt + (1/2)*0*dt^2;

  
for i = 1:n
  gamma_dd = -(l_OD * cos(gamma) * m0 * e * theta_d ^ 2 * sin(theta) + l_OM ^ 2 * cos(gamma) * K_my * gamma - I_xx * sin(beta) * beta_d * gamma_d - sin(beta) * gamma_d * I_yy * beta_d + I_zz * sin(beta) * beta_d * gamma_d - l_OG * cos(gamma) * m_tot * g + I_zz * beta_d * theta_d) / I_xx / cos(beta);
  beta_dd = (sin(beta) * sin(theta) * sin(gamma) * e * l_OD * m0 * theta_d ^ 2 + cos(beta) * cos(theta) * e * l_OD * m0 * theta_d ^ 2 + gamma * K_my * sin(beta) * sin(gamma) * l_OM ^ 2 - sin(beta) * gamma_d ^ 2 * I_xx * cos(beta) + I_zz * cos(beta) * sin(beta) * gamma_d ^ 2 - K_mx * cos(beta) * beta * l_OM ^ 2 - l_OG * sin(beta) * sin(gamma) * m_tot * g + I_zz * cos(beta) * gamma_d * theta_d) / I_yy;
  theta_dd = (beta_d * I_xx * cos(beta) * gamma_d - cos(beta) * gamma_d * I_yy * beta_d - I_zz * cos(beta) * beta_d * gamma_d + I_zz * sin(beta) * (l_OD * cos(gamma) * m0 * e * theta_d ^ 2 * sin(theta) + l_OM ^ 2 * cos(gamma) * K_my * gamma - I_xx * sin(beta) * beta_d * gamma_d - sin(beta) * gamma_d * I_yy * beta_d + I_zz * sin(beta) * beta_d * gamma_d - l_OG * cos(gamma) * m_tot * g + I_zz * beta_d * theta_d) / I_xx / cos(beta)) / I_zz;

%   gamma_d_np1 = gamma_d + gamma_dd*dt;
%   gamma_np1 = gamma + gamma_d*dt;
% 
%   beta_d_np1 = beta_d + beta_dd*dt;
%   beta_np1 = beta + beta_d*dt;
% 
%   theta_d_np1 = theta_d + theta_dd*dt;
%   theta_np1 = theta + theta_d*dt;

  gamma_np1 = 2*gamma - gamma_nm1 + gamma_dd*dt^2; 
  beta_np1  = 2*beta  - beta_nm1  + beta_dd*dt^2; 
  theta_np1 = 2*theta - theta_nm1 + theta_dd*dt^2;
  
  gamma_d_np1 = (gamma_np1 - gamma_nm1)/(2*dt);
  beta_d_np1  = (beta_np1 - beta_nm1)/(2*dt);
  theta_d_np1 = (theta_np1 - theta_nm1)/(2*dt);
 
  gamma_nm1 = gamma;
  beta_nm1 = beta;
  theta_nm1 = theta;
  
  gamma = gamma_np1;
  beta  = beta_np1;
  theta = theta_np1;
  
  gamma_d = gamma_d_np1;
  beta_d  = beta_d_np1;
  theta_d = theta_d_np1;
  
  % Update transformation matrices
  T_gamma = T_gf(gamma);
  T_beta  = T_bf(beta);
  T_theta = T_tf(theta);
  
  % Export iteration
  if mod(i,100) == 0
    n_exp = n_exp + 1;
    pos1(:,n_exp) = T_gamma.' * (T_beta.'* [0;0;l_OM] );
    gammaexp(n_exp) = gamma;
    betaexp(n_exp) = beta;
  end
  
end

figure()
plot(linspace(0,s,n_exp), gammaexp);hold on
plot(linspace(0,s,n_exp), betaexp)


figure()
subplot(1,2,1)
plot(linspace(0,s,n_exp),pos1(1,:)); grid on; hold on
plot(linspace(0,s,n_exp),pos1(2,:));

xlabel('Time [s]')
% legend('gamma','beta','theta')

subplot(1,2,2)
plot(pos1(1,:),pos1(2,:))
xlabel('Orbit')

