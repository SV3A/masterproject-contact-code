function debug1( t_total, y_total, te_total, s_total, s, cmod )
% 'animate' animates the numerical solution.

  % Simulation start and end timestamps
  %t_start = t_total(1);
  %t_end   = t_total(end);

  % Build orbit in I, get contact forces, retrieve the contact angle and get
  % relative radial velocity

  pos1     = zeros( 3, length(t_total) );
  alphas   = zeros( length(t_total), 1 );
  v_rel_r  = zeros( 3, length(t_total) );
  F_cxs    = zeros( length(t_total), 1 );
  F_cys    = zeros( length(t_total), 1 );
  deltas   = zeros( length(t_total), 1 );
  delta_ds = zeros( length(t_total), 1 );
  delta_d_ini = zeros( length(t_total), 1 );

  for i = 1:length(t_total)
    y_i = y_total(i,:)';

    T_gamma = s.T_gam(y_total(i,1));
    T_beta  = s.T_bet(y_total(i,3));
    pos1(:,i) = T_gamma' * (T_beta'*[0; 0; s.l_OC]);

    if s_total(i) == 0, state = 0; else, state = 1;  end

    [F_cxs(i), F_cys(i), ~, delta_d_ini(i)] = ...
      contact_force(y_i, s, cmod, state);

    alphas(i) = s.contact_ang(y_i);
    deltas(i) = s.calc_gap(y_i);
    [delta_ds(i), v_rel_r(:,i)] = s.pen_rate(y_i);
  end

  % Scale down velocity vector
  v_rel_r = v_rel_r*1/10;

  % Calculate resultant radial contact force
  F_r = sqrt(F_cxs.^2 + F_cys.^2);


  mainfig = figure('units','normalized','outerposition',[0 0 1 1]);
  % Rotor and stator circles
  subplot(2,2,[1,3]); hold on

  rotor  = plot(0,0,'r','LineWidth',1.5); grid on
  stator = plot(0,0,'k','LineWidth',1.2);

  % Center cross markers
  rotor_cross = plot(0,0,'r','LineWidth',1.5);

  stator_cross = plot(y_total(1,7), y_total(1,9),'k+',...
                      'LineWidth',1.5,'MarkerSize',20);

  % Revolution mark
  rev_mark = plot(0,0,'r','LineWidth',2);

  % Quiver plot
  v_vector = quiver( 0, 0, 0, 0, ...
                     'LineWidth', 1.2, 'AutoScale', 'on', ...
                     'MaxHeadSize', 0.9 );

  % Define axes
  axis(1*[-0.016 0.016 -0.016 0.016]);
  axis equal
  xticklabels([]); yticklabels([]);

  % Time stamp

  % Force plot
  fplt = subplot(2,2,2);
  force_plt = plot(0,0,'b','LineWidth',1.5); grid on

  xlabel('Time [s]')
  ylabel('Force [N]')
  fplt.YLim = [0 max(F_r)];


  % Force-indentation plot
  subplot(2,2,4);
  forcein_plt = plot(0,0,'k','LineWidth',1.5); grid on
  xlabel('Indentation [mm]')
  ylabel('Normal force [N]')

  % Add info text field
  info_field = uicontrol('style', 'text', ...
                         'backgroundcolor', [1 1 1], ...
                         'fontname', 'Monospaced', ...
                         'fontsize', 14, ...
                         'horizontalalignment', 'left', ...
                         'units', 'normalized', ...
                         'position', [0.01, 0.58, 0.12, 0.36]);

  % Add slider control
  uicontrol('style','text', 'string', 'Time step', 'fontsize', 14, ...
            'fontname', 'Monospaced', ...
            'horizontalalignment', 'left', ...
            'units','normalized', 'position',[0.015 0.04 0.1 0.05])

  h = uicontrol('style','slider',...
                'backgroundcolor', [1 1 1], ...
                'units','normalized',...
                'position',[0.01 0.04 0.4 0.01],...
                'min',1,'max',length(t_total),'val',1,...
                'SliderStep', [1/length(t_total) 20/length(t_total)]);

  % Create eventlistener for a change of slider value
  addlistener(h,'Value','PostSet', @(~,~) updatefig(h, s, y_total, t_total, ...
              te_total, pos1, rotor, stator, rotor_cross, stator_cross, ...
              rev_mark, force_plt, forcein_plt, fplt, F_r, v_vector, ...
              v_rel_r, F_cxs, F_cys, deltas, delta_ds, delta_d_ini, ...
              s_total, info_field));

  %set(mainfig,'KeyPressFcn', @(src, event) keypress(src, event, h, t_total),...
              %'KeyReleaseFcn', @(src, event) keypress(src, event, h, t_total));

end

