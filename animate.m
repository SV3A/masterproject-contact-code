function [] = animate(t_total,y_total,s_total,te_total,s,F_cxs,F_cys,deltas)
% 'animate' animates the numerical solution.

  % Simulation start and end timestamps
  t_start = t_total(1);
  t_end   = t_total(end);

  % Interpolate the data for efficient animation
  intr_points = 30000;
  %intr_points = length(t_total);
  t_plot = linspace(t_start, t_end, intr_points);
  y_plot = ones(length(t_plot), 14);

  for i = 1:14
    y_plot(:,i) = interp1(t_total, y_total(:,i), t_plot);
  end

  F_cx = interp1(t_total, F_cxs, t_plot);
  F_cy = interp1(t_total, F_cys, t_plot);
  F_r = sqrt(F_cx.^2 + F_cy.^2);


  % Build orbit in I
  pos1 = zeros(3, length(t_plot));

  for i = 1:length(t_plot)
    T_gamma = s.T_gam(y_plot(i,1));
    T_beta  = s.T_bet(y_plot(i,3));
    pos1(:,i) = T_gamma' * (T_beta'*[0; 0; s.l_OC]);
  end

  figure('units','normalized','outerposition',[0 0 1 1])
  subplot(2,2,[1,3]); hold on
  % Trajectory trail
  tjec = animatedline('MaximumNumPoints',4*intr_points*(11/1000),...
                      'LineStyle',':','LineWidth',1.75);

  % Rotor and stator circles
  rotor  = plot(0,0,'r','LineWidth',1.5); grid on
  stator = plot(0,0,'k','LineWidth',1.2);

  % Center cross markers
  %rotor_cross = plot(pos1(1,1), pos1(2,1),'r+',...
                      %'LineWidth',1.5,'MarkerSize',20);
  rotor_cross = plot(0,0,'r','LineWidth',1.5);

  stator_cross = plot(y_plot(1,7), y_plot(1,9),'k+',...
                      'LineWidth',1.5,'MarkerSize',20);

  % Revolution mark
  rev_mark = plot(0,0,'r','LineWidth',2);

  % Define axes
  axis(1*[-0.016 0.016 -0.016 0.016]);
  axis equal
  xticklabels([]); yticklabels([]);

  % Time stamp
  time_stamp = text(-14e-3, 20e-3, '', 'FontSize', 15);
  force_stamp = text(-14e-3, 15e-3, '', 'FontSize', 16);

  %dir_stamp = text(-14e-3, 15e-3, '', 'FontSize', 16);

  % Force plot
  fplt = subplot(2,2,2);
  force_plt = animatedline('MaximumNumPoints',300,'LineWidth',1.5, ...
                           'Color', 'blue'); grid on
  xlabel('Time [s]')
  ylabel('Force [N]')
  fplt.YLim = [0 max(F_r)];

  % Force-indentation plot
  subplot(2,2,4);
  forcein_plt = plot(0,0,'k','LineWidth',1.5); grid on
  xlabel('Indentation [mm]')
  ylabel('Normal force [N]')

  % Animation loop
  im_idx = 1; % Impact index

  for i = 1:length(t_plot)
    % Update trajectory trail
    addpoints(tjec, pos1(1,i), pos1(2,i));
    addpoints(force_plt, t_plot(i), F_r(i));

    % Determine new positions
    rotor_x  = s.r_r*cos(linspace(0,2*pi)) + pos1(1,i);
    rotor_y  = s.r_r*sin(linspace(0,2*pi)) + pos1(2,i);
    stator_x = s.r_s*cos(linspace(0,2*pi)) + y_plot(i,7);
    stator_y = s.r_s*sin(linspace(0,2*pi)) + y_plot(i,9);
    % Update rotor cross
    rotorcr_x = [pos1(1,i), 7e-4*cos(y_plot(i,5)       )+pos1(1,i), ...
                 pos1(1,i), 7e-4*cos(y_plot(i,5)+0.5*pi)+pos1(1,i), ...
                 pos1(1,i), 7e-4*cos(y_plot(i,5)+pi    )+pos1(1,i), ...
                 pos1(1,i), 7e-4*cos(y_plot(i,5)+1.5*pi)+pos1(1,i)];
    rotorcr_y = [pos1(2,i), 7e-4*sin(y_plot(i,5)       )+pos1(2,i), ...
                 pos1(2,i), 7e-4*sin(y_plot(i,5)+0.5*pi)+pos1(2,i), ...
                 pos1(2,i), 7e-4*sin(y_plot(i,5)+pi    )+pos1(2,i), ...
                 pos1(2,i), 7e-4*sin(y_plot(i,5)+1.5*pi)+pos1(2,i)];
    % Update revolution mark
    rev_x = [(s.r_r/1.3)*cos(y_plot(i,5))+pos1(1,i), ...
             s.r_r*cos(y_plot(i,5))+ pos1(1,i)];
    rev_y = [(s.r_r/1.3)*sin(y_plot(i,5))+ pos1(2,i), ...
             s.r_r*sin(y_plot(i,5))+ pos1(2,i)];

    % Update
    set(rotor, 'XData', rotor_x,  'YData', rotor_y);
    set(stator, 'XData', stator_x, 'YData', stator_y);
    set(rotor_cross, 'XData', rotorcr_x, 'YData', rotorcr_y);
    set(stator_cross, 'XData', y_plot(i,7), 'YData', y_plot(i,9));
    set(rev_mark, 'XData', rev_x, 'YData', rev_y);
    set(time_stamp,'String', ['t = ', num2str(t_plot(i)), ' s'])
    %set(dir_stamp,'String', ['dir = ', num2str(dirs(i))])


    % Force-indentation plots
    if te_total(im_idx+1) <= t_plot(i)
      start_idx = find(t_total == te_total(im_idx));
      end_idx = find(t_total == te_total(im_idx+1));

      forcein_x = deltas(start_idx:end_idx);
      forcein_y = sqrt(F_cxs(start_idx:end_idx).^2+F_cys(start_idx:end_idx).^2);

      set(forcein_plt, 'XData', forcein_x, 'YData', forcein_y)

      im_idx = im_idx + 2;
    end

    drawnow limitrate

    %pause(1)
  end

end

