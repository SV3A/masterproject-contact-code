function updatefig(h, s, y_total, t_total, te_total, pos1, rotor, stator, ...
                   rotor_cross, stator_cross, rev_mark, force_plt, ...
                   forcein_plt,  fplt, F_r, v_vector,  v_rel_r, ...
                   F_cxs, F_cys, deltas, delta_ds, delta_d_ini, ...
                   s_total, info_field)
% 'updatefig' is used by debug1 to update the debug plot

  % Get plot index from slider
  i = round(get(h, 'Value'));

  % Determine new positions
  rotor_x  = s.r_r*cos(linspace(0,2*pi)) + pos1(1,i);
  rotor_y  = s.r_r*sin(linspace(0,2*pi)) + pos1(2,i);
  stator_x = s.r_s*cos(linspace(0,2*pi)) + y_total(i,7);
  stator_y = s.r_s*sin(linspace(0,2*pi)) + y_total(i,9);

  % Update rotor cross
  rotorcr_x = [pos1(1,i), 7e-4*cos(y_total(i,5)       )+pos1(1,i), ...
               pos1(1,i), 7e-4*cos(y_total(i,5)+0.5*pi)+pos1(1,i), ...
               pos1(1,i), 7e-4*cos(y_total(i,5)+pi    )+pos1(1,i), ...
               pos1(1,i), 7e-4*cos(y_total(i,5)+1.5*pi)+pos1(1,i)];
  rotorcr_y = [pos1(2,i), 7e-4*sin(y_total(i,5)       )+pos1(2,i), ...
               pos1(2,i), 7e-4*sin(y_total(i,5)+0.5*pi)+pos1(2,i), ...
               pos1(2,i), 7e-4*sin(y_total(i,5)+pi    )+pos1(2,i), ...
               pos1(2,i), 7e-4*sin(y_total(i,5)+1.5*pi)+pos1(2,i)];

  % Update revolution mark
  rev_x = [(s.r_r/1.3)*cos(y_total(i,5))+pos1(1,i), ...
           s.r_r*cos(y_total(i,5))+ pos1(1,i)];
  rev_y = [(s.r_r/1.3)*sin(y_total(i,5))+ pos1(2,i), ...
           s.r_r*sin(y_total(i,5))+ pos1(2,i)];

  % Force plot data
  if i > 200 && i < (length(t_total) - 200)
    set(force_plt, 'XData', t_total(i-200:i), 'YData', F_r(i-200:i))
    set(fplt, 'XLim', [t_total(i)-1e-3 t_total(i)+1e-3])
  end

  % Update
  set(rotor, 'XData', rotor_x,  'YData', rotor_y);
  set(stator, 'XData', stator_x, 'YData', stator_y);
  set(rotor_cross, 'XData', rotorcr_x, 'YData', rotorcr_y);
  set(stator_cross, 'XData', y_total(i,7), 'YData', y_total(i,9));
  set(rev_mark, 'XData', rev_x, 'YData', rev_y);

  % Debug values
  if s_total(i) == 0, state = 'False'; else, state = 'True';  end

  info_str = sprintf([' DEBUG VALUES:\n\n time = %.5f s\n F_cx = %7.1f N\n', ...
                      ' F_cy = %7.1f N\n F_r  = %7.1f N\n\n', ...
                      ' delta   = %10.7f\n delta_d = %10.7f\n', ...
                      ' deltain = %10.7f\n\n Contact: %s'], ...
                      t_total(i), F_cxs(i), F_cys(i), F_r(i), deltas(i), ...
                      delta_ds(i), delta_d_ini(i), state);

  set(info_field, 'string', info_str);

  % Velocity vector
  % (XSTART, YSTART) (XDIR, YDIR)
  vector = [ pos1(1,i)  pos1(2,i)  v_rel_r(1,i) v_rel_r(2,i) ] ;

  set(v_vector,'XData',vector( :, 1 ), 'YData', vector( :, 2), ...
               'UData',vector( :, 3 ), 'VData', vector( :, 4))

  % Force-indentation plots
  % Find how many impact events that have occured so far
  im_idx = length(find(te_total <= t_total(i)));

  % If the number above is uneven we're currently in contact
  if mod(im_idx, 2) == 1

    % Start and end indices for the contact state change
    start_idx = find(t_total == te_total(im_idx));

    % Check weather simulation ended during impact
    if length(te_total) > im_idx
      end_idx = find(t_total == te_total(im_idx+1));
    else
      end_idx = length(t_total);
    end

    % Update plot
    forcein_x = deltas(start_idx:end_idx);
    forcein_y = F_r(start_idx:end_idx);
    set(forcein_plt, 'XData', forcein_x, 'YData', forcein_y)
  end

  drawnow
end

