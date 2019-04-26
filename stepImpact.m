function stepImpact(dir, sld_handle, t_total, te_total)
% 'stepImpact' is a callback function to jump to the next/previous impacts in a
% debug gui.

  % Get plot index from slider
  i = round(get(sld_handle, 'Value'));

  % Find number of impact events that have occured so far
  im_idx = length(find(te_total <= t_total(i)));

  % Determine steps, in terms of number of events, to increase/decrease
  if mod(im_idx, 2) == 0 || im_idx == 0
    step_idx = 1;
  else
    step_idx = 2;
  end

  % Jump forward or backward
  if strcmp(dir, 'Next')
    % If not at the end
    if length(te_total) > im_idx + step_idx
      new_value = find( t_total == te_total(im_idx + step_idx) );
    else
      uicontrol(sld_handle)
      return
    end
  elseif strcmp(dir, 'Previous')
    % If not at the beginning
    if im_idx >= 2
      new_value = find( t_total == te_total(im_idx - step_idx) );
    else
      uicontrol(sld_handle)
      return
    end
  end

  % Update slider
  sld_handle.Value = new_value;

  % Focus slider
  uicontrol(sld_handle)
end

