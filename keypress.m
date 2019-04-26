function keypress(~, event, sens_but_handle, sld_handle, t_total, te_total)
% 'keypress' is a function defining shortcut actions

  % Make the 'f' key toggle the sensitivity button
  if strcmp(event.Key, 'f')
    if sens_but_handle.Value == 1
      set(sens_but_handle, 'Value', 0)
    else
      set(sens_but_handle, 'Value', 1)
    end
  elseif strcmp(event.Key, 's')
    sld_handle.Value = 1;
  elseif strcmp(event.Key, 'n')
    stepImpact('Next', sld_handle, t_total, te_total)
  elseif strcmp(event.Key, 'p')
    stepImpact('Previous', sld_handle, t_total, te_total)
  end
end

