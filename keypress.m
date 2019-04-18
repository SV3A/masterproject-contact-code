function keypress(~, event, h, t_total)
  disp(event.Key);
  event

  if strcmp(event.Key, 'shift')
    if strcmp(event.EventName, 'KeyPress')
      set(h, 'SliderStep', [20/length(t_total) 20/length(t_total)]);

    elseif strcmp(event.EventName, 'KeyRelease')
      set(h, 'SliderStep', [1/length(t_total) 1/length(t_total)]);

    end
  end
end

