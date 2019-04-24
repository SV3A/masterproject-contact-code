function calcFft(t_total, tspan, pos1)
% FFT

  % Length of signal (round up to nearest thousand)
  L = 1000*(ceil(length(t_total)/1000));

  % Equidistant time and data
  ts = linspace(tspan(1), tspan(2), L);
  timesig = interp1(t_total, pos1(2,:), ts);

  % Sampling frequency
  Fs = L/(tspan(2)-tspan(1));

  % Frequency vector
  f = Fs*(0:(L/2))/L;

  Y = fft(timesig);
  P2 = abs(Y/L);
  P1 = P2(1:L/2+1);
  P1(2:end-1) = 2*P1(2:end-1);

  figure()
  plot(f, P1); grid on
  title('Single-Sided Amplitude Spectrum')
  xlabel('f (Hz)')
  ylabel('|P1(f)|')
  xlim([0 50])
end

