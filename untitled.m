`matlab
%% EEE 4186 - DSP Lab 2: Sampling, Quantization and Aliasing
% Requires: Signal Processing Toolbox (buttord, butter, freqs), MATLAB R2018b+ (xline/yline)
clear; clc; close all;

%% 1. Parameters
f0       = 900;          % tone frequency (Hz)
fs_safe  = 5000;         % safe sampling rate (Hz)
fs_alias = 1000;         % under-sampling rate (Hz)
fs_sim   = 50000;        % "continuous-time" simulation rate (Hz)
T        = 0.02;         % observation window (s)

t  = (0:round(T*fs_sim)-1)/fs_sim;     x  = sin(2*pi*f0*t);
n1 = (0:round(T*fs_safe)-1)/fs_safe;   x1 = sin(2*pi*f0*n1);
n2 = (0:round(T*fs_alias)-1)/fs_alias; x2 = sin(2*pi*f0*n2);

%% 2. Pre-lab calculations
fm = 3.5e3;  fprintf('Q1 Nyquist rate      = %g Hz\n', 2*fm);
fa = abs(900 - 1000*round(900/1000));  fprintf('Q2 Alias frequency   = %g Hz\n', fa);
fprintf('Q3 Nyquist frequency = %g Hz\n', 8000/2);

%% 3. Figure 1 - safe sampling (5 kHz)
figure(1);
plot(t*1e3, x, 'k', 'LineWidth', 1.2); hold on;
stem(n1*1e3, x1, 'b', 'MarkerFaceColor', 'b');
grid on; xlabel('Time (ms)'); ylabel('Amplitude');
title('Safe sampling: fs = 5 kHz');
legend('900 Hz tone', 'Samples at 5 kHz', 'Location', 'northeast');

%% 4. Figure 2 - under-sampling (1 kHz) with the inverted 100 Hz alias
figure(2);
plot(t*1e3, x, 'Color', [0.3 0.3 0.3]); hold on;
plot(t*1e3, -sin(2*pi*100*t), 'g--', 'LineWidth', 1.4);
stem(n2*1e3, x2, 'r', 'MarkerFaceColor', 'r');
grid on; xlabel('Time (ms)'); ylabel('Amplitude');
title('Under-sampling: fs = 1 kHz');
legend('900 Hz original', '100 Hz alias (inverted)', 'Samples at 1 kHz');
fprintf('Max |x2 - inverted 100 Hz sine| = %g\n', max(abs(x2 + sin(2*pi*100*n2))));

%% 5. Figure 3 - FFT of both sampled sequences
Nfft = 8192;
[fA, XA] = mag_spectrum(x1, fs_safe,  Nfft);
[fB, XB] = mag_spectrum(x2, fs_alias, Nfft);
[~, iA] = max(XA);  [~, iB] = max(XB);
fprintf('FFT peak (5 kHz) = %.1f Hz,  FFT peak (1 kHz) = %.1f Hz\n', fA(iA), fB(iB));

figure(3);
subplot(1,2,1); plot(fA, XA, 'b'); grid on; xlim([0 fs_safe/2]);
xlabel('Frequency (Hz)'); ylabel('Normalised |X(f)|');
title(sprintf('fs = 5 kHz (peak %.0f Hz)', fA(iA)));
subplot(1,2,2); plot(fB, XB, 'r'); grid on; xlim([0 fs_alias/2]);
xlabel('Frequency (Hz)'); ylabel('Normalised |X(f)|');
title(sprintf('fs = 1 kHz (peak %.0f Hz)', fB(iB)));

%% 6. Figure 4 - frequency-folding map for fs = 1 kHz
f   = 0:3000;
fal = abs(f - fs_alias*round(f/fs_alias));
figure(4);
plot(f, fal, 'r', 'LineWidth', 1.2); hold on;
plot(900, 100, 'ko', 'MarkerFaceColor', 'k');
yline(fs_alias/2, '--');
grid on; xlabel('True input frequency f_0 (Hz)'); ylabel('Apparent frequency f_a (Hz)');
title('Frequency folding for fs = 1 kHz');

%% 7. Quantization (mid-rise, full scale +/-1 V)
bits = [2 4 6 8 10 12];
rms_meas = zeros(size(bits));  rms_th = rms_meas;  sqnr = rms_meas;  step = rms_meas;
for k = 1:numel(bits)
    [q, d]      = quantize_midrise(x, bits(k), 1);
    e           = q - x;
    step(k)     = d;
    rms_meas(k) = sqrt(mean(e.^2));
    rms_th(k)   = d/sqrt(12);
    sqnr(k)     = 10*log10(mean(x.^2)/mean(e.^2));
end
T_q = table(bits(:), 2.^bits(:), step(:), rms_meas(:), rms_th(:), sqnr(:), ...
    'VariableNames', {'Bits','Levels','Step','RMS_measured','RMS_theory','SQNR_dB'});
disp(T_q);

% Figure 5 - 4-bit and 8-bit quantisation of the 5 kHz samples
figure(5);
bl = [4 8];
for p = 1:2
    b = bl(p);
    [q, d] = quantize_midrise(x1, b, 1);
    subplot(2,1,p);
    plot(t*1e3, x, 'Color', [0.6 0.6 0.6]); hold on;
    plot(n1*1e3, q, 's', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    plot([n1; n1]*1e3, [x1; q], 'r');       % red error lines
    xlim([0 4]); ylim([-1.2 1.2]); grid on; ylabel('Amplitude');
    title(sprintf('%d-bit quantiser, step = %.4f V', b, d));
end
xlabel('Time (ms), first 4 ms shown');

% Figure 6 - error waveform and RMS error versus bit depth
figure(6);
subplot(1,2,1); hold on;
for b = [4 8]
    q = quantize_midrise(x, b, 1);
    plot(t*1e3, q - x);
end
grid on; xlabel('Time (ms)'); ylabel('Error'); legend('4-bit', '8-bit');
title('Quantisation error');
subplot(1,2,2);
semilogy(bits, rms_meas, 'o-'); hold on;
semilogy(bits, rms_th, 'k--');
grid on; xlabel('Bits'); ylabel('RMS error'); legend('Measured', '\Delta/\surd12');
title('RMS error vs bit depth');

%% 8. Figure 7 - anti-aliasing filter (passband 300 Hz/1 dB, stopband 700 Hz/40 dB)
fp = 300; fstop = 700; Ap = 1; As = 40;
[N, Wn]  = buttord(2*pi*fp, 2*pi*fstop, Ap, As, 's');
[z,p,kk] = butter(N, Wn, 's');
[bA, aA] = zp2tf(z, p, kk);
ff = linspace(1, 1500, 3000);
H  = 20*log10(abs(freqs(bA, aA, 2*pi*ff)));
fprintf('Filter order N = %d, -3 dB freq = %.0f Hz, attenuation at 900 Hz = %.1f dB\n', ...
    N, Wn/(2*pi), -interp1(ff, H, 900));

figure(7);
plot(ff, H, 'b', 'LineWidth', 1.2); hold on; grid on;
xline(fp, ':g'); xline(fstop, ':r'); xline(fs_alias/2, 'k:');
ylim([-80 5]); xlabel('Frequency (Hz)'); ylabel('Gain (dB)');
title(sprintf('Butterworth anti-aliasing filter, order %d', N));

% Demonstration: filter the 900 Hz tone, then sample at 1 kHz
t_long = (0:round(0.2*fs_sim)-1)/fs_sim;
[zd, pd, kd] = butter(N, (Wn/(2*pi))/(fs_sim/2));
sos = zp2sos(zd, pd, kd);
y = sin(2*pi*f0*t_long);
for s = 1:size(sos,1)                         % cascade of biquads
    y = filter(sos(s,1:3), sos(s,4:6), y);
end
y_s = y(1:fs_sim/fs_alias:end);
fprintf('Peak after filter + 1 kHz sampling = %.5f (%.0f dB)\n', ...
    max(abs(y_s(end-9:end))), 20*log10(max(abs(y_s(end-9:end)))));

%% ---------------- local functions ----------------
function [fr, X] = mag_spectrum(seq, fs, Nfft)
    L = numel(seq);
    w = 0.5 - 0.5*cos(2*pi*(0:L-1)/(L-1));      % Hann window
    S = abs(fft(seq .* w, Nfft));
    X = S(1:Nfft/2+1) / max(S(1:Nfft/2+1));
    fr = (0:Nfft/2) * fs / Nfft;
end

function [q, d] = quantize_midrise(sig, b, vmax)
    L   = 2^b;
    d   = 2*vmax/L;
    idx = min(max(floor(sig/d), -L/2), L/2 - 1);
    q   = (idx + 0.5)*d;
end
```

This code box has no hidden characters, so it won't cause the error from before. If MATLAB shows a new error, send me a screenshot.