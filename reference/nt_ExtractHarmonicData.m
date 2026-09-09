
%% =========================================================================
%  Acoustic Emissions Analysis — FUS BBB Opening
%  Single-pass burst extraction: 2nd harmonic + wideband
% =========================================================================


%% 1. Load data
load('Mouse_Cntr_03_Target2.mat')


%% 2. Aliases
% raw  = data.posx.data.rawdata;
ns    = data.niscope;
trk   = data.niscope.Tracker;
awg1  = data.AWG.AWG1;
awg2  = data.AWG.AWG2;
raw   = data.niscope;

%% 3. Key parameters
fs        = ns.SampleRate;
f_axis    = ns.freqaxis;
t_axis    = ns.timeaxis;
df        = f_axis(2) - f_axis(1);
N_bursts  = length(raw.ff);
PRF       = data.posx.curPRF;
t_burst   = (0:N_bursts-1) / PRF;
HarmGoal  = data.posx.HarmGoal;
f_carrier = awg1.Frequency;
f_2nd     = 2 * f_carrier;
f_wb      = 1700000;

N_baseline = 5;       % bursts before bubble injection — adjust as needed
f2_half_win   = 20;       % ±20 bins = 41 points per window
WB_half_win   = 60;       % ±60 bins = 121 points per window

fprintf('========================================\n')
fprintf('  Sample rate:       %.2f MHz\n',   fs/1e6)
fprintf('  Freq resolution:   %.4f Hz\n',    df)
fprintf('  FFT length:        %d pts\n',     length(f_axis))
fprintf('  Nyquist:           %.2f MHz\n',   f_axis(end)/1e6)
fprintf('  Bursts:            %d\n',         N_bursts)
fprintf('  PRF:               %.1f Hz\n',    PRF)
fprintf('  Duration:          %.0f s\n',     t_burst(end))
fprintf('  AWG1:              %.4f MHz, %.3f V\n', awg1.Frequency/1e6, awg1.Voltage)
fprintf('  AWG2:              %.4f MHz, %.3f V\n', awg2.Frequency/1e6, awg2.Voltage)
fprintf('  Carrier:           %.4f MHz\n',   f_carrier/1e6)
fprintf('  2nd harmonic:      %.4f MHz\n',   f_2nd/1e6)
fprintf('  Wideband monitor:  %.4f MHz\n',   f_wb/1e6)
fprintf('  2f Window:         ±%d bins (%.1f Hz)\n', f2_half_win, f2_half_win*df)
fprintf('  WB Window:         ±%d bins (%.1f Hz)\n', WB_half_win, WB_half_win*df)
fprintf('  Baseline bursts:   1 to %d\n',   N_baseline)
fprintf('========================================\n')

%% 4. Define frequency windows (index-based)
[~, idx_2nd] = min(abs(f_axis - f_2nd));
[~, idx_wb]  = min(abs(f_axis - f_wb));

win_2nd = (idx_2nd - f2_half_win) : (idx_2nd + f2_half_win);
win_wb  = (idx_wb  - WB_half_win) : (idx_wb  + WB_half_win);

fprintf('2nd harmonic window: %.4f – %.4f MHz (%d pts)\n', ...
    f_axis(win_2nd(1))/1e6, f_axis(win_2nd(end))/1e6, length(win_2nd))
fprintf('Wideband window:     %.4f – %.4f MHz (%d pts)\n', ...
    f_axis(win_wb(1))/1e6,  f_axis(win_wb(end))/1e6,  length(win_wb))

%% 5. Single loop over bursts
peak_2nd = zeros(1, N_bursts);
area_2nd = zeros(1, N_bursts);
peak_wb  = zeros(1, N_bursts);
area_wb  = zeros(1, N_bursts);

for ii = 1:N_bursts
    ff_matrix(:,ii)=raw.ff(ii).data0;
    spectrum = abs(raw.ff(ii).data0);

    % 2nd harmonic
    seg_2nd      = spectrum(win_2nd);
    peak_2nd(ii) = max(seg_2nd);
    area_2nd(ii) = sum(seg_2nd) * df;

    % Wideband
    seg_wb      = spectrum(win_wb);
    peak_wb(ii) = max(seg_wb);
    area_wb(ii) = sum(seg_wb) * df;

end
ff_matrix = single(ff_matrix);

size(data.posx.data,2);
V_log   = data.posx.data(size(data.posx.data,2)).V;

V1_burst  = V_log(:, 1)';   
V2_burst  = V_log(:, 2)';   

%% 6. Baseline normalisation
bl = 1:N_baseline;
peak_2nd_norm = peak_2nd / mean(peak_2nd(bl));
area_2nd_norm = area_2nd / mean(area_2nd(bl));
peak_wb_norm  = peak_wb  / mean(peak_wb(bl));
area_wb_norm  = area_wb  / mean(area_wb(bl));

%% 7. Plots
figure('Name', 'Cavitation Metrics', 'Position', [100 100 900 700])
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
plot(t_burst, 20*log10(peak_2nd + eps), 'b', 'LineWidth', 1.5)
yline(0+20*log10(mean(peak_2nd(bl))+eps), 'k--', 'baseline'); grid on
xlabel('Time (s)'); ylabel('dB re baseline')
title('2nd harmonic — peak')
xlim([0 t_burst(end)])

nexttile
plot(t_burst, 20*log10(area_2nd + eps), 'b', 'LineWidth', 1.5)
yline(0+20*log10(mean(area_2nd(bl))+eps), 'k--', 'baseline'); grid on
xlabel('Time (s)'); ylabel('dB re baseline')
title('2nd harmonic — spectral area')
xlim([0 t_burst(end)])

nexttile
plot(t_burst, 20*log10(peak_wb + eps), 'r', 'LineWidth', 1.5)
yline(0+20*log10(mean(peak_wb(bl))+eps), 'k--', 'baseline'); grid on
xlabel('Time (s)'); ylabel('dB re baseline')
title('Wideband — peak')
xlim([0 t_burst(end)])

nexttile
plot(t_burst, 20*log10(area_wb + eps), 'r', 'LineWidth', 1.5)
yline(0+20*log10(mean(area_wb(bl))+eps), 'k--', 'baseline'); grid on
xlabel('Time (s)'); ylabel('dB re baseline')
title('Wideband — spectral area')
xlim([0 t_burst(end)])

sgtitle(sprintf('FUS BBB — Cavitation Metrics  |  %.4f MHz carrier  |  PRF %.0f Hz  |  %d bursts', ...
    f_carrier/1e6, PRF, N_bursts))

%% 8. Cumulative metrics
% Cumulative sum of normalised area over bursts
cum_area_2nd = cumsum(area_2nd_norm);
cum_area_wb  = cumsum(area_wb_norm);

% Cumulative sum of normalised peak
cum_peak_2nd = cumsum(peak_2nd_norm);
cum_peak_wb  = cumsum(peak_wb_norm);

% Cumulative sum of voltages
cum_V1 = cumsum(V1_burst);
cum_V2 = cumsum(V2_burst);

% Print summary scalars
fprintf('\n========================================\n')
fprintf('  Cumulative metrics (normalised)\n')
fprintf('  2nd harmonic — cumulative area: %.2f\n',  cum_area_2nd(end))
fprintf('  2nd harmonic — cumulative peak: %.2f\n',  cum_peak_2nd(end))
fprintf('  Wideband     — cumulative area: %.2f\n',  cum_area_wb(end))
fprintf('  Wideband     — cumulative peak: %.2f\n',  cum_peak_wb(end))
fprintf('  Harmonic/WB area ratio:         %.2f\n',  cum_area_2nd(end) / cum_area_wb(end))
fprintf('AWG1 voltage range: %.3f to %.3f V\n', min(V1_burst), max(V1_burst))
fprintf('========================================\n')

% Plot
figure('Name', 'Cumulative Cavitation Metrics', 'Position', [100 100 900 400])
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
plot(t_burst, cum_area_2nd, 'b', 'LineWidth', 1.5)
hold on
plot(t_burst, cum_area_wb,  'r', 'LineWidth', 1.5)
xlabel('Time (s)'); ylabel('Cumulative sum (normalised)')
title('Cumulative spectral area')
legend('2nd harmonic (stable)', 'Wideband (inertial)')
grid on; xlim([0 t_burst(end)])

nexttile
plot(t_burst, cum_peak_2nd, 'b', 'LineWidth', 1.5)
hold on
plot(t_burst, cum_peak_wb,  'r', 'LineWidth', 1.5)
xlabel('Time (s)'); ylabel('Cumulative sum (normalised)')
title('Cumulative peak')
legend('2nd harmonic (stable)', 'Wideband (inertial)')
grid on; xlim([0 t_burst(end)])

sgtitle(sprintf('Cumulative Cavitation Dose  |  %.4f MHz carrier  |  %d bursts', ...
    f_carrier/1e6, N_bursts))

%% 9. Save reduced data file

% ── Original filename (strip .mat if present, append _reduced) ────────────
% original_file   = FILES(ppp).name;
% [fpath, fname, ~] = fileparts(original_file);
% Freq_file    = fullfile(fpath, ['FreqOnly_' fname '.mat']);
% reduced_file    = fullfile(fpath, ['Reduced_' fname '.mat']);
% 
% % ── Save -----------------------------------------------------------------
% % save(Freq_file, 'ff_matrix')
% clear data raw ns trk awg1 awg2 data f_axis ff_matrix ns raw spectrum t_axis trk ff_matrix
% save(reduced_file)
% 
% clear area* cum* bl f2_* f_* fs HarmGoal peak_* seg_* V1_* V2_* V_* WB_* win_* idx* N_* posx*

