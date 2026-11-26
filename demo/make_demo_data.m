function truth = make_demo_data(outDir)
%MAKE_DEMO_DATA  Synthetic transient-absorption data for TRSpecAnalyzer.
%
%   make_demo_data             writes the files into the current folder
%   make_demo_data(outDir)     writes them into outDir
%   truth = make_demo_data(..) also returns every "true" quantity
%
%   Files written
%     demo_TA.csv                    the data set of the walkthrough
%     scans/demo_scan_01..05.csv     five repeated scans (Load & Average)
%     truth/demo_truth_spectra.csv   true species spectra, oscillation
%                                    profiles, t0(lambda), background
%     truth/demo_truth_kinetics.csv  true populations (no chirp)
%     demo_target.csv                second data set, for Target Analysis
%     truth/demo_target_truth_spectra.csv, demo_target_truth_kinetics.csv
%     demo_equilibrium.csv           third data set, for Target Analysis
%                                    of a reversible step
%     truth/demo_equilibrium_truth_spectra.csv,
%     truth/demo_equilibrium_truth_kinetics.csv
%     demo_solution.csv, demo_solvent.csv
%                                    a solute in cyclohexane and
%                                    cyclohexane alone, for the solvent
%                                    coherence and for Subtract data
%     truth/demo_solvent_truth_spectra.csv, demo_solvent_truth_modes.csv
%     demo_solution_raw.csv,         the same two measurements before the
%     demo_solvent_raw.csv           corrections: with chirp and a
%                                    background
%     truth/demo_solvent_raw_truth.csv
%     demo_aniso_par.csv, demo_aniso_perp.csv
%                                    one sample measured with parallel
%                                    and perpendicular polarization
%     truth/demo_aniso_truth_spectra.csv, demo_aniso_truth_kinetics.csv
%     demo_aniso_raw_par.csv,        the same sample as measured: with
%     demo_aniso_raw_perp.csv        background, pump scatter and chirp
%     demo_aniso_raw_perp_late.csv   the perpendicular measurement once
%                                    more, its time zero 0.05 ps late
%     truth/demo_aniso_raw_truth.csv
%     demo_solvation.csv             time-resolved emission of a dye in
%                                    a polar solvent: the band shifts to
%                                    the red while the solvent relaxes
%     truth/demo_solvation_truth_spectra.csv,
%     truth/demo_solvation_truth_kinetics.csv
%
%   File layout (the format TRSpecAnalyzer reads): first row = delays,
%   first column = wavelengths, the rest = signal (wavelength x delay).
%   Units: wavelength in nm, delay in ps, signal in mOD.
%
%   What is in the data
%     kinetics     sequential  A -> B -> C -> ground state
%                  tau = 0.60 ps, 15 ps, 700 ps
%     IRF          Gaussian, FWHM 0.060 ps
%     chirp        time zero depends on the wavelength: group delay of
%                  8 mm fused silica, t0 = 0 at 500 nm (blue is earlier)
%     oscillations 150 cm-1 (damping 1.0 ps), node at 580 nm
%                  250 cm-1 (damping 0.8 ps), node at 655 nm
%     background   constant offset and pump scatter at 500 nm, present
%                  at all delays (also before time zero)
%     noise        white, larger at both ends of the spectrum and at
%                  the pump wavelength
%   No coherent artifact is included.
%
%   What is in demo_target.csv (already corrected: no chirp, no background)
%     kinetics     A -> B (3 ps),  A -> C (6 ps),  B -> GS (40 ps),
%                  C -> GS (800 ps);  A decays in 2 ps, 67 % to B, 33 % to C
%     IRF          Gaussian, FWHM 0.10 ps, time zero 0
%     spectra      from 450 to 475 nm all three states show the same
%                  ground-state bleach and nothing else
%
%   What is in demo_equilibrium.csv (already corrected as well)
%     kinetics     A <-> B (5 ps forward, 15 ps back),  B -> C (100 ps),
%                  C -> GS (2000 ps);  equilibrium constant kf/kb = 3.
%                  The data show 3.71 ps, 135 ps and 2000 ps
%     IRF          Gaussian, FWHM 0.10 ps, time zero 0
%     spectra      B has no signal from 700 to 750 nm
%
%   What is in demo_solution.csv and demo_solvent.csv (no chirp, no
%   background; 5 fs steps up to 5 ps)
%     solute       A -> B (0.40 ps) -> ground state (80 ps); wavepacket
%                  modes at 215 and 590 cm-1; no signal above 725 nm
%     solvent      cyclohexane: a coherent artifact at time zero and its
%                  Raman modes at 384.1, 426.3, 801.3, 1028.3, 1157.6,
%                  1266.4 and 1444.4 cm-1 (the one at 801.3 cm-1 by far
%                  the strongest), node at 565 nm
%     IRF          Gaussian, FWHM 0.020 ps
%     solution     solute + 0.70 x the solvent signal, 0.003 ps later
%
%   What is in demo_solution_raw.csv and demo_solvent_raw.csv
%     the two measurements above as they would come from the instrument
%     (own noise), on the same delay axis:
%     chirp        the same for both files: 1.5 mm of fused silica, time
%                  zero 0 at 565 nm (-0.091 ps at 450 nm, +0.062 ps at
%                  750 nm)
%     background   constant in time, and not the same in the two files:
%                  an offset (0.25 and 0.10 mOD) and scattered pump light
%                  around 520 nm (-1.5 and -0.6 mOD at its maximum)
%
%   What is in demo_aniso_par.csv and demo_aniso_perp.csv (corrected data)
%     kinetics     A -> B (3 ps) -> ground state (800 ps)
%     anisotropy   0.38 for bleach and stimulated emission, -0.12 for the
%                  excited-state absorption; rotation time 50 ps
%     G            the perpendicular file is too small by the factor 1.15
%     IRF          Gaussian, FWHM 0.10 ps, time zero 0
%
%   The noise comes from a small generator written out below, not from
%   randn, so MATLAB and Octave produce the same files.
%
%   Runs in MATLAB (R2020b or later) and in GNU Octave; no toolboxes.

if nargin < 1 || isempty(outDir), outDir = pwd; end

% ------------------------------------------------------------------
% True parameters
% ------------------------------------------------------------------
P.tau      = [0.60 15 700];   % ps, A -> B -> C -> ground state
P.fwhm     = 0.060;           % ps, Gaussian IRF (FWHM)
P.chirpL   = 8;               % mm of fused silica
P.wlRef    = 500;             % nm, wavelength at which t0 = 0
P.oscFreq  = [150 250];       % cm-1
P.oscDamp  = [1.0 0.8];       % ps
P.oscAmp   = [0.50 0.60];     % mOD, before the IRF is applied
P.oscNode  = [580 655];       % nm, node of the amplitude profile
P.oscWidth = [22 30];         % nm
P.pumpWL   = 500;             % nm, pump scatter
P.bgOffset = 0.12;            % mOD, constant offset
P.bgScat   = -4.0;            % mOD, pump scatter at its maximum
P.noise    = 0.03;            % mOD, demo_TA.csv at 600 nm
P.noiseScan = 0.07;           % mOD, single scans at 600 nm

% ------------------------------------------------------------------
% Axes
% ------------------------------------------------------------------
wl = (450:2:750).';                               % nm
tPre  = (-20:-11) * 0.1;                          % -2.0 ... -1.1 ps
tFine = (-50:200) * 0.02;                         % -1.0 ... 4.0 ps
tLog  = logspace(log10(4), log10(3000), 71);      % 4 ps ... 3 ns
tLog  = tLog(2:end);
dec   = 10 .^ (floor(log10(tLog)) - 2);           % three significant digits
tLog  = round(tLog ./ dec) .* dec;
t = [tPre, tFine, tLog];
t = round(t * 1e6) / 1e6;
if any(diff(t) <= 0), error('Delay axis is not strictly increasing.'); end

% ------------------------------------------------------------------
% Species spectra (mOD)
% ------------------------------------------------------------------
g   = @(c, w) exp(-0.5 * ((wl - c) / w) .^ 2);
gsb = -(10.0 * g(520, 16) + 5.5 * g(487, 14));    % ground-state bleach
SA  = gsb - 6.0 * g(570, 18) + 8.0 * g(680, 45);  % A: hot excited state
SB  = gsb - 4.5 * g(585, 20) + 9.0 * g(650, 32);  % B: relaxed excited state
SC  = gsb + 6.5 * g(610, 24) + 2.0 * g(465, 15);  % C: long-lived state
S   = [SA, SB, SC];

% ------------------------------------------------------------------
% Time zero of every wavelength (chirp)
% ------------------------------------------------------------------
t0 = chirpCurve(wl, P.chirpL, P.wlRef);           % ps

% ------------------------------------------------------------------
% Oscillations: profile(lambda) * damped cosine, IRF applied numerically
% ------------------------------------------------------------------
nOsc = numel(P.oscFreq);
oscProf = zeros(numel(wl), nOsc);
for m = 1:nOsc
    x = (wl - P.oscNode(m)) / P.oscWidth(m);
    oscProf(:, m) = P.oscAmp(m) * x .* exp(0.5 - 0.5 * x .^ 2);   % +-1 at x = +-1
end
sig   = P.fwhm / (2 * sqrt(2 * log(2)));
dtf   = 0.001;                                    % 1 fs grid
tf    = (-1000:14000) * dtf;                      % -1 ... 14 ps
ker   = exp(-0.5 * (((-300:300) * dtf) / sig) .^ 2);
ker   = ker / sum(ker);
cLight = 2.99792458e-2;                           % cm / ps
oscT  = zeros(numel(tf), nOsc);
for m = 1:nOsc
    y = cos(2 * pi * cLight * P.oscFreq(m) * tf) .* exp(-tf / P.oscDamp(m));
    y(tf < 0) = 0;
    y(tf == 0) = 0.5;
    oscT(:, m) = conv(y, ker, 'same').';
end

% ------------------------------------------------------------------
% Background (independent of the delay)
% ------------------------------------------------------------------
bg = P.bgOffset + P.bgScat * g(P.pumpWL, 2.5);

% Noise level of every wavelength (relative to the level at 600 nm)
noiseRel = (1 + 2 * ((wl - 600) / 150) .^ 2) .* (1 + 4 * g(P.pumpWL, 3));

% ------------------------------------------------------------------
% Build the data sets
% ------------------------------------------------------------------
clean = @(shift, scale) buildClean(wl, t, S, P, t0 + shift, sig, ...
    oscProf, oscT, tf) * scale + bg;

D = clean(0, 1) + P.noise * noiseRel .* demoNoise(numel(wl), numel(t), 1001);

if ~exist(outDir, 'dir'), mkdir(outDir); end
scanDir  = fullfile(outDir, 'scans');
truthDir = fullfile(outDir, 'truth');
if ~exist(scanDir, 'dir'),  mkdir(scanDir);  end
if ~exist(truthDir, 'dir'), mkdir(truthDir); end

writeData(fullfile(outDir, 'demo_TA.csv'), wl, t, D);

% Five repeated scans. Scan 3: time zero 0.06 ps late. Scan 5: signal
% reduced to 65 % (a degraded sample), to be left out of the average.
scanShift = [0 0 0.06 0 0];
scanScale = [1 1 1 1 0.65];
for k = 1:5
    Dk = clean(scanShift(k), scanScale(k)) + ...
        P.noiseScan * noiseRel .* demoNoise(numel(wl), numel(t), 2000 + k);
    writeData(fullfile(scanDir, sprintf('demo_scan_%02d.csv', k)), wl, t, Dk);
end

% ------------------------------------------------------------------
% Truth files
% ------------------------------------------------------------------
fid = fopen(fullfile(truthDir, 'demo_truth_spectra.csv'), 'w');
fprintf(fid, ['wavelength_nm,S_A_mOD,S_B_mOD,S_C_mOD,osc_150cm-1_mOD,' ...
    'osc_250cm-1_mOD,t0_ps,background_mOD,noise_demo_TA_mOD\n']);
fprintf(fid, '%g,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', ...
    [wl, S, oscProf, t0, bg, P.noise * noiseRel].');
fclose(fid);

Ctrue = seqConc(t(:), P.tau, 0, sig);
fid = fopen(fullfile(truthDir, 'demo_truth_kinetics.csv'), 'w');
fprintf(fid, 'delay_ps,c_A,c_B,c_C\n');
fprintf(fid, '%.6g,%.6f,%.6f,%.6f\n', [t(:), Ctrue].');
fclose(fid);

writeTargetDemo(outDir, truthDir);
writeEquilibriumDemo(outDir, truthDir);
writeSolventDemo(outDir, truthDir);
writeAnisotropyDemo(outDir, truthDir);
writeAnisotropyRawDemo(outDir, truthDir);
writeSolvationDemo(outDir, truthDir);

if nargout > 0
    truth = P;
    truth.wl = wl;  truth.t = t;  truth.S = S;  truth.t0 = t0;
    truth.oscProf = oscProf;  truth.bg = bg;  truth.noiseRel = noiseRel;
    truth.scanShift = scanShift;  truth.scanScale = scanScale;
end
fprintf('Demo data written to %s\n', outDir);
end

% ======================================================================
function writeTargetDemo(outDir, truthDir)
% Second data set: a branching scheme, for Target Analysis.
%   A -> B : k1 (3 ps)      A -> C : k2 (6 ps)
%   B -> GS : k3 (40 ps)    C -> GS : k4 (800 ps)
tau  = [3 6 40 800];              % time constant of every step, ps
fwhm = 0.10;                      % ps
sig  = fwhm / (2 * sqrt(2 * log(2)));
wl = (450:2:750).';
tLin = (-20:40) * 0.05;                           % -1 ... 2 ps
tLog = logspace(log10(2), log10(5000), 61);
tLog = tLog(2:end);
dec  = 10 .^ (floor(log10(tLog)) - 2);
tLog = round(tLog ./ dec) .* dec;
t = round([tLin, tLog] * 1e6) / 1e6;
if any(diff(t) <= 0), error('Delay axis is not strictly increasing.'); end

g  = @(c, w) exp(-0.5 * ((wl - c) / w) .^ 2);
gsb = -(9.0 * g(470, 14) + 4.0 * g(445, 11));     % same bleach for all states
SA  = gsb - 5.0 * g(545, 18) + 7.0 * g(690, 35);
SB  = gsb + 8.0 * g(600, 26);
SC  = gsb + 6.0 * g(730, 30) + 2.5 * g(550, 18);
S   = [SA, SB, SC];

k  = 1 ./ tau;
kA = k(1) + k(2);
EA = expIRF(t(:), kA, sig);
C  = [EA, ...
      k(1) / (k(3) - kA) * (EA - expIRF(t(:), k(3), sig)), ...
      k(2) / (k(4) - kA) * (EA - expIRF(t(:), k(4), sig))];

noiseRel = 1 + 2 * ((wl - 600) / 150) .^ 2;
D = S * C.' + 0.03 * noiseRel .* demoNoise(numel(wl), numel(t), 3001);
writeData(fullfile(outDir, 'demo_target.csv'), wl, t, D);

fid = fopen(fullfile(truthDir, 'demo_target_truth_spectra.csv'), 'w');
fprintf(fid, 'wavelength_nm,S_A_mOD,S_B_mOD,S_C_mOD\n');
fprintf(fid, '%g,%.4f,%.4f,%.4f\n', [wl, S].');
fclose(fid);
fid = fopen(fullfile(truthDir, 'demo_target_truth_kinetics.csv'), 'w');
fprintf(fid, 'delay_ps,c_A,c_B,c_C\n');
fprintf(fid, '%.6g,%.6f,%.6f,%.6f\n', [t(:), C].');
fclose(fid);
end

% ======================================================================
function writeEquilibriumDemo(outDir, truthDir)
% Third data set: a reversible step, for Target Analysis.
%   A <-> B : kf (5 ps), kb (15 ps)
%   B -> C  : k2 (100 ps)
%   C -> GS : k3 (2000 ps)
tau  = [5 15 100 2000];           % time constant of kf, kb, k2, k3, ps
fwhm = 0.10;                      % ps
sig  = fwhm / (2 * sqrt(2 * log(2)));
wl = (450:2:750).';
tLin = (-20:40) * 0.05;                           % -1 ... 2 ps
tLog = logspace(log10(2), log10(5000), 61);
tLog = tLog(2:end);
dec  = 10 .^ (floor(log10(tLog)) - 2);
tLog = round(tLog ./ dec) .* dec;
t = round([tLin, tLog] * 1e6) / 1e6;
if any(diff(t) <= 0), error('Delay axis is not strictly increasing.'); end

g  = @(c, w) exp(-0.5 * ((wl - c) / w) .^ 2);
gsb = -(9.0 * g(495, 15) + 4.5 * g(463, 12));     % ground-state bleach
SA  = gsb - 4.0 * g(555, 16) + 8.0 * g(710, 30);  % A: band in the red
SB  = gsb + 9.0 * g(605, 22);                     % B: nothing above 700 nm
SC  = gsb + 6.0 * g(535, 18) + 2.5 * g(650, 45);  % C: product
S   = [SA, SB, SC];

% Populations, unit start in A. A and B exchange and B feeds C, so A
% and B decay with the two rates l1 > l2 of their 2 x 2 rate matrix,
% not with kf, kb or k2:
%   l1 + l2 = kf + kb + k2,   l1 * l2 = kf * k2
k  = 1 ./ tau;
kf = k(1);  kb = k(2);  k2 = k(3);  k3 = k(4);
sm = kf + kb + k2;
rt = sqrt(sm ^ 2 - 4 * kf * k2);
l1 = (sm + rt) / 2;
l2 = (sm - rt) / 2;
E1 = expIRF(t(:), l1, sig);
E2 = expIRF(t(:), l2, sig);
E3 = expIRF(t(:), k3, sig);
C  = [((l1 - kb - k2) * E1 + (kb + k2 - l2) * E2) / (l1 - l2), ...
      kf * (E2 - E1) / (l1 - l2), ...
      kf * k2 / (l1 - l2) * ((E2 - E3) / (k3 - l2) - (E1 - E3) / (k3 - l1))];

noiseRel = 1 + 2 * ((wl - 600) / 150) .^ 2;
D = S * C.' + 0.03 * noiseRel .* demoNoise(numel(wl), numel(t), 4001);
writeData(fullfile(outDir, 'demo_equilibrium.csv'), wl, t, D);

fid = fopen(fullfile(truthDir, 'demo_equilibrium_truth_spectra.csv'), 'w');
fprintf(fid, 'wavelength_nm,S_A_mOD,S_B_mOD,S_C_mOD\n');
fprintf(fid, '%g,%.4f,%.4f,%.4f\n', [wl, S].');
fclose(fid);
fid = fopen(fullfile(truthDir, 'demo_equilibrium_truth_kinetics.csv'), 'w');
fprintf(fid, 'delay_ps,c_A,c_B,c_C\n');
fprintf(fid, '%.6g,%.6f,%.6f,%.6f\n', [t(:), C].');
fclose(fid);
end

% ======================================================================
function writeSolventDemo(outDir, truthDir)
% Fourth data set: a solute in cyclohexane (demo_solution.csv) and
% cyclohexane alone (demo_solvent.csv), measured with short pulses.
% For the subtraction of the solvent coherence in the Coherence window
% and for "Subtract data...". No chirp, no background.
tau   = [0.40 80];                % ps, A -> B -> ground state
fwhm  = 0.020;                    % ps
sig   = fwhm / (2 * sqrt(2 * log(2)));
scale = 0.70;                     % part of the solvent signal in the solution
shift = 0.003;                    % ps, the solvent signal comes later there
% Raman modes of cyclohexane: shift (cm-1) and relative intensity
% (R. L. McCreery, standard spectra for Raman shift calibration,
% https://www.chem.ualberta.ca/~mccreery/ramanmaterials.html)
solvFreq = [384.1 426.3 801.3 1028.3 1157.6 1266.4 1444.4];
solvRel  = [2 3 95 15 6 14 12] / 95;
solvDamp = [1.0 1.0 2.0 1.2 1.0 1.0 0.9];          % ps
solvAmp  = 1.2;                   % mOD, 801 cm-1 mode before the IRF is applied
% wavepacket modes of the solute
oscFreq  = [215 590];             % cm-1
oscDamp  = [1.2 0.9];             % ps
oscAmp   = [0.50 0.45];           % mOD, before the IRF is applied
oscNode  = [560 650];             % nm
oscWidth = [20 22];               % nm
noise    = 0.02;                  % mOD at 600 nm

wl = (450:4:750).';                               % nm
tFine = (-60:1000) * 0.005;                       % -0.3 ... 5.0 ps, 5 fs steps
tLog  = logspace(log10(5), log10(300), 26);
tLog  = tLog(2:end);
dec   = 10 .^ (floor(log10(tLog)) - 2);
tLog  = round(tLog ./ dec) .* dec;
t = round([tFine, tLog] * 1e6) / 1e6;
if any(diff(t) <= 0), error('Delay axis is not strictly increasing.'); end

g   = @(c, w) exp(-0.5 * ((wl - c) / w) .^ 2);
gsb = -(8.0 * g(505, 16) + 3.5 * g(475, 13));     % ground-state bleach
SA  = gsb - 4.0 * g(560, 18) + 6.0 * g(650, 28);  % A: excited state
SB  = gsb + 5.0 * g(620, 26);                     % B: relaxed state
S   = [SA, SB];                                   % (no solute signal above 725 nm)
k   = 1 ./ tau;
E1  = expIRF(t(:), k(1), sig);
E2  = expIRF(t(:), k(2), sig);
C   = [E1, k(1) / (k(2) - k(1)) * (E1 - E2)];

% oscillations on a 0.5 fs grid, IRF applied numerically
dtf = 0.0005;
tf  = (-800:24000) * dtf;                         % -0.4 ... 12 ps
ker = exp(-0.5 * (((-120:120) * dtf) / sig) .^ 2);
ker = ker / sum(ker);
cLight = 2.99792458e-2;                           % cm / ps
nS = numel(solvFreq);
solvT = zeros(numel(tf), nS);
for m = 1:nS                                      % impulsive Raman: sine
    y = sin(2 * pi * cLight * solvFreq(m) * tf) .* exp(-tf / solvDamp(m));
    y(tf < 0) = 0;
    solvT(:, m) = conv(y, ker, 'same').';
end
nO = numel(oscFreq);
oscT = zeros(numel(tf), nO);
oscProf = zeros(numel(wl), nO);
for m = 1:nO                                      % excited-state wavepacket: cosine
    y = cos(2 * pi * cLight * oscFreq(m) * tf) .* exp(-tf / oscDamp(m));
    y(tf < 0) = 0;
    y(tf == 0) = 0.5;
    oscT(:, m) = conv(y, ker, 'same').';
    x = (wl - oscNode(m)) / oscWidth(m);
    oscProf(:, m) = oscAmp(m) * x .* exp(0.5 - 0.5 * x .^ 2);
end
% the solvent modes shift the probe spectrum: amplitude with a node at
% the middle of the probe spectrum (565 nm), the same for all modes
x = (wl - 565) / 75;
solvProf = solvAmp * x .* exp(0.5 - 0.5 * x .^ 2);         % +-solvAmp at x = +-1
% coherent artifact of the solvent around time zero (mOD)
sa   = 1.2 * sig;
art0 = 1.25 * (1 + 0.4 * sin((wl - 450) / 50));
art1 = -2.0 * (wl - 565) / 150;

% solvent signal at the delays tq (row vector)
solvent = @(tq) (solvProf * solvRel) * interp1(tf(:), solvT, tq(:), 'linear', 0).' ...
    + art0 * exp(-0.5 * (tq / sa) .^ 2) ...
    + art1 * ((tq / sa) .* exp(-0.5 * (tq / sa) .^ 2));

noiseRel = 1 + 2 * ((wl - 600) / 150) .^ 2;
Dsolv = solvent(t) + noise * noiseRel .* demoNoise(numel(wl), numel(t), 5002);
Dsol  = S * C.' + oscProf * interp1(tf(:), oscT, t(:), 'linear', 0).' ...
      + scale * solvent(t - shift) ...
      + noise * noiseRel .* demoNoise(numel(wl), numel(t), 5001);
writeData(fullfile(outDir, 'demo_solution.csv'), wl, t, Dsol);
writeData(fullfile(outDir, 'demo_solvent.csv'),  wl, t, Dsolv);

fid = fopen(fullfile(truthDir, 'demo_solvent_truth_spectra.csv'), 'w');
fprintf(fid, ['wavelength_nm,S_A_mOD,S_B_mOD,osc_215cm-1_mOD,osc_590cm-1_mOD,' ...
    'solvent_801cm-1_mOD,artifact_even_mOD,artifact_odd_mOD\n']);
fprintf(fid, '%g,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', ...
    [wl, S, oscProf, solvProf, art0, art1].');
fclose(fid);
fid = fopen(fullfile(truthDir, 'demo_solvent_truth_modes.csv'), 'w');
fprintf(fid, 'what,frequency_cm-1,amplitude_relative_to_801,damping_ps\n');
for m = 1:nS
    fprintf(fid, 'cyclohexane,%.1f,%.4f,%.2f\n', solvFreq(m), solvRel(m), solvDamp(m));
end
for m = 1:nO
    fprintf(fid, 'solute,%.1f,,%.2f\n', oscFreq(m), oscDamp(m));
end
fclose(fid);

% ---- The same two measurements before the corrections ----
% (demo_solution_raw.csv, demo_solvent_raw.csv): time zero depends on
% the wavelength (chirp of the probe), and each file has its own
% background. The chirp is that of the probe, so it is the same in both
% files; the solvent signal of the solution is still 0.003 ps late.
t0raw  = chirpCurve(wl, 1.5, 565);                % ps
bgSol  = 0.25 - 1.5 * g(520, 25);                 % mOD
bgSolv = 0.10 - 0.6 * g(520, 25);
DsolR  = zeros(numel(wl), numel(t));
DsolvR = zeros(numel(wl), numel(t));
for i = 1:numel(wl)
    ti = t - t0raw(i);                            % delays as this wavelength sees them
    e1 = expIRF(ti(:), k(1), sig);
    e2 = expIRF(ti(:), k(2), sig);
    ci = [e1, k(1) / (k(2) - k(1)) * (e1 - e2)];
    sv = solvent(ti);
    ss = solvent(ti - shift);
    DsolvR(i, :) = sv(i, :);
    DsolR(i, :)  = (ci * S(i, :).').' ...
        + oscProf(i, :) * interp1(tf(:), oscT, ti(:), 'linear', 0).' ...
        + scale * ss(i, :);
end
DsolR  = DsolR  + bgSol  + noise * noiseRel .* demoNoise(numel(wl), numel(t), 5003);
DsolvR = DsolvR + bgSolv + noise * noiseRel .* demoNoise(numel(wl), numel(t), 5004);
writeData(fullfile(outDir, 'demo_solution_raw.csv'), wl, t, DsolR);
writeData(fullfile(outDir, 'demo_solvent_raw.csv'),  wl, t, DsolvR);

fid = fopen(fullfile(truthDir, 'demo_solvent_raw_truth.csv'), 'w');
fprintf(fid, 'wavelength_nm,t0_ps,background_solution_mOD,background_solvent_mOD\n');
fprintf(fid, '%g,%.5f,%.4f,%.4f\n', [wl, t0raw, bgSol, bgSolv].');
fclose(fid);
end

% ======================================================================
function writeAnisotropyDemo(outDir, truthDir)
% Fifth data set: one sample measured with the probe polarized parallel
% (demo_aniso_par.csv) and perpendicular (demo_aniso_perp.csv) to the
% pump. Already corrected: no background, no chirp, time zero at 0.
tau    = [3 800];                 % ps, A -> B -> ground state
tauRot = 50;                      % ps, rotation of the molecule
r0     = [0.38 -0.12];            % anisotropy at time zero: bleach and
                                  % stimulated emission; excited-state absorption
G      = 1.15;                    % the perpendicular file is too small by this factor
fwhm   = 0.10;                    % ps
sig    = fwhm / (2 * sqrt(2 * log(2)));
wl = (450:2:750).';
tLin = (-20:40) * 0.05;                           % -1 ... 2 ps
tLog = logspace(log10(2), log10(5000), 61);
tLog = tLog(2:end);
dec  = 10 .^ (floor(log10(tLog)) - 2);
tLog = round(tLog ./ dec) .* dec;
t = round([tLin, tLog] * 1e6) / 1e6;
if any(diff(t) <= 0), error('Delay axis is not strictly increasing.'); end

g   = @(c, w) exp(-0.5 * ((wl - c) / w) .^ 2);
gsb = -(8.0 * g(505, 16) + 3.5 * g(475, 13));
% every state: [bleach + stimulated emission, excited-state absorption]
SA  = [gsb - 4.0 * g(565, 18), 7.0 * g(690, 38)];
SB  = [gsb - 3.0 * g(580, 20), 6.0 * g(650, 32)];

% populations, and populations x exp(-t / tauRot): the rate of every
% exponential grows by 1 / tauRot
k   = 1 ./ tau;
kr  = 1 / tauRot;
pop = @(dk) [expIRF(t(:), k(1) + dk, sig), ...
             k(1) / (k(2) - k(1)) * (expIRF(t(:), k(1) + dk, sig) - expIRF(t(:), k(2) + dk, sig))];
C   = pop(0);
CR  = pop(kr);
M   = [sum(SA, 2), sum(SB, 2)] * C.';             % magic angle
R   = [SA * r0(:), SB * r0(:)] * CR.';            % sum of (signal x anisotropy)
par  = M + 2 * R;
perp = (M - R) / G;

noiseRel = 1 + 2 * ((wl - 600) / 150) .^ 2;
par  = par  + 0.03 * noiseRel .* demoNoise(numel(wl), numel(t), 6001);
perp = perp + 0.03 * noiseRel .* demoNoise(numel(wl), numel(t), 6002);
writeData(fullfile(outDir, 'demo_aniso_par.csv'),  wl, t, par);
writeData(fullfile(outDir, 'demo_aniso_perp.csv'), wl, t, perp);

fid = fopen(fullfile(truthDir, 'demo_aniso_truth_spectra.csv'), 'w');
fprintf(fid, ['wavelength_nm,S_A_r0.38_mOD,S_A_r-0.12_mOD,' ...
    'S_B_r0.38_mOD,S_B_r-0.12_mOD\n']);
fprintf(fid, '%g,%.4f,%.4f,%.4f,%.4f\n', [wl, SA, SB].');
fclose(fid);
fid = fopen(fullfile(truthDir, 'demo_aniso_truth_kinetics.csv'), 'w');
fprintf(fid, 'delay_ps,c_A,c_B,c_A_x_rotation,c_B_x_rotation\n');
fprintf(fid, '%.6g,%.6f,%.6f,%.6f,%.6f\n', [t(:), C, CR].');
fclose(fid);
end

% ======================================================================
function writeAnisotropyRawDemo(outDir, truthDir)
% The sample of writeAnisotropyDemo as it would be measured: two files
% (demo_aniso_raw_par.csv, demo_aniso_raw_perp.csv) with a background,
% pump scatter and chirp, on the delay axis of demo_TA.csv. Spectra, time
% constants, anisotropies and G are those of writeAnisotropyDemo.
% A third file (demo_aniso_raw_perp_late.csv) is the perpendicular
% measurement made once more with its time zero 0.05 ps late: every
% feature sits 0.05 ps later than in the parallel file (own noise).
%   chirp        the same for both files: 8 mm of fused silica, time zero
%                0 at 500 nm
%   background   constant in time, and not the same in the two files:
%                an offset (0.12 and 0.20 mOD) and pump scatter at 480 nm
%                (-4.0 and -1.5 mOD at its maximum)
tau    = [3 800];
tauRot = 50;
r0     = [0.38 -0.12];
G      = 1.15;
fwhm   = 0.10;
sig    = fwhm / (2 * sqrt(2 * log(2)));
pumpWL = 480;
tLate  = 0.05;                    % ps, time zero of the third file
wl = (450:2:750).';
tPre  = (-20:-11) * 0.1;                          % -2.0 ... -1.1 ps
tFine = (-50:200) * 0.02;                         % -1.0 ... 4.0 ps
tLog  = logspace(log10(4), log10(3000), 71);      % 4 ps ... 3 ns
tLog  = tLog(2:end);
dec   = 10 .^ (floor(log10(tLog)) - 2);
tLog  = round(tLog ./ dec) .* dec;
t = round([tPre, tFine, tLog] * 1e6) / 1e6;
if any(diff(t) <= 0), error('Delay axis is not strictly increasing.'); end

g   = @(c, w) exp(-0.5 * ((wl - c) / w) .^ 2);
gsb = -(8.0 * g(505, 16) + 3.5 * g(475, 13));
SA  = [gsb - 4.0 * g(565, 18), 7.0 * g(690, 38)];
SB  = [gsb - 3.0 * g(580, 20), 6.0 * g(650, 32)];

t0 = chirpCurve(wl, 8, 500);                      % ps
k  = 1 ./ tau;
kr = 1 / tauRot;
% magic angle (M) and sum of (signal x anisotropy) (R)
[M, R]   = anisoRawSignal(t, t0, SA, SB, r0, k, kr, sig);
[Ml, Rl] = anisoRawSignal(t, t0 + tLate, SA, SB, r0, k, kr, sig);
bgPar  = 0.12 - 4.0 * g(pumpWL, 2.5);
bgPerp = 0.20 - 1.5 * g(pumpWL, 2.5);
noiseRel = (1 + 2 * ((wl - 600) / 150) .^ 2) .* (1 + 4 * g(pumpWL, 3));
par  = M + 2 * R + bgPar + 0.03 * noiseRel .* demoNoise(numel(wl), numel(t), 6101);
perp = (M - R) / G + bgPerp + 0.03 * noiseRel .* demoNoise(numel(wl), numel(t), 6102);
perpLate = (Ml - Rl) / G + bgPerp + 0.03 * noiseRel .* demoNoise(numel(wl), numel(t), 6103);
writeData(fullfile(outDir, 'demo_aniso_raw_par.csv'),  wl, t, par);
writeData(fullfile(outDir, 'demo_aniso_raw_perp.csv'), wl, t, perp);
writeData(fullfile(outDir, 'demo_aniso_raw_perp_late.csv'), wl, t, perpLate);

fid = fopen(fullfile(truthDir, 'demo_aniso_raw_truth.csv'), 'w');
fprintf(fid, 'wavelength_nm,t0_ps,background_par_mOD,background_perp_mOD,noise_mOD\n');
fprintf(fid, '%g,%.5f,%.4f,%.4f,%.4f\n', [wl, t0, bgPar, bgPerp, 0.03 * noiseRel].');
fclose(fid);
end

% ======================================================================
function [M, R] = anisoRawSignal(t, t0, SA, SB, r0, k, kr, sig)
% Magic-angle signal M and sum of (signal x anisotropy) R of the sample
% of writeAnisotropyRawDemo, with time zero t0 (one value per wavelength).
M = zeros(numel(t0), numel(t));
R = zeros(numel(t0), numel(t));
for i = 1:numel(t0)
    ti  = t(:) - t0(i);
    cA  = expIRF(ti, k(1), sig);
    cB  = k(1) / (k(2) - k(1)) * (cA - expIRF(ti, k(2), sig));
    cAr = expIRF(ti, k(1) + kr, sig);
    cBr = k(1) / (k(2) - k(1)) * (cAr - expIRF(ti, k(2) + kr, sig));
    M(i, :) = (sum(SA(i, :)) * cA + sum(SB(i, :)) * cB).';
    R(i, :) = ((SA(i, :) * r0(:)) * cAr + (SB(i, :) * r0(:)) * cBr).';
end
end

% ======================================================================
function writeSolvationDemo(outDir, truthDir)
% Sixth data set: time-resolved emission (intensity, not dA) of a dye in
% a polar solvent. The emission band shifts to the red while the solvent
% relaxes around the excited molecule (dynamic Stokes shift); the
% excited state itself decays with the fluorescence lifetime only.
% Not a particular molecule. No background, time zero at 0.
%
%   band      log-normal in wavenumber (width and asymmetry constant):
%             L(nu) = exp(-ln2 * (ln(1 + a) / gam)^2),
%             a = 2 * gam * (nu - nuPeak) / width   (0 where a <= -1)
%   peak      nuPeak(t) = nuInf + dNu * C(t),
%             C(t) = sum(aSolv .* exp(-t ./ tauSolv))
%   signal    photons per nm: N0 * exp(-t / tauFl) * L / area * dnu/dlambda,
%             so that the integral over the wavelength is the number of
%             excited molecules at every time (the dependence of the
%             radiative rate on the emission frequency is left out)
nuInf   = 18600;                  % cm-1, peak of the relaxed band (538 nm)
dNu     = 1800;                   % cm-1, shift from time zero to relaxed
aSolv   = [0.35 0.35 0.30];       % weights of the solvation times
tauSolv = [0.25 3.0 15];          % ps
width   = 3000;                   % cm-1, width parameter of the band
                                  % (its FWHM is width * sinh(gam) / gam
                                  % = 3031 cm-1)
gam     = -0.25;                  % asymmetry (tail to the red)
tauFl   = 4000;                   % ps, fluorescence lifetime
N0      = 10000;                  % integral over the wavelength at time zero
fwhm    = 0.15;                   % ps
sig     = fwhm / (2 * sqrt(2 * log(2)));
noise   = 0.5;                    % the same at all wavelengths

wl = (420:2:760).';               % nm
nu = 1e7 ./ wl;                   % cm-1
tLin = (-20:40) * 0.05;                           % -1 ... 2 ps
tLog = logspace(log10(2), log10(5000), 61);
tLog = tLog(2:end);
dec  = 10 .^ (floor(log10(tLog)) - 2);
tLog = round(tLog ./ dec) .* dec;
t = round([tLin, tLog] * 1e6) / 1e6;
if any(diff(t) <= 0), error('Delay axis is not strictly increasing.'); end

% area of the band (the same at all times), on a fine wavenumber grid
nuF  = (8000:2:32000).';
area = trapz(nuF, lognormBand(nuF, 20000, width, gam));
Cfun  = @(tt) exp(-tt(:) * (1 ./ tauSolv)) * aSolv(:);          % column
spec  = @(nuP) lognormBand(nu, nuP, width, gam) .* (1e7 ./ wl .^ 2) / area;

% the moving band on a fine time grid, IRF applied numerically; after
% 10 ps the IRF changes nothing and the model is taken directly
dtf = 0.002;
tf  = (-600:5500) * dtf;                          % -1.2 ... 11 ps
ker = exp(-0.5 * (((-150:150) * dtf) / sig) .^ 2);
ker = ker / sum(ker);
F = zeros(numel(wl), numel(tf));
for j = find(tf >= 0)
    F(:, j) = N0 * exp(-tf(j) / tauFl) * spec(nuInf + dNu * Cfun(tf(j)));
end
F(:, tf == 0) = 0.5 * F(:, tf == 0);
for i = 1:numel(wl)
    F(i, :) = conv(F(i, :), ker, 'same');
end
D = zeros(numel(wl), numel(t));
early = t <= 10;
D(:, early) = interp1(tf(:), F.', t(early).', 'linear', 0).';
for j = find(~early)
    D(:, j) = N0 * exp(-t(j) / tauFl) * spec(nuInf + dNu * Cfun(t(j)));
end
D = D + noise * demoNoise(numel(wl), numel(t), 7001);
writeData(fullfile(outDir, 'demo_solvation.csv'), wl, t, D);

% truth: the band at time zero and relaxed (for N0 molecules), and along
% the delays the population (with the IRF), the peak and C(t)
S0 = N0 * spec(nuInf + dNu);
S1 = N0 * spec(nuInf);
fid = fopen(fullfile(truthDir, 'demo_solvation_truth_spectra.csv'), 'w');
fprintf(fid, 'wavelength_nm,band_time_zero,band_relaxed\n');
fprintf(fid, '%g,%.4f,%.4f\n', [wl, S0, S1].');
fclose(fid);
pop = expIRF(t(:), 1 / tauFl, sig);
Ct  = Cfun(max(t, 0));
nuP = nuInf + dNu * Ct;
% maximum of the spectrum as the file holds it (per nm): a few nm to the
% blue of 1e7 / nuP, because of the factor dnu/dlambda
wlF = (420:0.01:760).';
lamMax = zeros(numel(t), 1);
for j = 1:numel(t)
    [~, k] = max(lognormBand(1e7 ./ wlF, nuP(j), width, gam) ./ wlF .^ 2);
    lamMax(j) = wlF(k);
end
Ct(t < 0) = NaN;  nuP(t < 0) = NaN;  lamMax(t < 0) = NaN;
fid = fopen(fullfile(truthDir, 'demo_solvation_truth_kinetics.csv'), 'w');
fprintf(fid, 'delay_ps,population,C_t,band_peak_cm-1,spectrum_maximum_nm\n');
fprintf(fid, '%.6g,%.6f,%.6f,%.1f,%.2f\n', [t(:), pop, Ct, nuP, lamMax].');
fclose(fid);
end

% ======================================================================
function L = lognormBand(nu, nuPeak, width, gam)
% Log-normal band shape with its maximum (1) at nuPeak, FWHM width and
% asymmetry gam (negative: tail to low wavenumbers).
a = 2 * gam * (nu - nuPeak) / width;
L = zeros(size(nu));
ok = a > -1;
L(ok) = exp(-log(2) * (log(1 + a(ok)) / gam) .^ 2);
end

% ======================================================================
function D = buildClean(wl, t, S, P, t0, sig, oscProf, oscT, tf)
% Signal without background and noise; t0 is the time zero of every row.
D = zeros(numel(wl), numel(t));
for i = 1:numel(wl)
    tt = t(:) - t0(i);
    C  = seqConc(tt, P.tau, 0, sig);
    y  = C * S(i, :).';
    for m = 1:size(oscProf, 2)
        y = y + oscProf(i, m) * interp1(tf(:), oscT(:, m), tt, 'linear', 0);
    end
    D(i, :) = y.';
end
end

% ======================================================================
function C = seqConc(t, tau, t0, sig)
% Populations of A -> B -> C -> ground state (unit start in A), with a
% Gaussian IRF of standard deviation sig applied.
k  = 1 ./ tau;
E  = zeros(numel(t), 3);
for j = 1:3
    E(:, j) = expIRF(t - t0, k(j), sig);
end
C = zeros(numel(t), 3);
C(:, 1) = E(:, 1);
C(:, 2) = k(1) / (k(2) - k(1)) * (E(:, 1) - E(:, 2));
C(:, 3) = k(1) * k(2) * ( ...
    E(:, 1) / ((k(2) - k(1)) * (k(3) - k(1))) + ...
    E(:, 2) / ((k(1) - k(2)) * (k(3) - k(2))) + ...
    E(:, 3) / ((k(1) - k(3)) * (k(2) - k(3))));
end

% ======================================================================
function y = expIRF(t, k, sig)
% exp(-k t) for t >= 0, convolved with a normalized Gaussian.
y = 0.5 * exp(-k * t + 0.5 * (k * sig) ^ 2) .* ...
    erfc((k * sig ^ 2 - t) / (sig * sqrt(2)));
y(~isfinite(y)) = 0;
end

% ======================================================================
function t0 = chirpCurve(wl, Lmm, wlRef)
% Time zero (ps) from the group delay of Lmm of fused silica; zero at
% wlRef. Wavelengths in nm. Sellmeier coefficients: Malitson (1965).
ng  = @(w) groupIndex(w / 1000);
c0  = 0.299792458;                 % mm / ps
t0  = -(Lmm / c0) * (ng(wl) - ng(wlRef));
end

function ng = groupIndex(lam)
% Group index of fused silica, lam in micrometres.
B  = [0.6961663 0.4079426 0.8974794];
Cc = [0.0684043 0.1162414 9.896161] .^ 2;
n  = @(x) sqrt(1 + B(1) * x.^2 ./ (x.^2 - Cc(1)) + ...
                   B(2) * x.^2 ./ (x.^2 - Cc(2)) + ...
                   B(3) * x.^2 ./ (x.^2 - Cc(3)));
h  = 1e-4;
ng = n(lam) - lam .* (n(lam + h) - n(lam - h)) / (2 * h);
end

% ======================================================================
function Z = demoNoise(nRow, nCol, seed)
% Standard normal numbers, nRow x nCol. Park-Miller generator with the
% Box-Muller transform: identical in MATLAB and Octave.
n   = nRow * nCol;
nH  = ceil(n / 2);
mod31 = 2147483647;
x   = seed;
for k = 1:20                       % warm-up
    x = mod(16807 * x, mod31);
end
u = zeros(2 * nH, 1);
for k = 1:2 * nH
    x = mod(16807 * x, mod31);
    u(k) = x / mod31;
end
r  = sqrt(-2 * log(u(1:nH)));
ph = 2 * pi * u(nH + 1:end);
z  = [r .* cos(ph); r .* sin(ph)];
Z  = reshape(z(1:n), nRow, nCol);
end

% ======================================================================
function writeData(fname, wl, t, D)
% First row: 0, delays. Following rows: wavelength, signal.
fid = fopen(fname, 'w');
if fid < 0, error('Cannot write %s', fname); end
fprintf(fid, '0');
fprintf(fid, ',%.6g', t);
fprintf(fid, '\n');
fmt = ['%g', repmat(',%.4f', 1, numel(t)), '\n'];
fprintf(fid, fmt, [wl, D].');
fclose(fid);
end
