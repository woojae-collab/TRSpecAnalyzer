TRSpecAnalyzer demo data
========================

Synthetic transient-absorption data with known answers, for a first run of
TRSpecAnalyzer.m. Everything in this folder is simulated; it contains no
measurement.


Files
-----
demo_TA.csv                     the data set of the walkthrough
                                (151 wavelengths x 331 delays)
scans/demo_scan_01..05.csv      five repeated scans, for Load & Average
truth/demo_truth_spectra.csv    true spectra of A, B and C, oscillation
                                profiles, time zero of every wavelength,
                                background, noise level
truth/demo_truth_kinetics.csv   true populations
demo_target.csv                 second data set, for Target Analysis
                                (151 wavelengths x 121 delays)
truth/demo_target_truth_*.csv   its true spectra and populations
demo_equilibrium.csv            third data set, for Target Analysis of a
                                step that runs in both directions
                                (151 wavelengths x 121 delays)
truth/demo_equilibrium_truth_*.csv   its true spectra and populations
demo_solution.csv               fourth data set: a solute in cyclohexane,
                                measured with short pulses
                                (76 wavelengths x 1086 delays)
demo_solvent.csv                cyclohexane alone, measured the same way
truth/demo_solvent_truth_*.csv  true spectra, profiles of the oscillations
                                and of the coherent artifact; list of the
                                modes
demo_solution_raw.csv           the fourth data set as measured: with
demo_solvent_raw.csv            background and chirp (same size)
truth/demo_solvent_raw_truth.csv   time zero and the two backgrounds of
                                every wavelength
demo_aniso_par.csv              fifth data set: one sample measured with
demo_aniso_perp.csv             parallel and with perpendicular
                                polarization (151 wavelengths x 121 delays)
truth/demo_aniso_truth_*.csv    its true spectra and populations
demo_aniso_raw_par.csv          the fifth data set as measured: with
demo_aniso_raw_perp.csv         background, pump scatter and chirp
                                (151 wavelengths x 331 delays)
demo_aniso_raw_perp_late.csv    the perpendicular measurement once more,
                                its time zero 0.05 ps late
truth/demo_aniso_raw_truth.csv  time zero, backgrounds and noise of
                                every wavelength
demo_solvation.csv              sixth data set: time-resolved emission
                                of a dye whose band shifts while the
                                solvent relaxes (171 wavelengths x
                                121 delays; intensity)
truth/demo_solvation_truth_*.csv   the band at time zero and relaxed;
                                population, C(t) and peak position
demo_preview.png                picture of demo_TA.csv and of its truth
demo_target_preview.png         the same for demo_target.csv
demo_equilibrium_preview.png    the same for demo_equilibrium.csv
demo_solvent_preview.png        picture of the cyclohexane pair
demo_aniso_preview.png          picture of the polarization pair
demo_solvation_preview.png      picture of demo_solvation.csv
make_demo_data.m                the script that wrote all files

File layout: first row = delays (ps), first column = wavelengths (nm),
the rest = signal (mOD). In the main window choose
    Time unit = ps        Spectral axis = Wavelength (nm)


What is in the data
-------------------
Kinetics      A -> B -> C -> ground state
              tau = 0.60 ps, 15 ps, 700 ps
IRF           Gaussian, FWHM 0.060 ps
Chirp         time zero runs from -0.26 ps at 450 nm to +0.56 ps at
              750 nm and is 0 at 500 nm (group delay of 8 mm fused silica)
Oscillations  150 cm-1, damping 1.0 ps: strongest at 558 and 602 nm,
              node at 580 nm
              250 cm-1, damping 0.8 ps: strongest at 625 and 685 nm,
              node at 655 nm
Background    offset of +0.12 mOD and pump scatter (-4 mOD) at 500 nm,
              at all delays
Noise         0.03 mOD at 600 nm, 0.09 mOD at 450 and 750 nm, five times
              larger at 500 nm
There is no coherent artifact. (The coherent-artifact terms of Global
Analysis, Target Analysis and Auto Analysis have nothing to fit here.)


Walkthrough
-----------
The expected numbers below were obtained by running the calculation code
of TRSpecAnalyzer on these files, step by step as written.

 1. Load...
      demo_TA.csv -> Load

 2. Background Corr...
      # of initial delay points to average = 10 -> Apply & Close
      (The offset and the pump scatter line disappear.)

 3. Chirp Corr...
      Auto-place (8 points) -> Fit & Apply
      The fitted curve follows the true time zero within 0.01 ps.
      The green ridge jumps below 465 nm and around 600 nm, where the
      signal at time zero is nearly zero. The eight points do not fall
      there. If you place points by hand, avoid these two regions.

 4. Crop...
      t range -1.5 to 2800 -> Apply Crop
      (Removes the empty delays that the chirp correction leaves at both
      ends.)

 5. Mask lambda...
      From 494, To 506, Mode "Set to NaN (mask out)" -> Add region
      (Pump scatter: the constant part is gone after step 2, the extra
      noise is not.)

 6. SVD...
      Three singular values stand far above the rest. The fourth and
      fifth, just above the noise, are the two oscillations.

 7. Global Analysis...
      Components = 3, start values 1 / 10 / 100, untick "Fixed" for
      IRF t_0 and IRF FWHM -> Run Fit
      Expected   tau = 0.60, 15.0, 701 ps      (true 0.60, 15, 700)
                 IRF FWHM = 0.060 ps           (true 0.060)
                 IRF t_0  = 0.000 ps
                 RMS = 0.065 mOD
      The EADS agree with the true spectra of A, B and C within 0.1 mOD.
      The RMS is twice the noise because the residual still holds the
      oscillations (step 9).

 8. Target Analysis...
      Model "Sequential, all decay:  A -> B -> C -> GS", "tau from GA"
      -> Run fit
      Expected: the same time constants and spectra as in step 7.
      Do step 7 first. The target window takes t_0 and FWHM from the
      global fit and keeps them fixed. Without step 7 it uses
      FWHM = 0.15, too wide for these data: tau = 0.56, 15.3, 697 ps and
      RMS 0.16.
      If t_0 and FWHM are fitted in the target window instead (tested
      with start values 1 / 10 / 100 and FWHM 0.15), the answer is the
      same as in step 7. The first search stops at a poorer solution
      with an FWHM far below the delay step (RMS 0.12); the
      program notices this and repeats the fit in two steps. The result
      box then says "The direct fit lost the IRF width; this is the
      result of a fit in two steps (IRF fixed, then free)."

 9. Coherence...
      Method "Use GA residual", Time window 0.15 to 4 -> Run Coherence
      Expected: two peaks, at 150 and 250 cm-1.
      150 cm-1: maxima near 556 and 606 nm, node at 580 nm
      250 cm-1: maxima near 620 and 682 nm, node at 656 nm
      (The 250 cm-1 mode is weak. Below 500 nm, where the noise is more
      than twice as large, single wavelengths of its profile come as
      high as these maxima, for example 474 nm: that is noise, not the
      mode.)

10. LDA...   (optional)
      Broad bands around 0.5, 15 and 600 to 700 ps. Their exact position
      and width depend on the regularization.


Auto Analysis (steps 2 to 9 in one run)
---------------------------------------
Load demo_TA.csv, then  Auto Analysis... -> Run.
The Report tab should show (numbers from the calculation code of
TRSpecAnalyzer run in GNU Octave):

  Background   mean of the first 43 delays (t <= -0.36 ps)
  Mask         none. 500 ... 502 nm is listed as "noisier than the
               neighbors, but kept" (noise up to 4.7 times, signal at
               least 26 times the noise)
  Chirp        time zero from -0.27 ps at 450 nm to +0.56 ps at 750 nm
               (within 0.01 ps of the true curve)
  Crop         t = -1.7 ... 2730 ps
  Artifact     coherent-artifact terms not needed (three terms would
               lower the RMS by 3 %, less than the 15 % asked for)
  Components   3:  tau = 0.604, 15.0, 700 ps;  IRF FWHM 0.060 ps;
               RMS 0.074
  Oscillations 150 and 249 cm-1, nodes near 580 and 654 nm

"Apply corrections to the data" writes the corrections into the main
window. "Send model to Global Analysis" then opens the Global Analysis
window with these values; "Run Fit" there gives the same result.

The numbers differ a little from the walkthrough because the corrections
differ: 43 instead of 10 background delays, another crop, and no mask.
The walkthrough masks 494 to 506 nm by hand. The automatic analysis
keeps these channels, because their signal is far above their noise.
With  Mask rule = "Every noisy channel"  it masks 494 ... 506 nm as
well and gives tau = 0.605, 15.0, 701 ps and RMS 0.065.

"Components used" shows "3  (automatic choice)" after the run. Tick
"Fit every number up to the largest" before Run and it lists 1 to 5;
picking another number switches the report and the plots to that fit
(nothing is fitted again). With 4, for example, the report shows an
extra time constant of 0.013 ps, a fifth of the IRF width, for an RMS
only 1.5 % lower: the reason the automatic choice stopped at 3.

For demo_target.csv the report proposes no correction and finds
tau = 2.00, 40.0, 799 ps. It does not find the branching: that takes the
Target Analysis steps further below.

For demo_equilibrium.csv the report proposes no correction and finds
tau = 3.71, 134.5, 1998 ps. These are the lifetimes the data show, not
the four steps of the scheme (see the equilibrium demo below).

For demo_solution.csv the report finds tau = 0.413, 80.1 ps and lists
oscillations at 215, 590, 801 and 1028 cm-1 with their nodes. It cannot
tell which of them belong to the solvent (see the solvent demo below).
It also proposes a chirp correction of 0.012 ps from end to end. The
data have no chirp: the coherent artifact changes its shape along the
spectrum and moves the apparent rise. Do not apply the corrections to
this file.

For demo_solvation.csv (loaded as Intensity) the report proposes no
correction and chooses 4 components, tau = 0.256, 2.72, 13.1, 4000 ps,
as step 5 of the solvation demo below. It cannot determine a time zero
along the spectrum: the late rise at the red wavelengths is the shift
of the band. The data need no chirp correction.


Load & Average (optional)
-------------------------
Load... and select the five files in scans/ together.

- Scan 5 has only 65 % of the signal (a degraded sample). Its Max |dA|
  is 10.9 against 13.7 to 14.1 for the others. Untick "Include".
- Scan 3 is 0.06 ps late. Set Ref lambda = 690, mark scan 1 as Ref and
  press "Auto-align to Ref": the Shift of scan 3 becomes -0.06.
  At the default Ref lambda (600 nm) the signal has no sharp rise and
  the shift is not found.
- "Load Average of Checked", then continue with step 2. The result is
  the same as for demo_TA.csv (tau = 0.60, 15.0, 700 ps).


Target analysis demo (demo_target.csv)
--------------------------------------
A second data set for a model that a sequential fit cannot describe.
It is already corrected: no background, no chirp, time zero at 0.

True model    A -> B (3 ps),  A -> C (6 ps),  B -> GS (40 ps),
              C -> GS (800 ps)
              A decays in 2.0 ps; 67 % goes to B, 33 % to C
IRF           Gaussian, FWHM 0.10 ps
Spectra       from 450 to 475 nm all three states show the same
              ground-state bleach and nothing else
Noise         0.03 mOD at 600 nm, 0.09 mOD at 450 and 750 nm

 1. Load...
      demo_target.csv -> Load

 2. Global Analysis...
      Components = 3, start values 1 / 10 / 100, untick "Fixed" for
      IRF t_0 and IRF FWHM -> Run Fit
      Expected   tau = 2.00, 40.0, 799 ps;  FWHM 0.100 ps;  RMS 0.052
      EADS 1 is the spectrum of A. EADS 2 is a mixture,
      0.67 x B + 0.33 x C, and EADS 3 is 0.32 x C: the sequential model
      knows nothing about the branching.

 3. Target Analysis...
      In the "Scheme (text)" tab type
          A -> B : k1
          A -> C : k2
          B -> GS : k3
          C -> GS : k4
      -> Apply scheme. Set the start values of k1 to k4 to 3, 5, 30, 500.
      Leave t_0 and FWHM fixed (they come from step 2). -> Run fit
      Expected: RMS 0.052 as in step 2, but k1 and k2 are marked
      "[not determined]". The data fix only their combination (A decays
      in 2.0 ps). How A divides between B and C, and with it the size
      of the spectra of B and C, depends on the start values
      (here: 71 % to B).

 4. Add two lines to the scheme
          equal: B = A 450 475
          equal: C = A 450 475
      -> Apply scheme -> Run fit
      Expected   k1 = 2.99, k2 = 6.02, k3 = 40.0, k4 = 799 ps
                 (true 3, 6, 40, 800);  67 % of A goes to B, 33 % to C
                 RMS 0.052
      The SAS agree with the true spectra within 0.2 mOD
      (truth/demo_target_truth_spectra.csv). The Diagram tab shows the
      scheme with the fitted values.
      The two lines say that between 450 and 475 nm B and C have the
      spectrum of A, the ground-state bleach. The bleach then counts
      the excited molecules, and that fixes the branching.

 5. Test of another model (optional)
      Replace the four steps by   A -> B -> C -> GS   and keep the two
      "equal" lines. Start values 3, 5, 30. -> Apply scheme -> Run fit
      Expected: RMS 0.27, five times higher. Without the two lines the
      sequential model fits as well as the true one (RMS 0.052); with
      them it is ruled out.

 6. Fit trace on a band (optional)
      Fit trace... -> lambda = 450, tick "Band, to:" and enter 475,
      keep "Integral" and Components = 3, untick "Fixed" for IRF t_0
      and IRF FWHM -> Run Fit
      Expected   tau = 2.5, 39.6, 794 ps      (true 2, 40, 800)
                 A_2 : A_3 = 68 : 32
      Between 450 and 475 nm every state shows the same bleach, so the
      band counts the excited molecules: the two slow amplitudes divide
      as A divides between B and C (67 : 33), without a model.
      One wavelength (462 nm, "Band, to:" unticked) gives the same
      values with uncertainties three times larger (39.7 +/- 1.0 ps
      against 39.6 +/- 0.3 ps).

Two things to know

- Start values matter in step 4. The two "equal" lines leave a second
  solution that fits equally well (RMS 0.052): A lives 40 ps and B only
  2 ps. The fit ends there when it starts from the values the window
  proposes (5, 34.7, 240, 1670). You recognize it by the spectrum of B,
  which comes out around 180 mOD, almost twenty times the data. Start
  values of the right order (a few ps, a few ps, tens of ps, hundreds
  of ps) lead to the true solution; 1/10/100/1000, 2/4/20/300 and
  5/10/60/1500 were tried as well.
- Do not use "tau from GA" in step 3. The global fit has three time
  constants, the scheme four rate constants, and the button would put
  40 ps on k2.


Equilibrium demo (demo_equilibrium.csv)
---------------------------------------
A third data set, with a step that runs in both directions. It is
already corrected: no background, no chirp, time zero at 0.

True model    A <-> B (5 ps forward, 15 ps back),  B -> C (100 ps),
              C -> GS (2000 ps)
              equilibrium constant K = kf/kb = 3
IRF           Gaussian, FWHM 0.10 ps
Spectra       B has no signal from 700 to 750 nm; A and C have
Noise         0.03 mOD at 600 nm, 0.09 mOD at 450 and 750 nm

The data do not show 5, 15 and 100 ps. A and B exchange faster than B
reacts on, and the three lifetimes in the data are
    3.71 ps   A and B reach their equilibrium (about 1/(kf + kb))
    135 ps    the pair, 26 % A and 74 % B, empties into C
              (about 100 ps / 0.74)
    2000 ps   C decays

 1. Load...
      demo_equilibrium.csv -> Load

 2. Global Analysis...
      Components = 3, start values 1 / 10 / 100, untick "Fixed" for
      IRF t_0 and IRF FWHM -> Run Fit
      Expected   tau = 3.71, 134.5, 1998 ps;  FWHM 0.100 ps;  RMS 0.053
      EADS 1 is the spectrum of A and EADS 3 that of C. EADS 2 is not
      B. It is 0.26 x A + 0.74 x B, the mixture at equilibrium, and
      still shows the band of A at 710 nm at a quarter of its height.

 3. Target Analysis...
      In the "Scheme (text)" tab type
          A <-> B : kf, kb
          B -> C : k2
          C -> GS : k3
      -> Apply scheme. Set the start values of kf, kb, k2, k3 to
      3, 10, 50, 1000. Leave t_0 and FWHM fixed (they come from
      step 2). -> Run fit
      Expected: RMS 0.053 as in step 2, but kf, kb and k2 are marked
      "[not determined]" (here 4.0, 50 and 124 ps; other start values
      give other numbers). "Observable lifetimes" lists 3.712, 134.5
      and 1998 ps: the data hold three numbers, the scheme asks for
      four.

 4. Add one line to the scheme
          zero: B 700 750
      -> Apply scheme -> Run fit
      Expected   kf = 4.99, kb = 15.0, k2 = 100.0, k3 = 1998 ps
                 (true 5, 15, 100, 2000);  RMS 0.053
                 K = 15.0 / 4.99 = 3.0
                 (-2.7 kJ/mol or -28 meV at 295 K)
      The SAS agree with the true spectra within 0.1 mOD
      (truth/demo_equilibrium_truth_spectra.csv). The Diagram tab
      shows the scheme with the fitted values.
      The line says that B has no signal from 700 to 750 nm. The
      signal there then belongs to A and C alone, and the share of A
      that is left after the fast step (26 %) gives kb against kf.
      "Branching" lists  B: -> A 86.9 %, -> C 13.1 %. A molecule in B
      goes back to A about seven times before it moves on to C.

 5. Test of another model (optional)
      Replace the three steps by   A -> B -> C -> GS   and keep the
      "zero" line. Start values 3, 10, 100. -> Apply scheme -> Run fit
      Expected: RMS 0.21, four times higher (tau = 5.6, 102, 2030 ps).
      Without the "zero" line the sequential model fits as well as the
      true one (RMS 0.053); with it, it is ruled out.

 6. The same answer from the equilibrium constant (optional)
      If K is known from another experiment, the back step can be tied
      to the forward step instead. Remove the "zero" line and type
          A <-> B : kf, kf/3
          B -> C : k2
          C -> GS : k3
      Start values 3, 50, 1000. -> Apply scheme -> Run fit
      Expected   kf = 5.00, k2 = 99.9, k3 = 1998 ps;  RMS 0.053
      The back step is listed under "Steps with a tied rate" with
      15.0 ps. (The result box writes kf/3 as 0.3333333333*kf.)

Things to know

- Start k2 shorter than k3. With the "zero" line a second solution
  fits equally well (RMS 0.053): kf = 4.99, kb = 14.5, k2 = 1486,
  k3 = 134.5 ps. There C forms slowly and decays fast, so that little
  of it is ever present, and its spectrum comes out near 97 mOD, ten
  times the data. The fit ends there when the start value of k2 is
  longer than that of k3 (seen with 3/10/1000/50 and 5/15/1670/240).
  Start values with k2 below k3 gave the true solution: 3/10/50/1000,
  1/10/100/1000, 2/5/30/500, 10/10/100/1000, 20/5/200/3000 and the
  values the window proposes (5, 34.7, 240, 1670).
- Do not use "tau from GA" in steps 3 and 4. The global fit has three
  time constants, the scheme four rate constants; the button puts
  134.5 ps on kb and 1998 ps on k2, and the fit ends in the second
  solution.
- A "zero" range has to be true. With  zero: B 650 750  (B has up to
  1.1 mOD there) the fit gives kf = 5.03, kb = 14.7, k2 = 97.0 ps and
  RMS 0.064 instead of 0.053. Compare the RMS with step 3: a
  constraint that is true costs almost nothing.
- An unchanged RMS does not show that a "zero" line is true. If B has
  a signal of the same shape as A in the range, the fit stays as good
  and only the rate constants come out wrong. In a test data set where
  B was given half of the 710 nm band of A, the same line
  zero: B 700 750  gave kf = 10.0, kb = 6.4, k2 = 50 ps at RMS 0.053.
  That B has no signal in the range has to be known from outside the
  fit.
- If the spectra of A and B overlap everywhere, no "zero" line is
  possible and one number has to come from elsewhere. In the same test
  data three things gave the true rate constants: the equilibrium
  constant (kf/3 as in step 6), one rate constant held with "Fixed"
  (k2 = 100 ps gave kf = 4.99, kb = 15.0 ps), and a known ratio of
  the two spectra in a range (equal: B = 0.5*A 700 750). Equal bleach
  (equal: B = A 450 470) does not help, because the sum of A and B
  does not change while they exchange.
- Do step 2 first, as for the other data sets. Without it the target
  window uses FWHM = 0.15: kf = 4.89, kb = 14.9, k2 = 100.8,
  k3 = 1993 ps and RMS 0.084.


Solvent demo (demo_solution.csv and demo_solvent.csv)
-----------------------------------------------------
A fourth data set, in two files: a solute in cyclohexane, and
cyclohexane alone measured the same way. The pulses are short (20 fs),
so the solvent oscillates as well and a coherent artifact stands at
time zero. Both files are already corrected: no background, no chirp,
time zero at 0. The delays run in steps of 0.005 ps up to 5 ps and in
wider steps from there to 300 ps.

Solute        A -> B -> ground state,  tau = 0.40 ps, 80 ps
              oscillations at 215 cm-1 (damping 1.2 ps, node at 560 nm)
              and 590 cm-1 (damping 0.9 ps, node at 650 nm)
              no signal above 725 nm
Cyclohexane   oscillations at its Raman lines 384, 426, 801, 1028,
              1158, 1266 and 1444 cm-1, with amplitudes in the ratio of
              the Raman intensities (2 : 3 : 95 : 15 : 6 : 14 : 12);
              one profile for all of them, with a node at 565 nm
              a coherent artifact at time zero
In the solution   the cyclohexane signal is 0.70 times as large as in
              demo_solvent.csv and comes 0.003 ps later
IRF           Gaussian, FWHM 0.020 ps
Noise         0.02 mOD at 600 nm, 0.06 mOD at 450 and 750 nm, in each
              file

 1. Load...
      demo_solution.csv -> Load

 2. Global Analysis...
      Components = 2, start values 1 / 10, untick "Fixed" for IRF t_0
      and IRF FWHM -> Run Fit
      Expected   tau = 0.41, 80.0 ps      (true 0.40, 80)
                 IRF FWHM = 0.024 ps      (true 0.020)
                 RMS = 0.118 mOD
      The RMS is several times the noise: the residual holds the
      oscillations and the coherent artifact. The artifact also widens
      the fitted IRF.

 3. Coherence...
      Method "Use GA residual", Time window 0.15 to 5 -> Run Coherence
      Expected: one line at 801 cm-1 far above everything else (maxima
      near 486 and 642 nm, node at 566 nm). In the spectrum summed over
      all wavelengths the next peaks are 215 cm-1 (0.12 of the 801
      peak), 590 cm-1 (0.02) and 1028 cm-1 (0.004).
      Nothing in this result says which peaks belong to the solute.

 4. Solvent coherence...   (button in the Coherence window)
      Load solvent data... -> demo_solvent.csv
      The band is set to the strongest mode of the solvent file,
      795 to 808 cm-1. -> Estimate factor + shift
      Expected   A = 0.700, time shift 0.0030 ps   (true 0.70, 0.003)
                 0.3 % of the power in the band is left
      -> Apply. The Coherence window runs again.
      Expected: no peak at 801 and 1028 cm-1; 215 and 590 cm-1 remain.
      215 cm-1: maxima near 538 and 582 nm, node at 558 nm
      590 cm-1: maxima near 630 and 678 nm, node at 650 nm
      The map is normalized to its largest value, so the 215 cm-1 line
      now fills the color scale.
      The data of the main window are not changed. The solvent is taken
      out of the residual only, before the Fourier transform.

      Before and after together (the bar above the plots):
      - The four curve plots show the curve from before the subtraction
        in gray behind the result ("Before" in gray).
      - Map shows: Before subtraction -> the map as it was, in the
        units of the map after the subtraction. Its largest value is
        2.25, at 801 cm-1.
      - Side by side... -> a window with both maps on one color scale
        and the two power spectra laid over each other.
      Color range (same bar; it is there without a solvent file as
      well): type 0.15 as the maximum -> the 590 cm-1 line, 0.14 of the
      215 cm-1 line, fills the color scale. Auto returns to the range
      of the map.

Things to know

- The time shift matters. "Estimate factor" alone (time shift 0) gives
  A = 0.63 and leaves 19 % of the power in the band: 0.003 ps is 7 % of
  a period of the 801 cm-1 mode.
- The noise of the solvent file comes with the subtraction. The floor
  of the spectrum rises by about one half (0.70 squared).
- Use a strong line of the solvent, with no mode of the solute in the
  band. A weak line gives a factor that is too small: the estimate is a
  least-squares factor, and the noise of the solvent file counts as
  solvent signal in it. The 1028 cm-1 line (band 1010 to 1045) gives
  A = 0.42 here.

The data without the solvent in the main window (after "Apply" in step 4)

 4c. View (main window) -> Solvent coherence subtracted
      "View" is on once "Apply" has been pressed in step 4. It offers
        Data | Solvent coherence subtracted | Solvent signal subtracted
      The Coherence window and the other analysis windows close; their
      setup and results are kept for the view "Data" and come back
      with it. In the main window:
        Signal: ... | solvent coherence subtracted
                      (A = 0.7002, 0.15 ... 5 ps)
        title of the 2D map: ... (solvent coh. subtracted)
      The 2D map, the kinetics, the export and every analysis now use
      the data without the oscillating part of the solvent, in the time
      window of the Coherence analysis. Before 0.15 ps and after 5 ps
      the data are as they were.

 5c. Global Analysis as in step 2
      Expected   tau = 0.41, 80.0 ps;  IRF FWHM = 0.024 ps;  RMS 0.091
      (0.118 in step 2: the solvent oscillation is no longer in the
      residual. The coherent artifact is still there.)

 6c. Coherence as in step 3
      Expected: the result of step 4. Beside the button "Solvent
      coherence...": out of the data in this view: A = 0.7002, shift
      0.003004 ps. Before and after are compared as in step 4;
      "before" is the residual with what the view took out of the data
      put back. "Solvent coherence..." does not open here: factor and
      shift are set in the view "Data".

 7c. View -> Solvent signal subtracted
        Signal: ... | solvent signal subtracted (A = 0.7002, all delays)
      The whole solvent signal is taken out, at all delays, with the
      factor and the shift of step 4: the oscillation next to time zero
      and the coherent artifact go as well.
      Global Analysis as in step 2:
        tau = 0.40, 80.1 ps;  IRF FWHM = 0.020 ps;  RMS 0.066

 8c. View -> Data
      The data as loaded. The fit of step 2 and the Coherence analysis
      of step 4 are back; the fits of the other two views come back
      when these are chosen again.
      The views end when nothing is subtracted in the Coherence window
      any more ("Remove subtraction" in the dialog, or "Reset" of the
      Coherence window in the view "Data").

The same with Subtract data (instead of step 4)

 4b. Subtract data...   (main window)
      demo_solvent.csv is already there if step 4 was done; otherwise
      Load reference... -> demo_solvent.csv
      Estimate in lambda from 725 to 750, t from -0.1 to 0.1
      -> Estimate scale + shift
      Expected   scale 0.703, time shift 0.0032 ps
      -> Apply
      The line at time zero above 700 nm, the coherent artifact,
      disappears from the map.
      In this range the solute has no signal and the solvent signal
      (its artifact) is far above the noise. Both are needed:
      - 725 to 750 nm at all delays gives a scale of 0.61. Most of
        these points hold little but noise, and noise in the reference
        lowers a least-squares factor.
      - The whole data (450 to 750 nm) give -0.53: the solute signal
        decides the result.

 5b. Global Analysis as in step 2
      Expected   tau = 0.40, 80.1 ps;  IRF FWHM = 0.020 ps;  RMS 0.066

 6b. Coherence as in step 3
      Expected: 215 and 590 cm-1 only.
      ("Solvent coherence" of step 4 is switched off by step 4b: its
      factor belonged to the data before the subtraction.)

Which of the three

- Subtract data changes the data. The coherent artifact goes as well,
  and the global fit then finds the true IRF. It needs a range in which
  only the solvent gives a signal, and it adds the noise of the solvent
  file to the data at all delays.
- Solvent coherence with Apply leaves the data alone and removes the
  solvent from the oscillations only. It needs a solvent mode away from
  the modes of the solute, but no range free of solute signal. The
  coherent artifact stays in the data.
- The views of the main window use the factor and the shift found
  there. "Solvent coherence subtracted" takes the oscillating part of
  the solvent signal out of the data, in the time window of the
  Coherence analysis: the map, the kinetics and every analysis are
  free of the solvent oscillation there; the coherent artifact stays.
  "Solvent signal subtracted" takes the whole solvent signal out at
  all delays, as Subtract data does, with a factor and a shift found
  on a solvent line instead of a range free of solute signal. The
  noise of the solvent file is then added at all delays, and the
  solvent file needs the background correction of the data (without
  one, its baseline is subtracted as well). Setup and results of the
  analysis windows are kept for each view.


The pair as measured (demo_solution_raw.csv and demo_solvent_raw.csv)

The same solution and solvent before any correction, with other random
numbers for the noise:
Background    solution: 0.25 mOD, with a band of -1.5 mOD at 520 nm
              solvent:  0.10 mOD, with a band of -0.6 mOD at 520 nm
Chirp         time zero runs from -0.091 ps at 450 nm to +0.062 ps at
              750 nm (0 at 565 nm), the same in both files
truth/demo_solvent_raw_truth.csv lists time zero and the two backgrounds
of every wavelength.

 1. Load...
      demo_solution_raw.csv -> Load

 2. Background Corr...
      # of initial delay points to average = 10 -> Apply & Close

 3. Chirp Corr...
      Auto-place (8 points) -> Fit & Apply
      Expected: fit RMS 0.005 ps. The curve lies within 0.015 ps of the
      true time zero at all wavelengths (0.004 ps rms).

 4. Global Analysis as in step 2 above
      Expected   tau = 0.41, 79.8 ps;  IRF FWHM = 0.023 ps;  RMS 0.115

 5. Coherence as in step 3 above
      Expected: the same picture, the 801 cm-1 line far above the rest.

 6. Solvent coherence...
      Load solvent data... -> demo_solvent_raw.csv
      The file is loaded as measured. The program gives it the
      corrections of the data: its own background (the mean of its
      first 10 delays), the same chirp curve, the same masks and crop.
      -> Estimate factor + shift
      Expected   A = 0.699, time shift 0.0030 ps   (true 0.70, 0.003)
                 0.3 % of the power in the band is left
      -> Apply: 0.2 % of the 801 cm-1 line is left in the map.
      Then View -> Solvent coherence subtracted: title of the 2D map
        ... (BG(N=10), chirp(Sellmeier), solvent coh. subtracted)
      and the global fit of step 4 then gives RMS 0.086.
      View -> Solvent signal subtracted:
        ... (BG(N=10), chirp(Sellmeier), solvent signal subtracted)
        tau = 0.40, 80.0 ps;  IRF FWHM = 0.020 ps;  RMS 0.068

Things to know

- The solvent file must be the measurement as it came from the
  instrument, with the same wavelengths and delays as the data. Do not
  correct it by hand first: it would be corrected twice.
- The 801 cm-1 line is as strong here as in the corrected pair above
  (within 1 %). The chirp correction moves every kinetic trace along
  the delay axis; the values between the measured delays come from a
  cubic spline wherever the delays are equally spaced (here up to
  5 ps), which keeps an oscillation with 8 delays in a period at its
  height. Where the delay step changes, and between unevenly spaced
  delays, the interpolation is linear: a spline would overshoot there.
  Linear interpolation damps a fast oscillation (8 delays in a period:
  9 % of its power on average), so keep the delays equally spaced as
  far as the oscillations are to be followed.


Anisotropy demo (demo_aniso_par.csv and demo_aniso_perp.csv)
------------------------------------------------------------
A fifth data set, in two files: one sample measured with the probe
polarized parallel and perpendicular to the pump. Both files are
already corrected: no background, no chirp, time zero at 0.

Kinetics      A -> B -> ground state,  tau = 3 ps, 800 ps
Rotation      50 ps: the anisotropy decays with this time constant, in
              A and in B
Anisotropy    at time zero 0.38 for the bleach and the stimulated
              emission, -0.12 for the excited-state absorption
G             the perpendicular file is too small by the factor 1.15,
              as if that polarization were detected less well
IRF           Gaussian, FWHM 0.10 ps
Noise         0.03 mOD at 600 nm, 0.09 mOD at 450 and 750 nm, in each
              file

 1. Load...
      demo_aniso_par.csv -> Load

 2. Anisotropy...
      Load perpendicular... -> demo_aniso_perp.csv
      G: choose "Tail matching" and keep the range proposed, 545 to
      5000 ps (the last decade of the delays). The rotation is over
      there, and the two files have to agree.
      Expected   G = 1.15       (true 1.15)
      The plot on the right shows the G of every single wavelength in
      this range, with its uncertainty, and the one G in use as a
      dashed line.
      Expected   the points scatter about 1.15, from 1.06 to 1.26;
                 title: "scatter: 0.96 x what the noise explains"
      A value near 1 says that one G holds for the whole spectrum.
      Grey points are wavelengths whose G the noise leaves open (the
      blue edge, 600 nm, above 710 nm). The plot is a check only:
      the calculation uses the one G.
      -> Apply
      The main window shows "View: Anisotropy r"; the color scale runs
      from -0.5 to 0.5.
      Expected at early delays   r = 0.38 from 470 to 570 nm
                                 r = -0.12 from 630 to 740 nm
      and r close to 0 after 200 ps.
      r is left empty before time zero and wherever the magic-angle
      signal is below 3 times its noise. Around 600 nm the bleach and
      the absorption cancel: the magic-angle signal goes through zero
      and r runs off to large values on both sides. That belongs to the
      definition of r and says nothing about the sample.
      With G = 1 (the value the window opens with) r levels off near
      0.05 instead of 0.

 3. Fit trace...   (View: Anisotropy r)
      lambda = 506, Components = 1 -> Run Fit
      Expected   tau = 48 ps, A_1 = 0.385       (true 50 ps, 0.38)
      lambda = 690:
      Expected   tau = 48 ps, A_1 = -0.122      (true 50 ps, -0.12)

 4. View: Magic angle, then Global Analysis...
      Components = 2, start values 1 / 10, untick "Fixed" for IRF t_0
      and IRF FWHM -> Run Fit
      Expected   tau = 3.00, 800 ps   (true 3, 800)
                 IRF FWHM = 0.100 ps;  RMS 0.045
      The magic-angle data hold the populations only. There is no
      50 ps in them.

 5. View: Parallel, then Global Analysis...   (optional)
      The same fit with 2 components gives tau = 12.5, 625 ps and
      RMS 0.31: two exponentials do not describe these data.
      With 3 components:
      Expected   tau = 3.01, 47.0, 799 ps;  RMS 0.053
      The rotation shows up as a time constant of its own. It belongs
      to no state: a fit of one polarization alone would count one
      state too many.
      (View: Perpendicular, 3 components: tau = 2.98, 46.7, 800 ps.)

 6. View: Magic angle again, then Global Analysis...
      The fit of step 4 is back: "Restored previous fit", 2 components,
      tau = 3.00, 800 ps. Setup and results of the analysis windows
      are kept for each view. A change of the view closes the windows;
      the results come back with the view. They are dropped when the
      data of the view change (a correction, a crop, Reset all, or a
      new G or time shift for every view but Parallel).

The same from files as measured
(demo_aniso_raw_par.csv and demo_aniso_raw_perp.csv)

The same sample, spectra and time constants, with what a measurement
adds, on the delay axis of demo_TA.csv (-2 ps to 3 ns):

Chirp         the same in both files: time zero from -0.26 ps at 450 nm
              to +0.56 ps at 750 nm
Background    not the same in the two files: an offset of 0.12 mOD
              (parallel) and 0.20 mOD (perpendicular), and pump scatter
              at 480 nm, -4.0 and -1.5 mOD at its maximum
Noise         as above, and higher at the pump scatter
truth/demo_aniso_raw_truth.csv lists time zero, the two backgrounds
and the noise for every wavelength.

 1. Load...
      demo_aniso_raw_par.csv -> Load

 2. The corrections of the first walkthrough
      Background Corr...  10 initial delay points -> Apply & Close
      Chirp Corr...       Auto-place (8 points) -> Fit & Apply
      Crop...             t range -1.5 to 2800 -> Apply Crop
      Mask lambda...      474 to 486 nm, "Set to NaN" (pump scatter)

 3. Anisotropy...
      Load perpendicular... -> demo_aniso_raw_perp.csv (the file as
      measured, not a corrected one)
      G: "Tail matching", keep the range proposed, 282 to 2730 ps
      -> Apply
      Expected   G = 1.15       (true 1.15)
                 r at 0.5 ps: 0.37 at 506 nm, -0.12 at 690 nm
                 plot on the right: "scatter: 1.6 x what the noise
                 explains"

 4. Fit trace...   (View: Anisotropy r), Components = 1
      lambda = 506:  tau = 48 ps, A_1 = 0.384    (true 50 ps, 0.38)
      lambda = 690:  tau = 50 ps, A_1 = -0.121   (true 50 ps, -0.12)

 5. View: Magic angle, then Global Analysis...
      Components = 2, IRF t_0 and IRF FWHM not fixed -> Run Fit
      Expected   tau = 3.01, 800 ps;  IRF FWHM = 0.102 ps

Things to know

- The perpendicular file gets the steps made so far by itself: its own
  average of the first 10 delays as background, the same chirp curve,
  the same crop and the same mask. The two files must have the same
  wavelengths and delays.
- The order does not matter. With the perpendicular file loaded right
  after step 1 and the corrections made afterwards (in the view
  Parallel), the data are the same to the last digit.
- The tail range proposed is the last decade of the delays: from one
  tenth of the last delay to the last delay. It is a proposal only and
  is right if the anisotropy has decayed there, as it has here. A range
  that starts too early gives a wrong factor: 5.84 to 2730 ps gives
  G = 1.48.
- The plot on the right tells the three cases apart. With the proposed
  range the G of the single wavelengths scatter 1.6 times as much as
  the noise of the tail explains. That is a little above 1 because the
  background taken from each trace is the mean of 10 delays and has
  their noise; it moves the whole trace, and with it the G of that
  wavelength. With the range from 5.84 ps the value is 20: the
  wavelengths do not agree on one G, because the anisotropy has not
  decayed. With the perpendicular file loaded before any correction it
  is 8.7: the two files have different backgrounds.
- The chirp curve found in step 2 is the true one within 0.01 ps. Both
  files are shifted with this one curve. A time zero that differs
  between the two measurements is not corrected by it; that is what
  the time shift of the Anisotropy window is for (next part).

A time zero that differs between the two measurements
(demo_aniso_raw_perp_late.csv)

The perpendicular measurement made once more, with its time zero
0.05 ps late: everything in this file happens 0.05 ps later than in
the parallel file (2.5 delay steps). Background, chirp and G are those
of demo_aniso_raw_perp.csv; the noise is its own.

 1. Steps 1 and 2 as above (demo_aniso_raw_par.csv and its corrections)

 2. Anisotropy...
      Load perpendicular... -> demo_aniso_raw_perp_late.csv
      G: "Tail matching", range 282 to 2730 ps -> Apply
      Kinetics at lambda: 506; tick "Show these delays only"
      Expected   G = 1.15 as before (the tail does not notice 0.05 ps)
                 r at 506 nm: 0.79 at 0 ps and 0.42 at 0.1 ps, where
                 0.38 is true; 0.38 at 0.5 ps, as before
                 middle plot: r rises above 1 just before time zero
                 left plot: the perpendicular curve rises later than
                 the parallel one

 3. Estimate   (range proposed: -0.5 to 0.5 ps)
      Expected   "Estimated time shift: -0.05032 ps +/- 0.00016 (from
                 the scatter between the wavelengths). It lowers the
                 sum of squares between the two data sets in the range
                 by 97 % ..."
                 (true -0.050: the perpendicular data have to go back
                 by 0.05 ps)
                 The value is in the field, and the preview shows it
                 at once: the two curves rise together, r stays near
                 0.4.
      -> Apply
      Expected   r at 506 nm: 0.40 at 0 ps, 0.38 at 0.1 ps and at
                 0.5 ps;  G = 1.15
                 main window: "... G = 1.15 (tail matching), time
                 shift -0.05032 ps"

 4. For comparison: Load perpendicular... -> demo_aniso_raw_perp.csv,
    then Estimate
      Expected   "No time shift is needed: the best value, -0.00029 ps
                 +/- 0.00013 ..., lowers the sum of squares between
                 the two data sets in the range by less than 1 % ...
                 The field was set to 0."
      -> Apply

Things to know

- The perpendicular data are moved along the delays after the
  corrections, before the two data sets are combined: perp(t - shift).
  A positive shift moves them to later delays. The values between the
  delays come from a cubic spline where the delays are equally spaced
  and from linear interpolation elsewhere, as for the time shift of
  "Subtract data..." and for the chirp correction. The parallel data
  are never touched.
- "Estimate" looks for the shift at which the perpendicular data
  follow the parallel data best in the range of delays given. Every
  wavelength gets its own factor and its own offset, because the two
  signals differ by a factor that depends on the wavelength (1 + 2 r
  against 1 - r) and each may sit on its own background. So it also
  works before the background is corrected: without the corrections of
  step 1 the value here is -0.0507 ps.
- The range has to hold the rise of the signal; without a rise there
  is nothing to align. The range proposed is 25 delays on either side
  of delay 0. The search goes as far as a quarter of the range on
  either side; a larger shift is reported as lying at the limit.
- The estimate assumes that within the range the two signals have the
  same shape at every wavelength. An anisotropy that decays while the
  signal rises breaks that: in tests on simulated data without any
  shift (IRF 0.10 ps), a rotation time of 0.5 ps gave -0.02 ps, one of
  5 ps -0.006 ps with a range of -1 to 1.25 ps. A range that ends soon
  after the rise makes this error smaller.
- The delays have to sample the rise, because the shift is made by
  interpolation. In the same tests the error was below 3 % of the
  shift with a delay step of half the IRF FWHM, about 12 % with a step
  equal to the FWHM, and more than half with a step of twice the FWHM.
- The +/- is the standard deviation from the scatter between the
  wavelengths (every wavelength pulls the shift its own way). It does
  not contain the two errors above.
- A shift that lowers the sum of squares by less than 1 % is reported
  as not needed, and the field is set to 0: shifting by interpolation
  changes the noise of the perpendicular data a little, so it is done
  only when there is something to correct.
- A new time shift changes every view but Parallel; what the analysis
  windows hold for those views is dropped, as after a new G.


Solvation demo (demo_solvation.csv)
-----------------------------------
A sixth data set: time-resolved emission of a dye in a polar solvent
(intensity, not dA). The emission band moves to the red while the
solvent relaxes around the excited molecule (dynamic Stokes shift). The
excited state itself only decays, with the fluorescence lifetime.
Not a particular molecule. No background, time zero at 0.

Band          log-normal in wavenumber, width parameter 3000 cm-1
              (FWHM 3031 cm-1), asymmetry -0.25 (tail to the red); the
              shape does not change
Peak          from 20400 cm-1 at time zero to 18600 cm-1 (relaxed):
              nu(t) = 18600 + 1800 x C(t)
              C(t) = 0.35 exp(-t/0.25) + 0.35 exp(-t/3)
                     + 0.30 exp(-t/15)        (t in ps)
              The maximum of the spectrum as the file holds it (per nm)
              moves from 487 to 533 nm.
Population    one decay, 4000 ps
Intensity     per nm, so that the integral over the wavelength is the
              number of excited molecules (10000 at time zero)
IRF           Gaussian, FWHM 0.15 ps
Noise         0.5 at all wavelengths (largest signal: 120)
Left out      the change of the radiative rate with the emission
              frequency. In a measurement the integrated emission falls
              somewhat while the band moves to the red.

 1. Load...
      demo_solvation.csv; in "Check axes orientation" set
      Signal = Intensity -> Load

 2. Look at the spectra
      Click on the map near 0.1, 1, 10 and 100 ps ("Pin" keeps a
      spectrum). The maximum moves: 493, 507, 525, 533 nm.

 3. Fit trace... at one wavelength
      Components = 3, untick "Fixed" for IRF t_0 and IRF FWHM
      lambda = 500 -> Run Fit
      Expected   tau = 4.3, 15.7, 3979 ps; all amplitudes positive
      lambda = 600 -> Run Fit
      Expected   tau = 1.6, 12.8, 3992 ps; the first two amplitudes
                 negative (the signal rises)
      The fast time constants change with the wavelength and belong to
      no state: at every wavelength the band slides past.

 4. Fit trace... on the whole band
      lambda = 420, tick "Band, to:" and enter 760, "Integral",
      Components = 1 -> Run Fit
      Expected   tau = 3987 ps, A_1 = 9988      (true 4000 ps, 10000)
                 IRF FWHM = 0.149 ps            (true 0.15)
      One exponential: the integral over the band counts the excited
      molecules, wherever the band is.

 5. Global Analysis...
      Components = 4, start values 1 / 10 / 100 / 1000, untick "Fixed"
      for IRF t_0 and IRF FWHM -> Run Fit
      Expected   tau = 0.26, 2.7, 13.1, 3998 ps;  RMS 0.51
      The four EADS are the same band with its maximum at 484, 504, 520
      and 532 nm: the sequential model describes the sliding band as a
      row of states. The first three time constants come out near the
      times of C(t) (0.25, 3, 15 ps) because C(t) is a sum of three
      exponentials here. They are not rates of steps between states.

 6. Band position...
      The window follows the band itself. With Signal = Intensity the
      box "Signal per nm: convert to a signal per cm-1" is ticked.
      taken as: "Log-normal fit";  t from 0.2;  Exponentials = 3
      -> Run
      Expected   tau_1 = 0.2452 +/- 0.0064 ps   A_1 = 653 +/- 13 cm-1
                 tau_2 = 2.918 +/- 0.077 ps    A_2 = 624.0 +/- 9.7 cm-1
                 tau_3 = 14.83 +/- 0.28 ps     A_3 = 551 +/- 12 cm-1
                 (true 0.25, 3, 15 ps;  630, 630, 540 cm-1)
                 shares 36, 34, 30 %           (true 35, 35, 30 %)
                 at long delays: 18599.60 +/- 0.57 cm-1   (true 18600)
                 at delay 0 (extrapolated): 20428 +/- 15 cm-1
                                                          (true 20400)
                 total shift: 1828.4 cm-1                 (true 1800)
                 average time: 5.553 ps                   (true 5.64)
      Plots      Display, Band at t: 0.5 -> "Band at t = 0.5 ps": the
                 data with the fitted band, position 19740.3 cm-1
                 (true 19741)
                 "Width of the band (FWHM)": 3018 to 3041 cm-1 from
                 0.2 ps on (true 3031): the band moves, its shape stays
      taken as: "First moment" -> Run
      Expected   tau = 0.255, 3.08, 14.88 ps;  total shift 1804.9 cm-1
                 at long delays: 18194.65 cm-1
      taken as: "Maximum (parabola)" -> Run
      Expected   tau = 0.2426, 2.85, 14.60 ps;  total shift 1835.9 cm-1
                 at long delays: 18564.61 cm-1

Things to know

- In step 4 the start value does not matter. With tau_init = 1, which
  the window proposes for one component, and the IRF free, the direct
  fit loses the decay: the IRF FWHM runs to 5001 ps, the range of the
  delays, and the RMS is 3601. The program takes an FWHM above half
  the range of the delays as lost and repeats the fit in two steps
  (IRF fixed, then free); the result text says so. Start values of 10,
  100 and 1000 give the same 3987 ps without a repeat.
- The band has to hold the whole emission. 480 to 620 nm with 2
  components: a rise of 0.75 ps with 15 % of the amplitude, then
  3978 ps.
- 0.2 % of the relaxed band lies beyond 760 nm. With 2 components the
  fit of step 4 shows it: 22 ps with an amplitude of 22, and 4001 ps.
- Step 6, the three ways to take the position: the times agree, the
  positions do not. The band is not symmetric, so at long delays its
  first moment lies 405 cm-1 and the top of the parabola 35 cm-1 below
  its maximum. As long as the shape of the band stays, that is a
  constant, and C(t) is the same.
- Step 6, the range of the band: the first moment needs the whole band
  inside the range at all delays. With 480 to 620 nm it gives a total
  shift of 1058 cm-1 instead of 1800, because the range cuts another
  part of the band at every delay. The log-normal fit still gives
  18600.3 cm-1 at long delays and a shift of 1820 cm-1 there.
- Step 6, the conversion: the file holds the intensity per nm. Taken
  as it is (box unticked), the maximum of the log-normal band at long
  delays is 18765 cm-1, 165 cm-1 too high: a signal per nm and a
  signal per cm-1 do not have their maxima at the same place. With the
  position given in nm and no conversion the maximum runs from 485.6
  to 532.7 nm, which is what the map of the main window shows.
- Step 6, the start of the fit: 0.2 ps lies after the instrument
  response (0.15 ps). A fit from 0 gives 0.31 ps for the first time
  constant: while the signal rises, the spectrum is a mixture of
  earlier and later times. The amplitudes are those at delay 0, found
  by extrapolation.
- Step 6, the number of exponentials: with 2 the RMS of the residual
  is 17.9 cm-1 instead of 3.2 (0.62 and 8.8 ps). With 4 the fourth
  time constant runs to the limit of the search; the result text says
  so, and C(t) is not drawn.
- "Mean" in step 4 gives the same decay with A_1 = 29.2, in the unit of
  the signal.
- truth/demo_solvation_truth_kinetics.csv lists C(t), the peak of the
  band and the maximum of the spectrum at every delay, to compare with
  what "Export..." of the Band position window saves (step 6).


Making the files again
----------------------
    make_demo_data             in MATLAB (R2020b or later) or GNU Octave
    make_demo_data('folder')   writes into that folder

The true parameters are listed at the top of the script. The noise comes
from a generator written out in the script, so the files are the same on
every run and in both programs, apart from possible differences in the
last printed digit.
