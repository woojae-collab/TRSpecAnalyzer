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
demo_preview.png                picture of demo_TA.csv and of its truth
demo_target_preview.png         the same for demo_target.csv
demo_equilibrium_preview.png    the same for demo_equilibrium.csv
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
                 IRF FWHM = 0.063 ps           (true 0.060)
                 IRF t_0  = 0.000 ps
                 RMS = 0.060 mOD
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
      RMS 0.15.
      If t_0 and FWHM are fitted in the target window instead (tested
      with start values 1 / 10 / 100 and FWHM 0.15), the answer is the
      same as in step 7. The first search stops at a poorer solution
      with an FWHM far below the delay step (0.01 ps, RMS 0.12); the
      program notices this and repeats the fit in two steps. The result
      box then says "The direct fit lost the IRF width; this is the
      result of a fit in two steps (IRF fixed, then free)."

 9. Coherence...
      Method "Use GA residual", Time window 0.15 to 4 -> Run Coherence
      Expected: two peaks, at 150 and 250 cm-1.
      150 cm-1: maxima near 556 and 606 nm, node at 580 nm
      250 cm-1: maxima near 620 and 682 nm, node at 656 nm

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
               neighbours, but kept" (noise up to 4.7 times, signal at
               least 26 times the noise)
  Chirp        time zero from -0.27 ps at 450 nm to +0.56 ps at 750 nm
               (within 0.01 ps of the true curve)
  Crop         t = -1.7 ... 2730 ps
  Artifact     coherent-artifact terms not needed (three terms would
               lower the RMS by 4 %, less than the 15 % asked for)
  Components   3:  tau = 0.604, 15.0, 700 ps;  IRF FWHM 0.063 ps;
               RMS 0.069
  Oscillations 150 and 249 cm-1, nodes near 580 and 654 nm

"Apply corrections to the data" writes the corrections into the main
window. "Send model to Global Analysis" then opens the Global Analysis
window with these values; "Run Fit" there gives the same result.

The numbers differ a little from the walkthrough because the corrections
differ: 43 instead of 10 background delays, another crop, and no mask.
The walkthrough masks 494 to 506 nm by hand. The automatic analysis
keeps these channels, because their signal is far above their noise.
With  Mask rule = "Every noisy channel"  it masks 494 ... 506 nm as
well and gives tau = 0.605, 15.0, 701 ps and RMS 0.060.

"Components used" shows "3  (automatic choice)" after the run. Tick
"Fit every number up to the largest" before Run and it lists 1 to 5;
picking another number switches the report and the plots to that fit
(nothing is fitted again). With 4, for example, the report shows an
extra time constant of 0.014 ps, a quarter of the IRF width, for an RMS
only 2 % lower: the reason the automatic choice stopped at 3.

For demo_target.csv the report proposes no correction and finds
tau = 2.00, 40.0, 799 ps. It does not find the branching: that takes the
Target Analysis steps further below.

For demo_equilibrium.csv the report proposes no correction and finds
tau = 3.71, 134.5, 1998 ps. These are the lifetimes the data show, not
the four steps of the scheme (see the equilibrium demo below).


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


Making the files again
----------------------
    make_demo_data             in MATLAB (R2020b or later) or GNU Octave
    make_demo_data('folder')   writes into that folder

The true parameters are listed at the top of the script. The noise comes
from a generator written out in the script, so the files are the same on
every run and in both programs, apart from possible differences in the
last printed digit.
