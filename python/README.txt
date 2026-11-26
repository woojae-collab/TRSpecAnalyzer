TRSpecAnalyzer - Python version 1.1.0
=====================================

The Python version of the MATLAB program TRSpecAnalyzer, for the analysis
of time-resolved spectroscopy data (TA, TR-PL, TR-Raman, TR-IR). It runs
without MATLAB. Version 1.1.0 has every window of the MATLAB program
v1.1.0 and calculates in the same way; the version number is that of the
MATLAB program whose functions it has. The user manual of the MATLAB
program describes the windows of this version as well.

   1. What is in it
   2. Install (once): Windows, macOS
   3. First steps with the demo data
   4. Updating
   5. Mouse in the plots
   6. Differences from the MATLAB program
   7. Where numbers can differ from the MATLAB program
   8. How it was checked
   9. If something goes wrong; removing; other systems; folders
  10. History


1. What is in it
----------------
  - Main window: load, 2D map, spectrum, kinetics, crosshair (click or
    drag), Scale (Linear / Log / Split), display range, colormap,
    Z range, units, pins.
  - Corrections: Background Corr, Chirp Corr, Crop, Set t=0, Mask,
    Reset all. Load & Average of several files (with Auto-align).
    Export Corrected Data. Save / Load Session.
  - Data tools (second row of buttons): Subtract data (a measurement of
    the solvent or substrate alone, with scale and time shift typed or
    estimated), Anisotropy (perpendicular data beside the loaded
    parallel data; views parallel, perpendicular, magic angle and
    anisotropy r), conversion dA <-> dT/T. What the numbers of a file
    are (dA, dT/T or Intensity) is chosen when it is loaded.
  - Analyses: Fit trace (one wavelength, or the integral or mean of a
    band), SVD, Global Analysis (with coherent-artifact
    terms, stretched components and the IRF FWHM handling of the MATLAB
    program), Target Analysis, LDA, MCR-ALS, Coherence (with
    Time-frequency, and with the subtraction of the vibrational
    coherence of the solvent measured alone), Auto Analysis.
  - The results of the analysis windows are kept when a window is closed
    and are stored in session files. They are discarded when the data of
    the main window change.


2. Install (once): Windows, macOS
---------------------------------
  Windows
  1. Install Python 3.13, 64-bit, from https://www.python.org/downloads/
     and tick "Add python.exe to PATH" on the first page of the
     installer. (3.12 works too. A version that has only just come out
     may not have the packages yet.)
  2. Double-click install_windows.bat. It makes a folder .venv inside
     this folder and downloads four packages into it (numpy, scipy,
     matplotlib, PySide6). This needs internet, takes a few minutes and
     uses about 500 MB. Nothing outside this folder is changed.
  3. Double-click run_windows.bat. The black window that opens with the
     program shows its messages; it closes when the program is closed.

  If "Windows protected your PC" appears: "More info" -> "Run anyway"
  (this can happen the first time a .bat file is run).

  macOS
  1. Install Python 3.13 for macOS from https://www.python.org/downloads/
  2. Double-click install_mac.command. A Terminal window opens; the file
     makes the folder .venv inside this folder and downloads the four
     packages into it (internet, a few minutes, about 500 MB).
  3. Double-click run_mac.command. The Terminal window that opens with
     the program shows its messages; close it after the program.

  The first time, macOS may refuse to open a .command file that came
  from the internet. Then either
    - open System Settings > Privacy & Security, scroll down and press
      "Open Anyway" for that file (on macOS 14 and earlier a right click
      on the file > "Open" does it too), or
    - use Terminal: type "cd " (with the blank), drag this folder into
      the Terminal window, press Return, then type
          bash install_mac.command
      and, to start the program,
          bash run_mac.command
  If it says that the file could not be executed because you do not have
  appropriate access privileges: in Terminal, in this folder,
          chmod +x install_mac.command run_mac.command

  The two macOS files were tried on Linux only (it has the same shell);
  they have not been run on a Mac.


3. First steps with the demo data
---------------------------------
  The steps of demo\README_demo.txt can be followed as they are. (In a
  copy of the GitHub repository the folder demo is next to the folder
  python, not inside it.) The first nine, in short:

   1. Load...             demo\demo_TA.csv -> Load
   2. Background Corr...  number of points to average 10 -> Apply & Close
   3. Chirp Corr...       Auto-place -> Fit & Apply
   4. Crop...             t range -1.5 to 2800 -> Apply Crop
   5. Mask...             From 494, To 506 -> Add region
   6. SVD...              three singular values stand far above the rest
   7. Global Analysis...  Components 3, starting values 1 / 10 / 100,
                          untick "Fixed" for IRF t_0 and IRF FWHM
                          -> Run Fit
                          expected: tau = 0.60, 15.0, 701 ps,
                          FWHM 0.063 ps, RMS 0.060
   8. Target Analysis...  Model "Sequential, all decay:
                          A -> B -> C -> GS", "tau from GA" -> Run fit
                          expected: the tau of step 7, RMS 0.060
                          (do step 7 first: t_0 and FWHM come from it)
   9. Coherence...        Method "Use GA residual", t_min 0.15, t_max 4
                          -> Run Coherence
                          expected: two peaks at 150 and 250 cm-1; the
                          150 cm-1 one is largest near 556 and 606 nm,
                          with a node at 580 nm

  Auto Analysis (steps 2 to 9 in one go): with demo\demo_TA.csv freshly
  loaded, Auto Analysis... -> Run.
     Expected on the Report tab: background from the first 43 delays; no
     mask (500 ... 502 nm "noisier than the neighbors, but kept");
     chirp from -0.27 ps (450 nm) to +0.56 ps (750 nm); crop
     t = -1.7 ... 2730 ps; artifact terms not needed (3 to 4 %);
     3 components, tau = 0.604, 15.0, 700 ps, IRF FWHM 0.063 ps,
     RMS 0.069; oscillations at 150 and 249 cm-1 (nodes near 580 and
     654 nm).
     "Apply corrections to the data" -> "Send model to Global Analysis"
     -> "Run Fit" there gives the same values. (The corrections differ
     from those made by hand in steps 2 to 5 - 43 background delays,
     another crop, no mask - so the numbers differ a little from step 7.)

  Average of several files: in Load..., select the five files of
  demo\scans together.

  Two more data sets for Target Analysis (their steps are in the later
  part of README_demo.txt):
     demo\demo_target.csv       branching model A -> B, A -> C. With the
                                two "equal:" lines: k1..k4 = 2.99, 6.02,
                                40.0, 799 ps, RMS 0.052
     demo\demo_equilibrium.csv  step in both directions, A <-> B. With
                                the line "zero: B 700 750": kf, kb, k2,
                                k3 = 4.99, 15.0, 100.0, 1998 ps,
                                RMS 0.053

  Two pairs of files for the functions added in 1.1.0 (their steps are
  in the parts "Solvent demo" and "Anisotropy demo" of README_demo.txt):
     demo\demo_solution.csv     a solute in cyclohexane, and cyclohexane
     demo\demo_solvent.csv      alone. Coherence window, "Solvent
                                coherence...": A = 0.700, time shift
                                0.0030 ps; the 801 cm-1 line of the
                                solvent goes, 215 and 590 cm-1 remain.
                                "View" of the main window then shows
                                the data without the solvent
     demo\demo_solution_raw.csv    the same two files as measured
     demo\demo_solvent_raw.csv     (background, chirp): correct the
                                solution file, then load the solvent
                                file in "Solvent coherence..."; it
                                gets the same corrections. A = 0.699,
                                time shift 0.0030 ps
     demo\demo_aniso_par.csv    one sample, two polarizations.
     demo\demo_aniso_perp.csv   "Anisotropy...", tail matching:
                                G = 1.15; r = 0.38 and -0.12 at early
                                delays, decaying with 48 ps
     demo\demo_aniso_raw_par.csv    the same two files as measured
     demo\demo_aniso_raw_perp.csv   (background, pump scatter, chirp):
                                correct the parallel file, then load
                                the perpendicular one; it gets the same
                                corrections
     demo\demo_aniso_raw_perp_late.csv   the perpendicular measurement
                                once more, its time zero 0.05 ps late.
                                "Anisotropy...", "Estimate": time shift
                                -0.050 ps

  Time-resolved emission (load it with Signal = Intensity; part
  "Solvation demo" of README_demo.txt):
     demo\demo_solvation.csv    a band that shifts to the red while the
                                solvent relaxes. Fit trace on the band
                                420 to 760 nm (Integral, 1 component,
                                start value 1000): 3987 ps, the
                                fluorescence lifetime alone.
                                "Band position...", Log-normal fit,
                                3 exponentials from 0.2 ps: 0.245,
                                2.92, 14.8 ps (made with 0.25, 3, 15)


4. Updating
-----------
  Unpack the new zip and copy the folders trspecgui and trspec and the
  README files over the old ones. The folder .venv stays, so nothing is
  installed again. (If the zip is unpacked into a new folder instead,
  run install_windows.bat there once.) Session files and settings of the
  earlier version keep working. A new version says so here if it needs
  another package.


5. Mouse in the plots
---------------------
  - Wheel: zoom in and out at the pointer. With Shift, or with the
    pointer on the numbers of the x-axis, only the x-axis. With Ctrl, or
    on the numbers of the y-axis, only the y-axis.
  - Drag with the right button: zoom into the rectangle.
  - Drag with the left button: move the zoomed plot. In plots where a
    click already means something (the 2D map with its crosshair, maps
    and spectra with a marker), press the wheel and drag instead.
  - Double click: back to the full view. In the plots where a click
    means something: right click, "Reset view". Turning the wheel back as
    far as it goes does the same.
  - A zoomed plot stays zoomed while the program draws it again: another
    delay or wavelength, a pinned curve, a correction on the same range.
    It returns to the full view when other data are loaded, with "Reset
    all", and when the unit, the scale (Linear / Log / Split) or the
    display range of that axis changes.
  - The two parts of a kinetics plot in the Split view keep one y-range;
    the two panels of a split map keep one wavelength range. The data,
    fit and residual maps of an analysis window are zoomed together.
  - Legends can be dragged and keep their place when the plot is drawn
    again.
  - Right click on a plot: copy or save the picture; hide the legend or
    put it back to its automatic place; "Reset view"; "Zoom with the
    mouse...", which lists the above.
  - On a trackpad (macOS): scrolling with two fingers is the wheel, a
    click with two fingers the right click. "Press the wheel and drag"
    needs a mouse with a wheel.


6. Differences from the MATLAB program
--------------------------------------
  - Session files are .trsession; they cannot be exchanged with the .mat
    sessions of the MATLAB program.
  - 'viridis' takes the place of the colormap 'parula'.
  - Zoom and pan with the mouse (section 5) instead of the toolbar that
    MATLAB shows on a plot; no data tips. By its code the MATLAB program
    sets the wavelength range of the main spectrum back whenever it
    draws it again; here the zoom stays.
  - Split view of the kinetics: the right part carries no numbers on its
    ordinate (both parts have one y-range), and the legend stands in the
    wide right part.
  - Dragging the crosshair from one panel of the split map into the other
    follows the pointer.
  - The windows are always light (the dark mode of the system is not
    followed), and the plots use the same font on every computer (DejaVu
    Sans, which comes with matplotlib).
  - On a small screen or with a large display scaling: titles and axis
    labels get smaller letters or a shorter wording when the room is
    short, rows of buttons continue on a second line, and panels scroll
    in a low window. The status text at the top right of the main window
    ends in "..." when it is too long; point at it to read all of it.
  - Target Analysis: "Export:" stands in front of its buttons instead of
    above them, and in a low window the scheme box and the parameter
    table get lower, so that the result box keeps its room.
  - Auto Analysis: long lines of the report continue on the next line;
    in a narrow window the legends of the two small plots and the label
    above the kinetics use shorter words; logarithmic axes carry plain
    numbers (0.03) instead of powers of ten.
  - A "Please wait" box is shown while a fit runs (Global Analysis: a few
    seconds to some tens of seconds).
  - Preferences (last folder, units, colormap, ...) are kept in
    %LOCALAPPDATA%\TRSpecAnalyzer\prefs_python.json  (Windows) or
    ~/Library/Preferences/TRSpecAnalyzer/prefs_python.json  (macOS).
  - Fitting and calculations are the same as in the MATLAB program.


7. Where numbers can differ from the MATLAB program
---------------------------------------------------
  The method is the same. What the data do not determine is left to the
  path the search takes, and that path is not the same in the two
  programs.

  - Auto Analysis
      * "three terms lower the RMS by ... %" (artifact terms that are
        only tried): 3 to 4 % for demo_TA, depending on the path. The
        decision "not needed" does not change.
      * "Coherent artifact: On" for data without an artifact: the fit has
        two neighboring solutions, and small changes decide which one
        the search finds. demo_TA with the default mask rule: IRF FWHM
        0.084 ps in both. With Mask rule "Every noisy channel":
        0.080 ps (RMS 0.0575) here, 0.084 ps (RMS 0.0569) with the
        MATLAB code. The time constants differ by 0.6 % or less.
      * The fit with one component too many (rejected in both) can be
        rejected for another reason.
      * Data without a rise (they start after time zero) for which a
        time zero is assumed all the same: in one synthetic test set the
        search drove the IRF width to zero here and not with the MATLAB
        code; the width is then fixed (the report says so) and the fit
        is worse. If the report shows "IRF FWHM ... (fixed)" for such
        data, crop or check the first delays.
  - Target Analysis: rate constants marked "[not determined]" can come
    out with other numbers (README_demo.txt: "71 % to B" in step 3 of
    the target demo; "4.0, 50 and 124 ps" in step 3 of the equilibrium
    demo). The lifetimes, the RMS and everything else the data fix are
    the same. For the same reason the note "The direct fit lost the IRF
    width ..." of README_demo.txt, step 8, does not appear here: the
    search reaches the same answer directly.
  - Coherence, "Per-pixel exp fit": as with every fit with free time
    constants, the search can end at another local fit at single
    wavelengths.
  - MCR-ALS, "Random" start: other random numbers than in MATLAB, so the
    same seed gives another start. Within this version a seed repeats a
    run exactly.


8. How it was checked
---------------------
  - Every window was compared with the code of the MATLAB program, run
    in GNU Octave: texts, tables and messages are identical; fitted
    values agree to 7 digits or more (Global and Target Analysis, where
    both searches end at the same fit), LDA and MCR-ALS results to 10
    digits or more, the maps of Coherence and Time-frequency to 12
    digits or more, the table of fits of Auto Analysis to 6 digits or
    more from the same corrected data.
  - The walkthroughs of demo\README_demo.txt give the numbers written
    there (section 3).
  - It was used on Windows by the author while it was built.
  - Not checked: macOS and Linux desktops (the program was only run
    on Linux without a screen, for the tests), and with them the two
    files for macOS.


9. If something goes wrong; removing; other systems; folders
------------------------------------------------------------
  - If installing or starting fails, copy the text of the black window.
  - If an error box appears inside the program, copy what "Show
    Details" shows.
  - Removing: delete this folder. (The settings file is at the place
    named in section 6.)
  - Linux, or macOS by hand:
        python3 -m venv .venv
        .venv/bin/python -m pip install -r requirements.txt
        .venv/bin/python -m trspecgui
  - Folders:  trspecgui\  the windows;  trspec\  the calculations;
              demo\  demo data and their description.


10. History
-----------
  1.1.0  (2026-11-26)
     The new functions of the MATLAB program 1.1.0, calculated in
     the same way (CHANGELOG.md of the repository has the details):
     - Signal kind: what the numbers of a file are (dA, dT/T or
       Intensity) is chosen when it is loaded; axes, color bar and
       export names follow. Conversion dA <-> dT/T by the exact
       relation (button "dA <-> dT/T...").
     - "Subtract data...": a reference measured on the same axes
       (solvent, substrate) is subtracted, data - scale x
       reference(t - shift); scale and shift typed or estimated by
       least squares. The reference is loaded as measured: it follows
       crop and the corrections applied from Auto Analysis.
     - "Anisotropy...": perpendicular data are loaded beside the
       (parallel) data and corrected alike; "View" chooses parallel,
       perpendicular, magic angle or anisotropy r, and the view shown
       is what the analysis windows work on. G typed or from tail
       matching; the range proposed for it is the last decade of the
       delays. Tail matching finds one G for all wavelengths; a third
       plot of the window shows the G of every single wavelength in
       the tail range with its uncertainty, as a check (display only).
       A time zero that differs between the two measurements: "Time
       shift of the perpendicular data", typed or found by "Estimate"
       (the shift for which the perpendicular data follow the parallel
       data best around the rise of the signal, at every wavelength
       with its own factor and offset; with a standard deviation from
       the scatter between the wavelengths).
     - Setup and results of the analysis windows are kept for each
       view: a change of the view closes the windows, and the results
       of the view chosen come back with it (also from a session
       file). Whatever changes the data of a view drops what is kept
       for it: a correction, crop or Reset all drops everything, a new
       G or a new time shift of the perpendicular data everything but
       the view Parallel, a new threshold the view r.
     - Coherence window, "Solvent coherence...": the vibrational
       coherence of the solvent, measured alone, is taken out of the
       residual before the FFT; factor and time shift are estimated on
       the real and imaginary parts of the Fourier transforms in a band
       of a solvent mode (after M. Liebel, thesis, Oxford 2014).
       While it is subtracted there ("Apply"), "View" of the main window
       offers the data without the solvent: "Solvent coherence
       subtracted" (its oscillating part, in the time window of the
       analysis) or "Solvent signal subtracted" (the whole solvent
       signal at all delays, with the factor and the shift found
       here). The 2D map, the export and every analysis use the view
       shown. The analysis windows close when the view changes; their
       setup and results are kept for each view and come back with it.
     - Coherence window, before and after: with a solvent subtraction
       the window keeps the result from before it as well. A bar above
       the plots chooses which map is shown, the curve plots carry the
       curve from before in gray, and "Side by side..." opens a window
       with the two maps on one color scale. The exports write what
       the plots show.
     - Coherence window, "Color range": the color range (z scale) of
       the |FFT|^2 map can be typed, for instance a small maximum to
       bring out weak modes; "Auto" returns to the automatic range.
     - Fit trace: a band of wavelengths can be fitted instead of one
       wavelength ("Band, to:"), as its integral over the spectral axis
       or as its mean. In the view Anisotropy r, r is calculated from
       the band values of the two polarizations.
     - "Band position..." (new window): where a band lies at every
       delay, for a spectral shift such as a dynamic Stokes shift. The
       position is taken as the first moment, as the maximum of a
       parabola through the top, or as the maximum of a fitted
       log-normal band (chosen in the window), in cm-1 or nm; an
       intensity per nm is converted to one per cm-1. The position
       along the delays is fitted with 1 to 4 exponentials and a
       constant; the window shows the band with its model, the
       position with the fit, C(t) and the width, and exports them.
     - Interpolation along the delays: the chirp correction and the
       time shifts ("Subtract data...", perpendicular data) take the
       values between the delays from a cubic spline inside every
       stretch of equally spaced delays and from linear interpolation
       elsewhere. Up to 1.0.2 the chirp correction was linear
       everywhere, which damped fast oscillations (8 delays in a
       period: 9 % of the power on average, 0 to 14 % depending on
       the wavelength). Chirp-corrected data therefore differ
       slightly from those of 1.0.2. For scripts: trspec.chirp has
       delay_runs, spline_curv, interp_runs and chirp_shift.
     - Details of these functions: tail
       matching finds G by orthogonal regression (noise in both data
       sets is allowed for); in the view Anisotropy r the automatic
       color range of the map is at most -0.5 ... 0.5; a subtraction
       applied, changed or removed in the main window switches the
       solvent subtraction of the Coherence window off.
     - Demo data for them: demo_solution.csv with demo_solvent.csv
       (and the same pair as measured, demo_solution_raw.csv with
       demo_solvent_raw.csv), demo_aniso_par.csv with
       demo_aniso_perp.csv (and the same pair as measured,
       demo_aniso_raw_*.csv), and demo_solvation.csv
       (time-resolved emission with a spectral shift, for the fit of a
       band); the steps and the expected numbers are in README_demo.txt.
     - Fits with a free IRF width (Fit trace, Global Analysis, Target
       Analysis): an FWHM above half the range of the delays counts as
       lost, like one the delays do not resolve, and the fit is
       repeated in two steps (IRF fixed, then free). Before, a decay
       much slower than its start value (a ns decay, start 1 ps) could
       end with the FWHM at the range of the delays and a meaningless
       time constant. A width that stays there is marked as not
       determined.
     - Auto Analysis: the status line after a run names the steps
       that could not be decided ("Not determined: chirp") instead of
       saying that no correction was needed.
     - Plots: each of the three maps of the Global Analysis window
       (Experimental, Global fit, Residual) has a color bar with the
       name of the signal; the Global fit map had none. On
       a zoomed logarithmic y-axis, numbers of more than five
       characters are written with a common factor (x 10^n), so that
       the label of the axis stays inside the picture. The LOF axis
       of the convergence plot of the MCR-ALS window shows plain
       numbers (1.37, 1.38) when the LOF changes by less than a
       decade; before, these numbers were cut off at the left edge.
     - Plots in rows and columns: in every window with several
       plots, the plots of a row have exactly the same upper and
       lower edge and the plots of a column the same left edge. A
       plot without a color bar does not leave the room of the bar
       of a map free: in the Coherence, Target and Lifetime Density
       Analysis windows, in the rows under the maps of the Global
       Analysis window and on the Residual tab of the Auto Analysis
       window, the room right of such a plot is about as wide as the
       blank left of the label of its y-axis. These stay exactly as
       wide as the map above or below them: the spectrum under the
       map of the Lifetime Density Analysis window (it has the
       wavelength axis of the map) and of the Time-Frequency window,
       the plots of the SVD window and of the Corrections tab of the
       Auto Analysis window. The DADS and the EADS of its Components
       tab are equally wide.
     The residual of the per-pixel fit of the Coherence window is kept
     until its settings change. Sessions of earlier versions load as
     before. For scripts, trspec.coherence has the new functions
     coh_spectrum, solvent_vc, solvent_shifted, solvent_band and
     fit_solvent_coherence.

  1.0.2  (2026-11-22)
     First version with a number: the program of preview 12, under the
     number of the MATLAB program whose functions it has. Its texts
     (buttons, tooltips, messages, reports) use American spelling, as
     those of the MATLAB program 1.0.2 do. With install_mac.command and
     run_mac.command for macOS.

  Preview builds (November 2026)
     12  Main window before data are loaded, Split view: no numbers and
         no "Wavelength (nm)" between the two panels of the map.
     11  Auto Analysis keeps its last run: shown again when the window
         is reopened ("[restored]"), stored in session files.
     10  Zoom and pan with the mouse in every plot.
      9  The status text of the main window stays on the first row.
      8  Auto Analysis. Fixed: the evenly spaced delay grid of
         Coherence, Time-frequency and Auto-align could end a rounding
         error after the last delay.
      7  Coherence and Time-frequency.
      6  LDA and MCR-ALS.
      5  Target Analysis.
      4  Fixes for legends that had been dragged.
      3  Legends can be dragged; checkbox frames in the dark mode of
         Windows; texts fitted to small screens.
   1, 2  Main window, corrections, Load & Average, sessions, Fit trace,
         SVD, Global Analysis.


License: MIT (see LICENSE).  Author: Woojae Kim, ULS Lab, Department of
Chemistry, Yonsei University.  The code was written with Claude Code.
