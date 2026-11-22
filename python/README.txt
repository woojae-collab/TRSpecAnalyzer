TRSpecAnalyzer - Python version 1.0.2
=====================================

The Python version of the MATLAB program TRSpecAnalyzer, for the analysis
of time-resolved spectroscopy data (TA, TR-PL, TR-Raman, TR-IR). It runs
without MATLAB. Version 1.0.2 has every window of the MATLAB program
v1.0.2 and calculates in the same way; the version number is that of the
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
  - Analyses: Fit trace, SVD, Global Analysis (with coherent-artifact
    terms, stretched components and the IRF FWHM handling of the MATLAB
    program), Target Analysis, LDA, MCR-ALS, Coherence (with
    Time-frequency), Auto Analysis.
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
  - Coherence window: t_min / t_max, the averaging range and the display
    range each take one row instead of two.
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
