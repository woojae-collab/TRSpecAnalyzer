# Changelog

Version numbers are major.minor.patch:

- patch (1.0.x): bug fixes and changes of wording
- minor (1.x.0): new features
- major (x.0.0): changes that break older session files or the way the
  program is used

The version of a copy of the program is shown in the title of its main
window and returned by `TRSpecAnalyzer.Version`.

## 1.1.0 (2026-11-26)

Six new functions: conversion ΔA ↔ ΔT/T, subtraction of a solvent or
substrate measurement, anisotropy from two polarizations, removal of the
solvent coherence in the Coherence window (with views of the main
window that show the data without the solvent), fits of a
band of wavelengths in the Fit trace window, and a window that follows
the position of a band along the delays. The existing analyses calculate as in 1.0.2,
with two exceptions (see below): fits with a free IRF width that ran to
the range of the delays, and the chirp correction, which now
interpolates with a cubic spline where the delays are equally spaced.
Chirp-corrected data, and the results obtained from them, therefore
differ slightly from those of 1.0.2.

- Signal kind: what the numbers of a file are (ΔA, ΔT/T or Intensity) is
  chosen under `Signal` when data are loaded, in `Check axes orientation`
  and in `Load & Average Datasets`. The names on the axes, the color bars
  and the plots of all windows, the heading `Max |…|` of the Load &
  Average table and the file name suggested by `Export Corrected Data...`
  follow it; the label "ΔA or Intensity" is gone. The choice is
  remembered for the next file.
- Main window: a new second row of buttons, the data tools:
  `Subtract data...`, `Anisotropy...`, `View` and `ΔA <-> ΔT/T...`, with
  a text that says what the data shown are. The analysis buttons are now
  the third row, the display range, export, session and unit controls the
  fourth.
- `ΔA <-> ΔT/T...`: transient-absorption data are converted by the exact
  relation ΔT/T = 10^(−ΔA) − 1, ΔA = −log10(1 + ΔT/T), with ΔA in OD or
  mOD. Like a crop, the conversion takes the data with all corrections
  made so far as the new starting point. The program asks before it
  converts numbers that do not fit the unit chosen; values with
  ΔT/T ≤ −1 are left empty. Intensity data are not converted.
- `Subtract data...`: a reference data set measured at the same
  wavelengths and delays (solvent or substrate alone) is subtracted,
  data − scale × reference(t − shift), before the background and chirp
  corrections. Scale and time shift are typed or estimated by least
  squares in a range of wavelengths and delays; a preview shows data,
  reference and difference. The values between the delays come from a
  cubic spline inside every stretch of equally spaced delays, so that a
  feature only a few delay steps wide keeps its height, and from linear
  interpolation elsewhere (see "Interpolation along the delays" below).
  A background taken from the first N
  delays is recalculated after the subtraction.
- Second data sets follow the corrections: the program records the steps
  that make corrected data the starting point (crop, conversion,
  `Apply corrections to the data` of Auto Analysis). A reference,
  perpendicular or solvent file is loaded as measured at any time and is
  brought to the state of the data; it follows every later correction.
- `Anisotropy...`: the loaded data are taken as the measurement with
  parallel polarization, and the data measured with perpendicular
  polarization are loaded beside them. Both receive every correction
  alike. `View` chooses what the main window shows and every analysis
  window works on: parallel, perpendicular (× G), magic angle
  (par + 2 G perp) / 3, or the anisotropy r = (par − G perp) /
  (par + 2 G perp). G is typed or taken from tail matching: the factor
  between the two data sets in a range of delays, found by orthogonal
  regression through the origin, which allows for the noise of both.
  The range proposed is the last decade of the delays (from one tenth of
  the last delay to the last one); it is right only if the anisotropy has
  decayed there.
  Tail matching finds one G for all wavelengths. As a check, a third
  plot of the window shows the G of every single wavelength in the tail
  range with its uncertainty (from the scatter of the points of that
  range about the line par = G × perp), and the G in use as a line. Its
  title gives the scatter of these values about the G in use in units of
  that uncertainty: about 1 when one G holds for the whole spectrum,
  larger when the anisotropy has not decayed in the range or the two
  files differ in more than a factor (a background, for example).
  Wavelengths whose G is determined to less than 10 % are drawn in grey
  and not counted. The plot is display only: the calculation uses the
  one G.
  A time zero that differs between the two measurements is corrected
  with `Time shift of the perpendicular data`: the perpendicular data
  are moved along the delays, perp(t − shift), interpolated as in
  `Subtract data...`, after the corrections and before the two data sets
  are combined. The shift is typed, or `Estimate` finds it: the shift
  for which the perpendicular data follow the parallel data best in a
  range of delays around the rise of the signal (proposed: 25 delays on
  either side of delay 0), at every wavelength with its own factor and
  offset. The value comes with a standard deviation from the scatter
  between the wavelengths and with the change of the sum of squares; a
  shift that lowers the sum of squares by less than 1 % is reported as
  not needed. The estimate assumes that the two signals have the same
  shape in the range, so it is biased when the anisotropy decays while
  the signal rises, and it needs delays that sample the rise (tests on
  simulated data: error below 3 % of the shift with a delay step of
  half the IRF FWHM, about 12 % with a step equal to it). The preview
  follows the shift typed; `Show these delays only` limits its two
  time plots to the range.
  r is left empty where |par + 2 G perp| is below 3 times its noise (the
  factor can be changed). In the view r the automatic color range of the
  2D map is at most −0.5 … 0.5. While perpendicular data are loaded,
  `Subtract data...` is not available, and in the view r the corrections
  of Auto Analysis cannot be applied.
- Results per view: setup and results of the analysis windows (Global
  Analysis, Fit trace, Target Analysis, SVD, LDA, MCR-ALS, Coherence,
  Time-Frequency, Band position, Auto Analysis) are kept for each of the
  four views. A change of the view closes the analysis windows, stores
  what they hold for the view that is left and brings back what is kept
  for the view chosen; a view that was not analyzed yet starts with the
  setup of the view before and without results. Results that existed before the
  perpendicular data were loaded are those of the view Parallel, and
  stay with the data when the perpendicular data are removed. What is
  kept for a view is dropped when its data change: a correction, a crop,
  `Reset all`, a unit switch or new data drop everything, as before; a
  new G or a new time shift of the perpendicular data drops everything
  but the view Parallel, a new threshold only the view r. Sessions store
  the results of all views.
- Interpolation along the delays (chirp correction; time shift of
  `Subtract data...` and of the perpendicular data): inside every stretch
  of equally spaced delays (at least 4 delays whose steps agree within
  10 %) the values between the delays come from a cubic spline through
  the delays of that stretch; everywhere else they are interpolated
  linearly, as the chirp correction did everywhere up to 1.0.2.
  Linear interpolation flattens what changes within a few delays: an
  oscillation with 8 delays in a period (801 cm⁻¹ at 5 fs steps) lost
  9 % of its power on average and between 0 and 14 % depending on the
  wavelength, because the loss follows the part of a delay step by which
  time zero is shifted; that put a ripple on the spectral profile of the
  mode (30 periods across the spectrum of `demo_solution_raw.csv`). With
  the spline the same line keeps 99.8 % at every wavelength, a 1444 cm⁻¹
  line (4.6 delays in a period) 98 % instead of 73 %, and the coherent
  artifact at time zero is reproduced to 0.002 mOD instead of 0.08 mOD
  (exact model of the solvent demo). A spline is not used across a
  change of the delay step, because it overshoots there: one spline over
  all delays put errors of up to 2.8 mOD on the same demo data where the
  step goes from 5 fs to 0.89 ps. Between unevenly spaced delays a fast
  oscillation is therefore still damped; keep the delays equally spaced
  as far as oscillations are to be followed. The noise of
  chirp-corrected data is a little higher than before (linear
  interpolation averaged it down to between 50 and 100 % of its variance,
  depending on the wavelength; the spline to between 75 and 100 %), and
  the RMS of a fit to such data rises accordingly.
  Results of the interpolations are kept between redraws, so that the
  main window stays as fast as before.
- Coherence window, `Solvent coherence...`: the vibrational coherence of
  the solvent, measured alone, is taken out of the residual before the
  FFT: residual − A × solvent coherence(t − shift). The solvent file is
  the reference data set of `Subtract data...`; the corrections of the
  data are applied to it. A and the time shift are typed or estimated by
  least squares on the real and imaginary parts of the Fourier transforms
  of the two residuals, over a range of wavelengths, in a frequency band
  that holds a solvent mode (after M. Liebel, thesis, University of
  Oxford, 2014, Section 4.2); `Strongest solvent band` proposes the band.
  The time shift is made by band-limited interpolation. With `Apply`
  the data of the main window are not changed. Not available while
  perpendicular data are loaded. A subtraction applied, changed or
  removed with `Subtract data...` switches it off, because its factor
  was found on the data as they were.
- Solvent views of the main window: while the solvent coherence is
  subtracted in the Coherence window (`Apply`), `View` of the main
  window offers `Data`, `Solvent coherence subtracted` and `Solvent
  signal subtracted`, so that the 2D map, the kinetics, the export and
  every analysis window can work on data without the solvent.
  `Solvent coherence subtracted` takes out what the Coherence analysis
  subtracts from its residual, A × solvent coherence(t − shift) in the
  time window of the analysis (the solvent data of that window minus
  their smooth part); the delays outside the window and the coherent
  artifact stay as they are. `Solvent signal subtracted` takes out the
  solvent data themselves at all delays, with the corrections of the
  data and the same factor and shift; the coherent artifact goes as
  well, and the noise of the solvent file is added at all delays. The
  analysis windows close when the view changes; their setup and results
  are kept for each view and come back with it, as for the views of two
  polarizations. In a solvent view the Coherence window subtracts
  nothing more and compares with the residual before the view's
  subtraction; `Solvent coherence...` opens in the view `Data` only,
  where factor, shift and time window are set. The views follow that
  analysis and end when it subtracts nothing any more (`Remove
  subtraction`, `Reset` there, another solvent file, a change in
  `Subtract data...`). While a solvent view is shown, another solvent
  file cannot be loaded, `Subtract data...` cannot be changed and
  perpendicular data cannot be loaded. A crop, a conversion or the
  corrections of Auto Analysis make the view shown part of the data,
  like the other corrections, and `Reset all` returns to the view
  `Data`. A color range typed in for the 2D map stays when the view
  changes. The title of the 2D map, the suggested name of the exported
  file (`_solvcoh`, `_solvsig`) and the text beside `Signal:` name the
  view; that text now also names a subtraction made with `Subtract
  data...` (`reference subtracted (scale …)`). Stored in sessions, with
  what is kept for the views not shown.
- Coherence window, before and after the solvent subtraction: when a
  solvent coherence is subtracted (in the analysis, or by a solvent view
  of the main window), a
  run keeps the result from before the subtraction as well: the residual
  with the solvent coherence still in it, transformed with the same
  settings and given in the units of the map after the subtraction
  (divided by the same normalization). A bar above the plots chooses
  what the map shows (`Map shows: After subtraction / Before
  subtraction`). The four curve plots (power spectrum at λ, average over
  the λ range, spectral profile with its phase, residual trace) carry
  the curve from before in gray behind the result (`"Before" in gray`).
  `Side by side...` opens a window with the two maps on one color scale
  and the two power spectra laid over each other. Each spectrum is
  exactly as wide as the map above it, so a line of the map stands
  straight above its peak. The window follows the cursor, the display
  range and the scale of the Coherence window, and a click on one of its
  maps moves the cursor. The exports write what the plots show: the map
  on display, and in the three curve exports further columns with the
  values from before while the gray curves are on.
- Coherence window, `Color range`: the color range (z scale) of the
  |FFT|² map can be set by hand in the same bar, for instance a small
  maximum to bring out weak modes. The two fields show the range in
  use, in the units of the map also with the log scale. A typed value
  fixes the range, which then stays when the display range or the map
  shown changes; `Auto` returns to the automatic range. The side-by-side
  window has its own fields. Range, map shown and the state of the gray
  curves are kept with the result.
- `Fit trace...`: a band of wavelengths can be fitted instead of one
  wavelength. With `Band, to:` ticked, the trace is the integral of the
  signal over the spectral axis of the band (`Integral`, the default) or
  its mean (`Mean`). Wavelengths without data (masked) are left out, and a
  delay at which one of the others has no value stays empty, so that the
  same wavelengths contribute at every delay. In the view Anisotropy r, r
  is calculated from the band values of the two polarizations. Plot
  titles, result text, suggested file name and parameter file name the
  band; its settings are stored in sessions.
- `Band position...` (new window, third row of buttons): where a band
  lies at every delay, for a spectral shift such as the dynamic Stokes
  shift of an emission band. The band is a range of the spectral axis,
  taken as positive or negative. Its position is given in cm⁻¹ or nm and
  taken in one of three ways, chosen in the window: `First moment` (the
  mean of the axis weighted with the signal; no assumption about the
  shape, but the band has to lie inside the range at all delays),
  `Maximum (parabola)` (the top of a parabola through the points above a
  share of the maximum, 0.8 by default), or `Log-normal fit` (the
  maximum of a log-normal band fitted to the whole range: height,
  position, width, asymmetry). For a signal that is a density along the
  spectral axis, as an emission intensity per nm, the signal is
  converted to one per unit of the position axis (× λ² / 10⁷); this is
  the default for intensity data on a wavelength axis and is left off
  for ΔA. A delay gets a position only when its band reaches a share of
  the largest signal (5 % by default).
  The position along the delays is fitted with a sum of 1 to 4
  exponentials and a constant, position(t) = ν∞ + Σ Aᵢ exp(−t/τᵢ), from
  a start delay that should lie after the instrument response. The
  result gives τᵢ and Aᵢ with standard deviations, the position at long
  delays and the one extrapolated to delay 0, the total shift and the
  average time Σ Aᵢτᵢ / Σ Aᵢ. Four plots show the band at one delay
  with its model, the position with the fit,
  C(t) = (position − ν∞) / total shift, and the width of the band.
  Starting values are found by the program (three starts, the best fit
  is kept) or typed. A time constant that runs to the limit of the
  search is marked, and C(t) is then not drawn. `Export...` saves one
  row per delay: position, height, width, area, asymmetry, fit, C(t).
  The result is kept when the window is closed, per anisotropy view and
  in sessions, and dropped when the data change, as for the other
  analysis windows.
  With `demo_solvation.csv` (band made with its maximum moving from
  20400 to 18600 cm⁻¹ with 0.25, 3 and 15 ps, shares 35, 35, 30 %), the
  log-normal fit and 3 exponentials from 0.2 ps give 0.245 ± 0.006,
  2.92 ± 0.08 and 14.8 ± 0.3 ps, shares 36, 34, 30 %, 18599.6 ± 0.6 cm⁻¹
  at long delays and a total shift of 1828 cm⁻¹. The first moment gives
  0.255, 3.08, 14.9 ps and the parabola 0.243, 2.85, 14.6 ps; their
  positions lie 405 and 35 cm⁻¹ below the maximum, because the band is
  not symmetric.
- Fits with a free IRF width (`Fit trace...`, Global Analysis, Target
  Analysis): a fitted IRF FWHM above half the range of the delays now
  counts as lost, as a width the delays do not resolve already did, and
  the fit is repeated in two steps (time constants with the IRF held at
  its start values, then everything free). Before, a decay much slower
  than its start value could send the search the other way: with one
  component, the start value 1 ps and the IRF free, the 4000 ps decay of
  `demo_solvation.csv` ended with τ = 1 ps and an FWHM of 5001 ps, the
  range of the delays (RMS 3601 instead of 15). The result text says
  when the result comes from the repeat. A width that stays above half
  the range, as for data without a step, is marked as not determined
  (`[not determined by the data]`, in Target Analysis `[not determined]`)
  with the advice to fix it. Fits that worked before give the same
  numbers.
- Coherence window: the two limits of the time window, of the averaging
  range and of the display range share one row each, which makes room for
  the new button and for the result text on small screens. The residual
  of `Per-pixel exp fit` is kept until its settings change, so a further
  run with another time window is fast.
- Color bars of the split map of the main window and of the three maps
  of the Global Analysis window (MATLAB version): the label of the bar
  was cut off at the edge of its panel, in the main window when the tick
  labels were three or four characters wide (data in the hundreds, as in
  `demo_solvation.csv`), in the Global Analysis window always. The bar
  now sits at the left edge of its own space, and the tick labels and the
  label have the rest of the width; in the Global Analysis window that
  space is 75 pixels wide instead of 50.
- Global Analysis window: the `Global fit` map shows a color bar as the
  two other maps do. It had none, because it is drawn with the color
  range of the `Experimental` map; the `Residual` map has its own range.
- Main window, `Split` view (MATLAB version): the upper and the lower
  panel of the 2D map are exactly equally wide. Before, the lower panel
  was a few pixels narrower, because the last number of its spectral
  axis reaches beyond the plot. The color bar reaches from the lower
  edge of the lower panel to the upper edge of the upper one. The trace
  plot and the residual plot of the `Fit trace...` window are equally
  wide in the same way.
- Plots in rows and columns: in every window with several plots, the
  plots of a row have exactly the same upper and lower edge and the
  plots of a column the same left edge. A plot without a color bar does
  not leave the room of the bar of a map free: in the Coherence, Target
  and Lifetime Density Analysis windows, in the rows under the maps of
  the Global Analysis window and on the Residual tab of the Auto
  Analysis window, the room right of such a plot is about as wide as
  the blank left of the label of its y-axis. These stay exactly as wide
  as the map above or below them: the spectrum under the map of the
  Lifetime Density Analysis window (it has the wavelength axis of the
  map) and of the Time-Frequency window, the plots of the SVD window
  and of the Corrections tab of the Auto Analysis window. The DADS and
  the EADS of its Components tab are equally wide; the kinetics above
  the EADS are narrower by the room of their legend, which stands
  beside the plot.
  In the MATLAB version the plot boxes are now
  placed at fixed distances from the edges of their cells, instead of
  taking the room that the labels of each plot happen to leave; the
  plots are slightly smaller for it. The Python version had the same
  differences wherever a color bar or a legend stood beside a plot.
  Not in line by design: the 2D map of the main window against the
  spectrum and the kinetics beside it (each panel has its own row of
  controls), and the singular values plot in the setup panel of the SVD
  window.
- Auto Analysis window: the status line after a run names the steps
  that could not be decided. `Done. No correction was needed.` was also
  shown when the chirp, the background or the mask had not been
  determined at all (for example emission data whose red side rises
  slowly: no time zero along the whole spectrum). The line then reads
  `Done. No correction proposed. Not determined: chirp (see Report).`;
  the reason is in part 5 of the report, as before. The line may take
  two lines of the window now.
- Sessions store the new items: signal kind, reference and perpendicular
  data, the recorded steps, and the settings of the three tools. Session
  files of earlier versions load as before; their data are taken as ΔA.
  In a session saved by an earlier version after a crop, the steps made
  before the crop are not recorded, and a second data set can be loaded
  only after `Reset all`.
- Demo data: new synthetic files with known answers. `demo_solution.csv`
  with `demo_solvent.csv` (a solute in cyclohexane and cyclohexane alone,
  for `Solvent coherence...` and `Subtract data...`;
  `demo_solution_raw.csv` and `demo_solvent_raw.csv` are the same pair
  as measured, with background and chirp),
  `demo_aniso_par.csv` with `demo_aniso_perp.csv` (one sample, two
  polarizations; `demo_aniso_raw_par.csv` and `demo_aniso_raw_perp.csv`
  are the same pair as measured, with background, pump scatter and
  chirp; `demo_aniso_raw_perp_late.csv` is the perpendicular measurement
  once more with its time zero 0.05 ps late, for the time shift), and
  `demo_solvation.csv` (time-resolved emission of a
  dye whose band shifts to the red while the solvent relaxes, for the fit
  of a band). `README_demo.txt` gives the steps and the numbers to
  expect; `make_demo_data.m` writes the files.
- Manuals: a new Section 13, `Band position`; Auto Analysis, saving
  and "Things to know" are now Sections 14, 15 and 16. Sections 1 to 6,
  10, 12, 15, 16 and the glossary describe the other new functions.
  Section 2 has the steps for the new demo files. Every window now has a
  picture in which its plots are marked (a), (b), ...: 19 screenshots of
  the Python version, taken with the demo data, the same in both
  manuals. The windows of the MATLAB version are arranged in the same
  way.
- References: the new file `REFERENCES.md` names, part by part, the
  publications the methods of the program rest on. `README.md` links
  to it.
- Python version: the same functions, calculated in the same way.
  For scripts, `trspec.coherence` has the new functions `coh_spectrum`,
  `solvent_vc`, `solvent_shifted`, `solvent_band` and
  `fit_solvent_coherence`, the result of `compute_coherence` has the new
  field `norm_scale` (what the map was divided by), the new module
  `trspec.bandpos` has
  `band_position`, `fit_band_shift`, `lognorm_shape` and `lm_fit` (the
  Levenberg-Marquardt search these two use), and `trspec.chirp` has
  `delay_runs`, `spline_curv`, `interp_runs` and `chirp_shift` (the
  interpolation along the delays; `Dataset.correct_chirp` uses it).
- Python version, plots: the color bars of the three maps of the Global
  Analysis window carry the name of the signal, as in the MATLAB
  version. On a zoomed logarithmic y-axis,
  numbers of more than five characters (300000, 0.0002) are written with
  a common factor (× 10^n); before, they pushed the label of the axis off
  the picture (singular values of data with large numbers). The text
  `N=...` of the singular-value plot no longer lies on the numbers of the
  axis when the view is zoomed away from its line. In the MCR-ALS
  window, the LOF axis of the convergence plot shows plain numbers
  (1.37, 1.38) when the LOF changes by less than a decade, which is the
  usual case; before, these numbers were cut off at the left edge and
  the label of the axis was not on the picture.
- Python version, Auto Analysis window: the status line has room for
  three lines. On a small screen with a scaled display (1280 x 672
  points at 150 %) the longest text, which names three steps that could
  not be determined, did not fit into two.

## 1.0.2 (2026-11-22)

Auto Analysis keeps its last run, as the other analysis windows do; two
small fixes; American spelling in the texts. No calculation is changed:
results are the same as with 1.0.1.

This version was not published as a release of its own; everything in
it is part of 1.1.0.

- Auto Analysis: the last run is kept when the window is closed, as in the
  other analysis windows. A window opened again shows the report, the
  plots, the options of that run, `Components used` and `Kinetics at` as
  they were, with `[restored]` in the status line. The run is stored in
  session files. It is discarded by a new run, and when the data of the
  main window change in any other way than by `Apply corrections to the
  data`. Session files of earlier versions load as before; a session saved
  with 1.0.2 also loads in earlier versions, which leave the kept run
  aside.
- MCR-ALS window: the export button `C` (concentration profiles) ended
  with an error message, because the header of the file had one column
  name too many. Fixed.
- Main window before data are loaded, `Split` view: the upper panel of the
  2D map no longer carries numbers and a label of its own on the spectral
  axis. The axis is labeled once, under the lower panel, as when data are
  shown.
- Texts: buttons, tooltips, messages and reports of the program, and the
  English manual, use American spelling throughout (for example
  `neighbors`, `center`, `Not analyzed`).
- Python version: the repository now holds a Python version of the program
  in the folder `python`. It has the functions of the MATLAB version
  1.0.2, runs without MATLAB and carries the same version number.
  `python/README.txt` describes it. The manuals name it in the
  introduction and in a new part `Python version` of section 2.

## 1.0.1 (2026-11-22)

Display fixes. No calculation is changed: results are the same as with
1.0.0.

- `Reset all`: the status text no longer keeps `[cropped]` or the message
  of an earlier time-zero shift.
- Crop window: the numbers of kept points in the title of the preview are
  counted on the data that are being cropped. They were wrong after an
  earlier crop or time-zero shift.
- Main window: a wavelength or delay typed into its edit field now also
  moves the marker line in the other small plot (spectrum or kinetics).
- Crop and Chirp windows: the white limit lines and the white `+` of the
  saved points have a dark outline, so that they show on the white center
  of the BWR, RdBu and seismic colormaps.
- Fit trace window: changing `Components` keeps the values typed into the
  table. They used to go back to those of the last run.

## 1.0.0 (2026-11-21)

First public version.
