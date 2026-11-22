# Changelog

Version numbers are major.minor.patch:

- patch (1.0.x): bug fixes and changes of wording
- minor (1.x.0): new features
- major (x.0.0): changes that break older session files or the way the
  program is used

The version of a copy of the program is shown in the title of its main
window and returned by `TRSpecAnalyzer.Version`.

## 1.0.2 (2026-11-22)

Auto Analysis keeps its last run, as the other analysis windows do; two
small fixes; American spelling in the texts. No calculation is changed:
results are the same as with 1.0.1.

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
