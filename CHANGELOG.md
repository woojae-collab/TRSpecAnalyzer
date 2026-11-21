# Changelog

Version numbers are major.minor.patch:

- patch (1.0.x): bug fixes and changes of wording
- minor (1.x.0): new features
- major (x.0.0): changes that break older session files or the way the
  program is used

The version of a copy of the program is shown in the title of its main
window and returned by `TRSpecAnalyzer.Version`.

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
  saved points have a dark outline, so that they show on the white centre
  of the BWR, RdBu and seismic colormaps.
- Fit trace window: changing `Components` keeps the values typed into the
  table. They used to go back to those of the last run.

## 1.0.0 (2026-11-21)

First public version.
