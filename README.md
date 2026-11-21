# TRSpecAnalyzer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24420778.svg)](https://doi.org/10.5281/zenodo.24420778)

A MATLAB app for analysing time-resolved spectroscopy data: transient
absorption (fs/ns TA), time-resolved photoluminescence, time-resolved
Raman and IR, or any other signal measured as a delay x wavelength matrix.

Loading, corrections, fitting and decomposition are done in one program,
which is a single file, `TRSpecAnalyzer.m`.

![Synthetic demo data set with its true spectra and populations](demo/demo_preview.png)

*The first demo data set (`demo/demo_TA.csv`, synthetic) and the answer the
analysis should recover.*

## Requirements

- MATLAB R2020b or later
- No additional toolboxes

## Getting started

1. Download `TRSpecAnalyzer.m` (or clone this repository) and put it in a
   folder on the MATLAB path.
2. In MATLAB, run

   ```matlab
   TRSpecAnalyzer
   ```

3. Press `Load...` and open `demo/demo_TA.csv`. Then either follow the
   walkthrough in [`demo/README_demo.txt`](demo/README_demo.txt), or press
   `Auto Analysis...` and `Run` for a first overview.

## Data format

Text files (CSV / DAT / TXT; comma, tab or space delimited):

```
corner , t1 , t2 , ... , tN        first row    : delay times
l1     , signal matrix (M x N)     first column : wavelength or wavenumber
...
lM     ,
```

Several scans of the same measurement can be loaded together, inspected
and averaged.

## What it does

**Data handling**

- Load one file, or load and average several scans with per-file quality
  check, outlier detection and time-shift alignment
- Background subtraction, time-zero shift, chirp correction, cropping, and
  masking of wavelength ranges (for example pump scatter)
- Time unit fs / ps / ns / us; spectral axis in nm or cm-1
- Sessions: the whole working state, including every result, is saved to
  and restored from one `.mat` file

**Analysis**

- Single-wavelength fit: sum of exponentials (optionally stretched)
  convolved with a Gaussian IRF
- SVD
- Global Analysis: DADS and EADS, with 1-sigma uncertainties of the fitted
  parameters and optional coherent-artifact terms at time zero
- Target Analysis: a compartment model of your own, written as text
  (`A -> B`, `A <-> B : kf, kb`, `B -> GS`), with shared, tied, fixed or
  limited rate constants and spectral constraints (zero ranges, equal
  spectra, non-negative spectra). Parameters that the data do not
  determine are flagged
- Lifetime Density Analysis (LDA)
- MCR-ALS
- Coherence analysis of the residual oscillations and time-frequency maps
- Auto Analysis: proposes the corrections and the number of components and
  writes a report that lists every decision with the numbers behind it

Corrected data and analysis results can be exported as CSV files.

## Documentation

| | PDF | Word |
|---|---|---|
| User manual (English) | [PDF](docs/TRSpecAnalyzer_User_Manual_EN.pdf) | [docx](docs/TRSpecAnalyzer_User_Manual_EN.docx) |
| 사용 설명서 (한국어) | [PDF](docs/TRSpecAnalyzer_User_Manual_KO.pdf) | [docx](docs/TRSpecAnalyzer_User_Manual_KO.docx) |

## Demo data

The [`demo`](demo) folder holds three synthetic data sets with known
answers, the true spectra and populations they were made from, and the
script that generated them (`make_demo_data.m`). Nothing in it is a
measurement.

| File | Kinetics | Used for |
|---|---|---|
| `demo_TA.csv` (+ `scans/`) | A -> B -> C -> ground state, with chirp, background, pump scatter and oscillations | the full walkthrough, Load & Average, Auto Analysis |
| `demo_target.csv` | A -> B and A -> C in parallel, B -> ground state, C -> ground state | Target Analysis of a model that a sequential fit cannot describe |
| `demo_equilibrium.csv` | A <-> B, B -> C -> ground state | Target Analysis with a step that runs in both directions |

`demo/README_demo.txt` gives the steps and the numbers to expect at each
one, so the program can be checked on your own installation.

## Versions

The version number (major.minor.patch) is shown in the title of the main
window and returned by `TRSpecAnalyzer.Version`. It is also written into
saved sessions and Auto Analysis reports. Changes are listed in
[CHANGELOG.md](CHANGELOG.md), and every version stays available under
[Releases](https://github.com/woojae-collab/TRSpecAnalyzer/releases).
When you report results, please state the version you used.

## 한국어 안내

TRSpecAnalyzer는 시간분해 분광 데이터(TA, TR-PL, TR-Raman, TR-IR)를
불러와 보정하고, 피팅과 분해 분석까지 한 프로그램에서 수행하는 MATLAB
도구입니다. `TRSpecAnalyzer.m` 파일 하나만 있으면 되고, MATLAB R2020b
이상에서 추가 toolbox 없이 실행됩니다. 자세한 사용법은
[한국어 사용 설명서](docs/TRSpecAnalyzer_User_Manual_KO.pdf)를 보세요.

## Feedback

Questions and bug reports are welcome under *Issues*. For a bug report,
the text of the error dialog and, if possible, a small data file that
reproduces the problem help most.

## How to cite

If TRSpecAnalyzer was useful for your work, please cite the version you
used. For v1.0.0:

> Kim, W. (2026). TRSpecAnalyzer: a MATLAB app for the analysis of
> time-resolved spectroscopy data (v1.0.0). Zenodo.
> https://doi.org/10.5281/zenodo.24420778

## Author

Woojae Kim, ULS Lab, Department of Chemistry, Yonsei University.
The code was written with Claude Code.

## License

MIT License. See [LICENSE](LICENSE).
