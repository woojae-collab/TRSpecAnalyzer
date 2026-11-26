# References

TRSpecAnalyzer brings established methods for the analysis of
time-resolved spectra together in one program. This file names, part by
part, the publications these methods rest on. Where the program goes its
own way (for example in the automatic analysis, or in the rule for the
interpolation along the delays), the manual and
[`CHANGELOG.md`](CHANGELOG.md) describe what it does.

If a result of one of these analyses goes into a publication, please
consider citing the papers of that method beside TRSpecAnalyzer itself
(see *How to cite* in [`README.md`](README.md)).

## Global analysis, target analysis and fits of single traces

The model: a sum of exponential decays convolved with a Gaussian
instrument response; decay-associated and evolution-associated spectra
(DADS, EADS); compartment schemes with species-associated spectra (SAS),
including constraints on the spectra.

1. I. H. M. van Stokkum, D. S. Larsen, R. van Grondelle, "Global and
   target analysis of time-resolved spectra," *Biochim. Biophys. Acta -
   Bioenergetics* **1657**, 82-104 (2004).
   https://doi.org/10.1016/j.bbabio.2004.04.011
2. A. R. Holzwarth, "Data analysis of time-resolved measurements," in
   *Biophysical Techniques in Photosynthesis*, edited by J. Amesz and
   A. J. Hoff, Advances in Photosynthesis and Respiration, vol. 3
   (Kluwer Academic Publishers, Dordrecht, 1996), pp. 75-92.
   https://doi.org/10.1007/0-306-47960-5_5

The spectra enter the model linearly and are solved for at every step of
the fit, so that only the time constants and the parameters of the
instrument response are searched (separable least squares, variable
projection).

3. G. H. Golub, V. Pereyra, "The differentiation of pseudo-inverses and
   nonlinear least squares problems whose variables separate," *SIAM J.
   Numer. Anal.* **10**, 413-432 (1973).
   https://doi.org/10.1137/0710036
4. K. M. Mullen, I. H. M. van Stokkum, "The variable projection algorithm
   in time-resolved spectroscopy, microscopy and mass spectrometry
   applications," *Numer. Algorithms* **51**, 319-340 (2009).
   https://doi.org/10.1007/s11075-008-9235-2

The search itself is the simplex method of Nelder and Mead in the form
of MATLAB's `fminsearch` (the Python version uses the same algorithm).

5. J. A. Nelder, R. Mead, "A simplex method for function minimization,"
   *Comput. J.* **7**, 308-313 (1965).
   https://doi.org/10.1093/comjnl/7.4.308
6. J. C. Lagarias, J. A. Reeds, M. H. Wright, P. E. Wright, "Convergence
   properties of the Nelder-Mead simplex method in low dimensions," *SIAM
   J. Optim.* **9**, 112-147 (1998).
   https://doi.org/10.1137/S1052623496303470

Uncertainties of the fitted parameters come from the linear
approximation of the model at the optimum. Non-negative spectra in
target analysis are found by non-negative least squares.

7. D. M. Bates, D. G. Watts, *Nonlinear Regression Analysis and Its
   Applications* (Wiley, New York, 1988).
8. C. L. Lawson, R. J. Hanson, *Solving Least Squares Problems*, Classics
   in Applied Mathematics 15 (SIAM, Philadelphia, 1995; first published
   by Prentice-Hall, 1974).

Stretched exponential decays, and the quadrature rule used to convolve
them with the instrument response.

9. R. Kohlrausch, *Ann. Phys. Chem. (Poggendorff)* **91**, 56 and 179
   (1854).
10. G. Williams, D. C. Watts, "Non-symmetrical dielectric relaxation
    behaviour arising from a simple empirical decay function," *Trans.
    Faraday Soc.* **66**, 80-85 (1970).
    https://doi.org/10.1039/TF9706600080
11. G. H. Golub, J. H. Welsch, "Calculation of Gauss quadrature rules,"
    *Math. Comp.* **23**, 221-230 (1969).
    https://doi.org/10.2307/2004418

Glotaran, a widely used program for global and target analysis; the terms
DADS, EADS and SAS are used here as they are used there.

12. J. J. Snellenburg, S. P. Laptenok, R. Seger, K. M. Mullen,
    I. H. M. van Stokkum, "Glotaran: a Java-based graphical user
    interface for the R package TIMP," *J. Stat. Softw.* **49**(3), 1-22
    (2012). https://doi.org/10.18637/jss.v049.i03

## Chirp of the probe and signals at time zero

Time zero that depends on the probe wavelength (chirp of a
supercontinuum probe) and its correction; signals that exist only while
pump and probe overlap (cross-phase modulation, two-photon absorption,
stimulated Raman amplification), which the fits describe by a Gaussian
and its first and second time derivative (see also reference 1).

13. S. A. Kovalenko, A. L. Dobryakov, J. Ruthmann, N. P. Ernsting,
    "Femtosecond spectroscopy of condensed phases with chirped
    supercontinuum probing," *Phys. Rev. A* **59**, 2369-2384 (1999).
    https://doi.org/10.1103/PhysRevA.59.2369
14. M. Lorenc, M. Ziolek, R. Naskrecki, J. Karolczak, J. Kubicki,
    A. Maciejewski, "Artifacts in femtosecond transient absorption
    spectroscopy," *Appl. Phys. B* **74**, 19-27 (2002).
    https://doi.org/10.1007/s003400100750

## Singular value decomposition

15. E. R. Henry, J. Hofrichter, "Singular value decomposition:
    application to analysis of experimental data," *Methods Enzymol.*
    **210**, 129-192 (1992).

## Lifetime density analysis

The data are described by a continuous distribution of lifetimes on a
fixed grid; the amplitudes are found with Tikhonov regularization, and
the regularization parameter is proposed at the corner of the L-curve
(the point of largest curvature).

16. G. F. Dorlhiac, C. Fare, J. J. van Thor, "PyLDM - An open source
    package for lifetime density analysis of time-resolved spectroscopic
    data," *PLoS Comput. Biol.* **13**(5), e1005528 (2017).
    https://doi.org/10.1371/journal.pcbi.1005528
17. C. Slavov, H. Hartmann, J. Wachtveitl, "Implementation and evaluation
    of data analysis strategies for time-resolved optical spectroscopy,"
    *Anal. Chem.* **87**, 2328-2336 (2015).
    https://doi.org/10.1021/ac504348h
18. A. N. Tikhonov, V. Y. Arsenin, *Solutions of Ill-Posed Problems*
    (Winston, Washington, D.C., 1977).
19. P. C. Hansen, "Analysis of discrete ill-posed problems by means of
    the L-curve," *SIAM Rev.* **34**, 561-580 (1992).
    https://doi.org/10.1137/1034115
20. P. C. Hansen, D. P. O'Leary, "The use of the L-curve in the
    regularization of discrete ill-posed problems," *SIAM J. Sci.
    Comput.* **14**, 1487-1503 (1993). https://doi.org/10.1137/0914086

## Multivariate curve resolution (MCR-ALS)

Alternating least squares with non-negativity and unimodality
constraints; the lack of fit (LOF) as the measure of convergence.

21. R. Tauler, "Multivariate curve resolution applied to second order
    data," *Chemom. Intell. Lab. Syst.* **30**, 133-146 (1995).
    https://doi.org/10.1016/0169-7439(95)00047-X
22. A. de Juan, R. Tauler, "Multivariate curve resolution (MCR) from
    2000: progress in concepts and applications," *Crit. Rev. Anal.
    Chem.* **36**, 163-176 (2006).
    https://doi.org/10.1080/10408340600970005
23. J. Jaumot, R. Gargallo, A. de Juan, R. Tauler, "A graphical
    user-friendly interface for MCR-ALS: a new tool for multivariate
    curve resolution in MATLAB," *Chemom. Intell. Lab. Syst.* **76**,
    101-110 (2005). https://doi.org/10.1016/j.chemolab.2004.12.007
24. A. de Juan, Y. Vander Heyden, R. Tauler, D. L. Massart, "Assessment
    of new constraints applied to the alternating least squares method,"
    *Anal. Chim. Acta* **346**, 307-318 (1997).
25. C. Ruckebusch, M. Sliwa, P. Pernot, A. de Juan, R. Tauler,
    "Comprehensive data analysis of femtosecond transient absorption
    spectra: a review," *J. Photochem. Photobiol. C* **13**, 1-27 (2012).
    https://doi.org/10.1016/j.jphotochemrev.2011.10.002

## Vibrational coherence

Fourier transform of the residual that a kinetic fit leaves; the
subtraction of the coherence of the solvent, measured alone (thesis,
Fig. 4.6); window functions; the Morlet wavelet transform of the
time-frequency window and its cone of influence.

26. M. Liebel, C. Schnedermann, T. Wende, P. Kukura, "Principles and
    applications of broadband impulsive vibrational spectroscopy,"
    *J. Phys. Chem. A* **119**, 9506-9517 (2015).
    https://doi.org/10.1021/acs.jpca.5b05948
27. M. Liebel, *Understanding Molecular Dynamics with Coherent
    Vibrational Spectroscopy in the Time-domain*, DPhil thesis
    (University of Oxford, 2014).
    https://ora.ox.ac.uk/objects/uuid:e0289d80-f6e3-4e6f-817e-f8dd55d15bc4
28. F. J. Harris, "On the use of windows for harmonic analysis with the
    discrete Fourier transform," *Proc. IEEE* **66**, 51-83 (1978).
    https://doi.org/10.1109/PROC.1978.10837
29. C. Torrence, G. P. Compo, "A practical guide to wavelet analysis,"
    *Bull. Am. Meteorol. Soc.* **79**, 61-78 (1998).
    https://doi.org/10.1175/1520-0477(1998)079%3C0061:APGTWA%3E2.0.CO;2

## Anisotropy

Anisotropy and magic-angle signal from parallel and perpendicular
measurements; the factor G and its determination by tail matching.

30. J. R. Lakowicz, *Principles of Fluorescence Spectroscopy*, 3rd ed.
    (Springer, New York, 2006).
    https://doi.org/10.1007/978-0-387-46312-4

## Position of a band along the delays

The position of a band from its first moment, its maximum or a
log-normal fit; the normalized shift C(t) of solvation dynamics.

31. M. Maroncelli, G. R. Fleming, "Picosecond solvation dynamics of
    coumarin 153: the importance of molecular aspects of solvation,"
    *J. Chem. Phys.* **86**, 6221-6239 (1987).
    https://doi.org/10.1063/1.452460
32. M. L. Horng, J. A. Gardecki, A. Papazyan, M. Maroncelli,
    "Subpicosecond measurements of polar solvation dynamics: coumarin
    153 revisited," *J. Phys. Chem.* **99**, 17311-17337 (1995).
    https://doi.org/10.1021/j100048a004
33. D. B. Siano, D. E. Metzler, "Band shapes of the electronic spectra of
    complex molecules," *J. Chem. Phys.* **51**, 1856-1861 (1969).
    https://doi.org/10.1063/1.1672270

## Interpolation along the delays

Cubic splines with the not-a-knot end condition (chirp correction and
time shifts).

34. C. de Boor, *A Practical Guide to Splines*, Applied Mathematical
    Sciences 27 (Springer, New York, 1978; revised edition 2001).

## Software the Python version is built on

35. C. R. Harris, K. J. Millman, S. J. van der Walt, et al., "Array
    programming with NumPy," *Nature* **585**, 357-362 (2020).
    https://doi.org/10.1038/s41586-020-2649-2
36. P. Virtanen, R. Gommers, T. E. Oliphant, et al. (SciPy 1.0
    Contributors), "SciPy 1.0: fundamental algorithms for scientific
    computing in Python," *Nat. Methods* **17**, 261-272 (2020).
    https://doi.org/10.1038/s41592-019-0686-2
37. J. D. Hunter, "Matplotlib: a 2D graphics environment," *Comput. Sci.
    Eng.* **9**(3), 90-95 (2007). https://doi.org/10.1109/MCSE.2007.55
