%% 
classdef TRSpecAnalyzer < handle
    %TRSPECANALYZER  Time-Resolved Spectroscopy data analysis GUI.
    %
    %   A general-purpose tool for any 2D time-resolved data
    %   (delay x wavelength), including but not limited to:
    %     * Transient Absorption (TA, fsTA, nsTA)
    %     * Time-Resolved Fluorescence (TRPL)
    %     * Time-Resolved Raman / TR-IR
    %
    %   Usage:
    %       >> TRSpecAnalyzer
    %
    %   Data file format (CSV / DAT / TXT, comma/tab/space delimited):
    %       First row  : [corner, t1, ..., tN]      delay time
    %       First col  : [corner, l1, ..., lM]      wavelength / energy
    %       Interior   : signal matrix (M x N)
    %
    %   Features:
    %     Data handling
    %     - Single- or multi-file load (multi-file -> Load & Average window)
    %     - Load & Average window with per-file QC, outlier detection, a
    %       "Show legend" toggle (compact Scan-N labels), click+drag Ref
    %       lambda / Ref t markers on the preview plots, and per-file time
    %       shift with manual entry or auto-alignment (cross-correlation
    %       to a user-selected reference) before averaging
    %     - Crop / Background / t=0 correction tools. Crop bakes all
    %       prior corrections into the new baseline (so e.g. the NaN edge
    %       left by chirp correction really disappears); "Reset all"
    %       returns to the as-loaded raw data
    %     - Chirp correction with a click+drag candidate crosshair and an
    %       explicit "Save point" commit (>=4 points -> fit), plus a ridge
    %       auto-suggest helper (max|dA/dt| / max|dA| / zero-cross)
    %     - Mask wavelength regions (e.g. pump scattering)
    %     - Time-axis unit switching (fs / ps / ns / us) with unit-aware
    %       defaults (e.g. chirp window, GA tau bounds)
    %     - Spectral-axis unit switching: wavelength (nm, default) or
    %       wavenumber (cm-1) for TR-Raman / TR-IR; relabels every axis,
    %       title, control and export header (values are not converted)
    %     - Session save / load: the full working state (data, all
    %       corrections, fit settings, every analysis result and its
    %       settings, view options) round-trips through a .mat file so an
    %       analysis can be resumed across MATLAB restarts
    %     - Preferences persist across MATLAB restarts (prefdir): last
    %       folder used in any file dialog, time unit, spectral unit,
    %       colormap, delay scale mode and split threshold are restored
    %       at startup
    %
    %     Visualization
    %     - 2D map with drag-based crosshair: click and drag to smoothly
    %       move the monitored (lambda, t), with synchronized vertical
    %       lines on the spectrum and kinetics side panels
    %     - Linear / Log / Split scale on the delay axis; split-view
    %       kinetics with synchronized y-range and decade ticks
    %     - Spectrum & kinetics side panels with pin overlay
    %     - Multiple colormaps with consistent crosshair styling across
    %       all child windows
    %
    %     Kinetic / spectral analysis
    %       (every analysis window supports the same click+drag
    %        inspection pattern as the main window, and keeps its result
    %        and settings when closed and reopened. Any change to the
    %        working data -- new load, crop, reset, background / chirp /
    %        mask correction, t=0 shift, or a unit switch -- closes the
    %        analysis windows and discards their stored results, so a
    %        result is never shown for data it was not computed from)
    %     - Single-wavelength kinetic fit (multi-exp + stretched-exp +
    %       IRF), with an IRF-free "None" mode and a fit-time-range
    %       restriction (fit only t within [t_min, t_max] without
    %       cropping; display zooms to the fitted range)
    %     - Singular Value Decomposition (SVD) explorer
    %     - Global Analysis (DADS / EADS, sequential model; with
    %       stretched components the SAS and populations come from a
    %       population-conserving sequential chain in which any step may
    %       be stretched -- with a choice of whether the clock of a
    %       stretched step starts when the population enters that state
    %       or at the excitation), with the same IRF-free mode
    %       and fit-time-range restriction (maps and kinetics zoom to the
    %       fitted range), an EADS population-kinetics panel c_i(t) with
    %       CSV export, map/kinetics time-scale sync, and an "SVD?"
    %       button that estimates the number of components from the
    %       singular values above the noise floor
    %     - Target analysis: fit a compartment model of your own. The
    %       scheme is written as text (A -> B, A <-> B, A -> GS, B -> 2 C)
    %       or edited in a K-matrix table, and drawn as a box-and-arrow
    %       diagram (rates on the arrows; after a fit also the time
    %       constants and branching fractions; a click enlarges it in a
    %       window of its own); compartment names are free
    %       and generic presets are provided (sequential, parallel,
    %       branching, reversible, plus a few worked examples).
    %       First-order steps. Rate parameters can be named, shared
    %       between steps, tied by a factor (2*k1), fixed, or limited.
    %       Spectral constraints: zero ranges, equal / proportional
    %       spectra, non-negative spectra. Gaussian IRF or IRF-free,
    %       fit-time-range restriction. Shows SAS or DAS of the same
    %       fit, populations, kinetics with the contribution of every
    %       compartment, the residual map and the SVD of the residual;
    %       reports observable lifetimes and branching fractions. The
    %       probe position of the kinetics panel is a click+drag marker
    %       on the spectra panel and on the residual map
    %     - Fit parameter uncertainties: GA, target and single-wavelength
    %       fits report 1-sigma errors (numerical-Jacobian covariance)
    %       next to each free tau / beta / t0 / FWHM. A parameter the
    %       data do not determine (e.g. two rate constants that only act
    %       through their sum) is flagged instead of given a number
    %     - IRF model: "Numerical conv." includes the Gaussian IRF (plain
    %       exponentials by the closed-form exp (x) Gaussian, stretched
    %       exponentials by per-delay Gauss-Legendre quadrature); "None"
    %       fits pure exponentials without an IRF
    %     - Fitted IRF width (GA, target and single-wavelength fits): a
    %       FWHM close to or below half the delay step around time zero
    %       is not resolved by the delays. A fit that ends there is
    %       repeated in two steps (IRF fixed, then free, with the FWHM
    %       kept above that width) and the better result is shown; a
    %       width that stays there is marked "[not resolved]"
    %     - Coherent-artifact terms (Global Analysis, Target analysis,
    %       Auto Analysis): a signal that exists only while pump and
    %       probe overlap is fitted with up to three terms at time zero
    %       -- a Gaussian and its first and second time derivative, each
    %       with a free spectrum. Their width is that of the IRF or one
    %       of their own (fitted or fixed). The artifact spectra are
    %       drawn dashed, exported with the DADS / SAS, and are not part
    %       of the DADS / EADS / SAS
    %     - Lifetime Density Analysis (LDA, Tikhonov L2 / 1st-derivative
    %       smoothness, automatic L-curve corner finder, click+drag
    %       inspection), with an optional no-IRF basis and a fit-time-
    %       range restriction (mind the tau-grid lower bound when
    %       restricting; the L-curve uses the same window)
    %     - Multivariate Curve Resolution by ALS (MCR-ALS, non-negativity
    %       and unimodality constraints, SVD/random/GA-based init -- the
    %       random start takes a new seed every run, reports it, and can
    %       be repeated from it; the GA start uses the EADS / SAS or the
    %       DADS of the stored global fit -- log / linear time scale,
    %       wavelength marker on S spectra)
    %
    %     Vibrational coherence
    %     - Coherence map: per-pixel FFT of the kinetic residual (from
    %       GA fit or per-pixel exp fit), with selectable time window,
    %       apodization, zero-padding, detrending, and wavenumber/THz
    %       output. Heatmap z-axis (and 1D y-axis) supports linear/log; a
    %       display nu range zooms the frequency axes and renormalizes
    %       the colour scale to the visible band
    %     - Click+drag inspection on the |FFT|^2 map updates four linked
    %       plots: power spectrum at the clicked lambda; power spectrum
    %       averaged over a user-defined lambda range (dashed bounds on
    %       the heatmap); spectral profile at the clicked wavenumber,
    %       split into |FFT|^2 and phase (masked below 5% of peak);
    %       residual time trace at the clicked lambda
    %     - CSV exports: the full map (frequency on rows, wavelength on
    %       columns), P(nu) at lambda, the lambda-averaged P(nu), and the
    %       profile at nu including a phase column
    %     - Time-frequency analysis window (opened from Coherence): nu-t
    %       map of the coherence residual at one lambda or a lambda-range
    %       average by STFT (sliding-window FFT) or Morlet CWT (FFT-based,
    %       no Wavelet Toolbox; peak-response normalisation so equal-
    %       amplitude modes have equal power, centre-frequency scale
    %       mapping so peaks sit at the true nu). Shows the cone of
    %       influence, power-vs-time at a chosen nu (dephasing), the
    %       instantaneous spectrum at a chosen t, and the input trace with
    %       the effective window; click+drag cursor and CSV export
    %
    %     Automatic analysis
    %     - Auto Analysis window: a report on the data as they are shown.
    %       Proposes the corrections -- background from the delays found
    %       to lie before the signal, mask for channels with excess
    %       noise (only where the signal is lost in it, or every noisy
    %       channel), chirp curve through the steepest rise of every
    %       wavelength, crop of what the corrections leave empty -- and
    %       summarizes the kinetics without a kinetic model: singular
    %       values, global fits with 1, 2, ... components (several
    %       starting values each) and the number the data support, time
    %       constants with uncertainties, a component that does not
    %       decay, DADS / EADS, kinetics at a few wavelengths (data and
    %       fit), and oscillations left in the residual.
    %       Every decision is listed with the numbers behind it, plus the
    %       points to check by eye. Coherent-artifact terms are added
    %       when the data call for them (a time constant far below the
    %       IRF width, or a clear gain in RMS). The automatic choice of
    %       the number of components can be overruled: every number that
    %       was fitted can be picked ("Components used"). The data are
    %       changed only by "Apply corrections"; the model in use can be
    %       sent to Global Analysis
    %
    %   Author: Woojae Kim (written with Claude Code)
    %           ULS Lab, Department of Chemistry, Yonsei University
    %   Last updated: 2026-11-21
    %
    %   Requires MATLAB R2020b or later. No additional toolboxes.
    
    % =====================================================================
    properties (Access = public)
        % ---- Data (current working) ----
        wavelength  = []
        delay       = []
        deltaA      = []
        deltaA_raw  = []
        
        % ---- Original (snapshot of loaded data; used for crop/revert) ----
        originalWavelength = []
        originalDelay      = []
        originalDeltaA     = []
        dataSourceDesc     = ''
        
        % ---- Correction state ----
        bgApplied    = false
        bgSpectrum   = []
        bgN          = 0
        
        chirpApplied = false
        chirpPts     = zeros(0,2)
        chirpParams  = []            % [a, b, c, d] for t0 = a*sqrt((b*w^2-1)/(c*w^2-1)) + d
        chirpFitRMS  = NaN
        chirpViewTMin = []           % display window on axMap (defaults to first delay)
        chirpViewTMax = 2            % default 2 ps
        
        % ---- Ridge helper ----
        ridgeMethod  = 'max|dA/dt|'
        ridgeTmin    = -0.5
        ridgeTmax    =  1.5
        ridgeSmoothN = 5
        ridgeWL      = []
        ridgeT       = []
        
        % ---- Current selection ----
        selWL = []
        selT  = []
        
        % ---- Overlays ----
        specOverlays = zeros(0,1)
        kinOverlays  = zeros(0,1)
        
        % ---- Delay-axis scale mode ----
        delayScaleMode = 'split'
        splitThreshold = 1.0
        
        % ---- Main 2D map display range ----
        mainViewTMin = []
        mainViewTMax = []
        
        % ---- 2D map color settings ----
        mapColormap = 'BWR'
        mapZMin     = []
        mapZMax     = []
        
        % ---- Main GUI handles ----
        fig
        pnl2D; pnlSpec; pnlKin
        axCont2D; axContSpec; axContKin
        ax2D       = gobjects(0)
        axSpec     = gobjects(0)
        axKin      = gobjects(0)
        hVLine     = gobjects(0)
        hHLine     = gobjects(0)
        cbarAx     = []   % dedicated axes to host colorbar in split mode
        
        edSelWL; edSelT
        ddScale; edThreshold
        edViewTMin; edViewTMax
        ddColormap; edZMin; edZMax
        lblStatus
        
        % ---- Child windows ----
        bgFig    = []
        chirpFig = []
        cropFig  = []
        loadFig  = []
        gaFig    = []
        maskFig  = []
        kinFitFig = []
        svdFig    = []
        ldaFig    = []
        targetFig = []         % target analysis (kinetic scheme) window
        mcrFig    = []
        coherenceFig = []
        tfFig        = []      % time-frequency (STFT / CWT) window
        autoFig      = []      % automatic analysis window
        
        % ---- Main 2D map drag state ----
        dragging2D = false
        drag2DAx   = []
        
        % ---- Single-trace kinetic fit settings (persist across opens) ----
        kfNComp     = 3
        kfTauInit   = [1 10 100]
        kfTauFixed  = [false false false]
        kfBetaInit  = [1 1 1]            % stretched exponent per component
        kfBetaFixed = [true true true]   % default fixed=1 -> ordinary exp
        kfStretchOn = [false false false]% per-component stretched ON/OFF
        kfHasInf    = false
        kfT0        = 0
        kfT0Fixed   = true
        kfFWHM      = 0.15
        kfFWHMFixed = true
        kfIRFTreat  = 'numerical'       % IRF model: 'numerical' | 'none'
        
        % ---- Zero-time shift (for data without chirp, e.g. ns-TA) ----
        tZeroShift = 0   % cumulative applied shift; delay = delay_raw - this
        
        % ---- Wavelength region masks (e.g. for pump scattering) ----
        % Nx3 cell array: {wlMin, wlMax, 'nan'|'zero'}. Applied in
        % recompute after BG + chirp corrections. Persists through
        % unrelated corrections; cleared by Reset Corrections, Load,
        % and Crop.
        maskedRegions = cell(0, 3)
        
        % ---- Global analysis state ----
        gaNComp       = 3
        gaTauInit     = [1 10 100]
        gaTauFixed    = [false false false]
        gaBetaInit    = [1 1 1]
        gaStretchOn   = [false false false]
        gaIRFTreat    = 'numerical'     % IRF model: 'numerical' | 'none'
        gaCA          = []              % coherent-artifact terms of the fit
                                        % (see caOptions): order, mode, fwhm
        gaResultCA    = []              % artifact part of the stored result:
                                        % struct order, spec (M x order),
                                        % fwhm, fwhmStd, mode; [] = none
        gaSeqModel    = 'age'           % stretched steps in the SAS: clock
                                        % from entry ('age') or from t0 ('t0')
        gaHasInf      = false
        gaT0          = 0
        gaT0Fixed     = true
        gaFWHM        = 0.15
        gaFWHMFixed   = true
        
        % ---- Global analysis results ----
        gaResultTau      = []
        gaResultBeta     = []
        gaResultT0       = []
        gaResultFWHM     = []
        gaResultHasInf   = false
        gaResultDADS      = []      % M x Nc amplitudes (wavelength x component)
        gaResultEADS     = []      % M x Nc sequential EADS (sorted by tau asc)
        gaResultFit      = []      % M x N fitted data
        gaResultRMS      = NaN
        gaResultTauStd   = []      % 1-sigma per tau (NaN where fixed)
        gaResultStretchOn = []     % stretched flags used in the stored fit
        gaResultIRF      = ''      % IRF model used in the stored fit
        gaResultTRange   = []      % fit time range used in it ([] = full)
        gaFitTMin = []             % GA fit time range (empty = full)
        gaFitTMax = []
        kfFitTMin = []             % single-wl fit time range
        kfFitTMax = []
        
        % Persistent stores for child analysis windows, so results
        % survive the window being closed and reopened. Each holds a
        % struct with both the settings used and the computed result
        % (empty [] = nothing run yet).
        cohSaved = []      % vibrational coherence
        ldaSaved = []      % lifetime density analysis
        mcrSaved = []      % MCR-ALS
        svdSaved = []      % SVD
        kfSaved  = []      % single-wavelength kinetic fit
        tfSaved  = []      % time-frequency (STFT / CWT) analysis
        tgSaved  = []      % target analysis: last result
        tgSet    = []      % target analysis: scheme and fit settings
        
        % Persistent app preferences (stored in prefdir, survive restarts)
        prefLastDir = ''   % last folder used in a file dialog
        
        % true while a session is being restored: recompute() must not
        % treat the freshly restored data as "changed" and wipe the
        % analysis results that were just loaded with it.
        holdDataChange = false
        quadGL = []                % cached Gauss-Legendre nodes / weights
        
        % Spectral-axis unit: 'nm' (wavelength) or 'cm-1' (wavenumber,
        % for TR-Raman / TR-IR). Display + export labels only; the
        % loaded numeric values are never converted.
        specUnit = 'nm'
        ddSpecUnit                    % dropdown handle (main toolbar)
        lblSelWL                      % 'lambda (nm):' label, main kinetics
        
        % ---- Time-axis unit ----
        timeUnit             = 'ps'   % 'fs' | 'ps' | 'ns' | 'us'
        ddTimeUnit                    % dropdown handle (main toolbar)
        lblSplitThreshold             % toolbar label handles whose Text
        lblViewTMin                   % is unit-dependent -- refreshed on
        lblViewTMax                   % timeUnit change so static UI text
        lblSelT                       % stays consistent with axes
    end
    
    % =====================================================================
    methods
        % ========== Constructor ==========
        function app = TRSpecAnalyzer()
            app.createMainGUI();
        end
        
        % =================================================================
        % GUI construction
        % =================================================================
        function createMainGUI(app)
            % Stored preferences first (properties only), so every
            % control and axes below is created with them.
            app.loadPrefs();
            
            app.fig = uifigure('Name', 'TR Spectroscopy Analyzer', ...
                'Position', [80 60 1400 820], ...
                'CloseRequestFcn', @(~,~) app.onClose());
            app.applyLightTheme(app.fig);
            % Drag handlers are permanently attached; they no-op until a
            % click on the 2D map flips the dragging flag.
            app.fig.WindowButtonMotionFcn = @(~,~) app.on2DMotion();
            app.fig.WindowButtonUpFcn     = @(~,~) app.on2DUp();
            
            gl = uigridlayout(app.fig, [5, 2]);
            gl.RowHeight     = {38, 38, 38, '1x', '1x'};
            gl.ColumnWidth   = {'2x', '1x'};
            gl.RowSpacing    = 6;
            gl.ColumnSpacing = 8;
            gl.Padding       = [10 10 10 10];
            
            % ---- Toolbar row 1: load / corrections / scale ----
            tb = uigridlayout(gl, [1, 13]);
            tb.Layout.Row = 1;  tb.Layout.Column = [1 2];
            tb.ColumnWidth = {80, 90, 80, 120, 110, 85, 130, ...
                              'fit', 85, 110, 60, 'fit', '1x'};
            tb.ColumnSpacing = 4;  tb.Padding = [0 0 0 0];
            
            uibutton(tb, 'Text', 'Load...', ...
                'ButtonPushedFcn', @(~,~) app.loadData());
            uibutton(tb, 'Text', 'Reset all', ...
                'Tooltip', ['Restore the raw data taken at load time. ' ...
                    'Undoes crop, background, chirp, mask, and t=0 ' ...
                    'shift -- equivalent to re-loading the file.'], ...
                'ButtonPushedFcn', @(~,~) app.resetCorrections());
            uibutton(tb, 'Text', 'Crop...', ...
                'ButtonPushedFcn', @(~,~) app.openCropWindow());
            uibutton(tb, 'Text', 'Background Corr...', ...
                'ButtonPushedFcn', @(~,~) app.openBackgroundWindow());
            uibutton(tb, 'Text', 'Chirp Corr...', ...
                'ButtonPushedFcn', @(~,~) app.openChirpWindow());
            uibutton(tb, 'Text', 'Set t=0', ...
                'Tooltip', ['Make the currently selected crosshair ' ...
                    'position on the main 2D map the new t=0. Click ' ...
                    'anywhere on the map first to move the crosshair ' ...
                    'to your desired zero-time location.'], ...
                'ButtonPushedFcn', @(~,~) app.applyZeroTimeAtCrosshair());
            uibutton(tb, 'Text', app.uiStr('Mask \lambda...'), ...
                'Tooltip', ['Mask wavelength regions (e.g., pump ' ...
                    'scattering). Can set them to NaN (excluded from ' ...
                    'fits) or to 0.'], ...
                'ButtonPushedFcn', @(~,~) app.openMaskWindow());
            
            uilabel(tb, 'Text', 'Scale:', ...
                'HorizontalAlignment', 'right');
            app.ddScale = uidropdown(tb, ...
                'Items', {'Linear', 'Log', 'Split'}, ...
                'Value', app.scaleModeLabel(), ...
                'ValueChangedFcn', @(src,~) app.onScaleModeChange(src.Value));
            
            app.lblSplitThreshold = uilabel(tb, 'Text', ...
                sprintf('Split thr. (%s):', app.tUnitTxt()), ...
                'HorizontalAlignment', 'right');
            app.edThreshold = uieditfield(tb, 'numeric', ...
                'Value', app.splitThreshold, ...
                'Limits', [1e-6 Inf], ...
                'ValueChangedFcn', @(src,~) app.onThresholdChange(src.Value));
            
            uilabel(tb, 'Text', '');   % spacer
            app.lblStatus = uilabel(tb, 'Text', 'No data loaded', ...
                'HorizontalAlignment', 'right', ...
                'FontColor', [0.3 0.3 0.3]);
            
            % ---- Toolbar row 2: analyses (SVD, GA, Target, LDA, MCR-ALS,
            %      Coherence, Auto Analysis) ----
            tbA = uigridlayout(gl, [1, 8]);
            tbA.Layout.Row = 2;  tbA.Layout.Column = [1 2];
            tbA.ColumnWidth   = {110, 160, 150, 110, 110, 130, 140, '1x'};
            tbA.ColumnSpacing = 6;  tbA.Padding = [0 0 0 0];
            
            uibutton(tbA, 'Text', 'SVD...', ...
                'BackgroundColor', [0.20 0.45 0.70], 'FontColor', 'w', ...
                'Tooltip', ['Singular Value Decomposition: shows ' ...
                    'singular-value spectrum, spectral and temporal ' ...
                    'basis vectors, and reconstruction with N ' ...
                    'components. Useful for choosing the number of ' ...
                    'exponentials in Global Analysis.'], ...
                'ButtonPushedFcn', @(~,~) app.openSvdWindow());
            uibutton(tbA, 'Text', 'Global Analysis...', ...
                'BackgroundColor', [0.55 0.30 0.70], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~,~) app.openGlobalAnalysisWindow());
            uibutton(tbA, 'Text', 'Target Analysis...', ...
                'BackgroundColor', [0.60 0.25 0.45], 'FontColor', 'w', ...
                'Tooltip', ['Target analysis: fit a compartment model ' ...
                    'of your own (any number of compartments; ' ...
                    'sequential, parallel, branching or reversible ' ...
                    'steps; shared or tied rate constants; spectral ' ...
                    'constraints) and get the species-associated ' ...
                    'spectra and populations it implies.'], ...
                'ButtonPushedFcn', @(~,~) app.openTargetWindow());
            uibutton(tbA, 'Text', 'LDA...', ...
                'BackgroundColor', [0.20 0.55 0.40], 'FontColor', 'w', ...
                'Tooltip', app.uiStr(['Lifetime Density Analysis: solve for ' ...
                    'a continuous distribution A(\lambda, \tau) of ' ...
                    'amplitudes over a log-spaced lifetime grid ' ...
                    'with Tikhonov regularization. Useful when the ' ...
                    'kinetics are too complex for a discrete ' ...
                    'few-exponentials fit.']), ...
                'ButtonPushedFcn', @(~,~) app.openLdaWindow());
            uibutton(tbA, 'Text', 'MCR-ALS...', ...
                'BackgroundColor', [0.85 0.50 0.20], 'FontColor', 'w', ...
                'Tooltip', ['Multivariate Curve Resolution by ' ...
                    'Alternating Least Squares: decompose the data ' ...
                    'as D = C * S^T with constraints ' ...
                    '(non-negativity, unimodality). Most model-free ' ...
                    'soft-modeling approach.'], ...
                'ButtonPushedFcn', @(~,~) app.openMcrWindow());
            uibutton(tbA, 'Text', 'Coherence...', ...
                'BackgroundColor', [0.30 0.65 0.75], 'FontColor', 'w', ...
                'Tooltip', ['Vibrational coherence analysis: FFT of ' ...
                    'kinetic residual at each wavelength gives a 2D ' ...
                    'oscillation map (wavelength x wavenumber). ' ...
                    'Reveals coherent vibrational modes of excited ' ...
                    'or ground state.'], ...
                'ButtonPushedFcn', @(~,~) app.openCoherenceWindow());
            uibutton(tbA, 'Text', 'Auto Analysis...', ...
                'BackgroundColor', [0.20 0.50 0.45], 'FontColor', 'w', ...
                'Tooltip', ['Automatic analysis: proposes the corrections ' ...
                    '(background, mask, chirp, crop), estimates the number ' ...
                    'of components and their time constants by global ' ...
                    'fits, and looks for oscillations in the residual. ' ...
                    'The result is a report with the plots it rests on; ' ...
                    'the data are changed only when you apply the ' ...
                    'corrections.'], ...
                'ButtonPushedFcn', @(~,~) app.openAutoWindow());
            uilabel(tbA, 'Text', '');   % spacer (right edge)
            
            % ---- Toolbar row 3: display t range + time unit ----
            tb2 = uigridlayout(gl, [1, 13]);
            tb2.Layout.Row = 3;  tb2.Layout.Column = [1 2];
            tb2.ColumnWidth = {'fit', 80, 'fit', 80, 100, 160, ...
                               105, 105, ...
                               'fit', 70, 'fit', 165, '1x'};
            tb2.ColumnSpacing = 5;  tb2.Padding = [0 0 0 0];
            
            app.lblViewTMin = uilabel(tb2, 'Text', ...
                sprintf('Display t_min (%s):', app.tUnitTxt()), ...
                'HorizontalAlignment', 'right');
            app.edViewTMin = uieditfield(tb2, 'numeric', 'Value', 0, ...
                'ValueChangedFcn', @(src,~) app.onMainViewChange('min', src.Value));
            app.lblViewTMax = uilabel(tb2, 'Text', ...
                sprintf('t_max (%s):', app.tUnitTxt()), ...
                'HorizontalAlignment', 'right');
            app.edViewTMax = uieditfield(tb2, 'numeric', 'Value', 0, ...
                'ValueChangedFcn', @(src,~) app.onMainViewChange('max', src.Value));
            uibutton(tb2, 'Text', 'Reset display', ...
                'ButtonPushedFcn', @(~,~) app.resetMainView());
            uibutton(tb2, 'Text', 'Export Corrected Data...', ...
                'BackgroundColor', [0.30 0.55 0.80], 'FontColor', 'w', ...
                'ButtonPushedFcn', @(~,~) app.exportCorrected());
            uibutton(tb2, 'Text', 'Save Session...', ...
                'Tooltip', ['Save the entire working state (data, ' ...
                    'corrections, fit settings, GA results, view ' ...
                    'options) to a .mat file so you can resume ' ...
                    'later exactly where you left off.'], ...
                'ButtonPushedFcn', @(~,~) app.saveSession());
            uibutton(tb2, 'Text', 'Load Session...', ...
                'Tooltip', 'Restore a session saved with "Save Session".', ...
                'ButtonPushedFcn', @(~,~) app.loadSessionFile());
            
            uilabel(tb2, 'Text', 'Time unit:', ...
                'HorizontalAlignment', 'right', ...
                'FontColor', [0.2 0.2 0.2]);
            app.ddTimeUnit = uidropdown(tb2, ...
                'Items', {'fs', 'ps', 'ns', [char(956) 's']}, ...
                'ItemsData', {'fs', 'ps', 'ns', 'us'}, ...
                'Value', app.timeUnit, ...
                'Tooltip', ['Switch time-axis labels between ' ...
                    'femtoseconds, picoseconds, nanoseconds, and ' ...
                    'microseconds. Numeric values in the loaded data ' ...
                    'are NOT converted; only display labels and ' ...
                    'export headers change.'], ...
                'ValueChangedFcn', @(src,~) app.onTimeUnitChange(src.Value));
            uilabel(tb2, 'Text', 'Spectral axis:', ...
                'HorizontalAlignment', 'right', ...
                'FontColor', [0.2 0.2 0.2]);
            app.ddSpecUnit = uidropdown(tb2, ...
                'Items', {'Wavelength (nm)', 'Wavenumber (cm-1)'}, ...
                'ItemsData', {'nm', 'cm-1'}, ...
                'Value', app.specUnit, ...
                'Tooltip', ['Label the spectral axis as wavelength (nm; ' ...
                    'electronic TA / TRPL) or wavenumber (cm-1; ' ...
                    'TR-Raman / TR-IR). Numeric values in the loaded ' ...
                    'data are NOT converted; only labels and export ' ...
                    'headers change.'], ...
                'ValueChangedFcn', @(src,~) app.onSpecUnitChange(src.Value));
            uilabel(tb2, 'Text', 'full range shown when tmin >= tmax', ...
                'FontColor', [0.45 0.45 0.45], 'FontAngle', 'italic');
            
            % ---- 2D map panel ----
            app.pnl2D = uipanel(gl, ...
                'Title', app.uiStr('2D Map  -  click to select (\lambda, t)'));
            app.pnl2D.Layout.Row = [4 5];
            app.pnl2D.Layout.Column = 1;
            pnl2DGrid = uigridlayout(app.pnl2D, [2, 1]);
            pnl2DGrid.RowHeight = {32, '1x'};
            pnl2DGrid.Padding   = [5 5 5 5];
            pnl2DGrid.RowSpacing = 3;
            
            % Color / Z-range control strip
            ctrl2D = uigridlayout(pnl2DGrid, [1, 8]);
            ctrl2D.Layout.Row = 1;
            ctrl2D.ColumnWidth = {'fit', 110, 'fit', 90, 'fit', 90, 90, '1x'};
            ctrl2D.Padding     = [0 0 0 0];
            ctrl2D.ColumnSpacing = 5;
            
            uilabel(ctrl2D, 'Text', 'Colormap:', ...
                'HorizontalAlignment', 'right');
            app.ddColormap = uidropdown(ctrl2D, ...
                'Items', {'turbo', 'parula', 'jet', 'hot', 'cool', ...
                          'gray', 'RdBu', 'BWR', 'seismic'}, ...
                'Value', app.mapColormap, ...
                'ValueChangedFcn', @(src,~) app.onColormapChange(src.Value));
            uilabel(ctrl2D, 'Text', 'Z min:', ...
                'HorizontalAlignment', 'right');
            app.edZMin = uieditfield(ctrl2D, 'numeric', 'Value', 0, ...
                'ValueDisplayFormat', '%.3g', ...
                'ValueChangedFcn', @(src,~) app.onZChange('min', src.Value));
            uilabel(ctrl2D, 'Text', 'Z max:', ...
                'HorizontalAlignment', 'right');
            app.edZMax = uieditfield(ctrl2D, 'numeric', 'Value', 0, ...
                'ValueDisplayFormat', '%.3g', ...
                'ValueChangedFcn', @(src,~) app.onZChange('max', src.Value));
            uibutton(ctrl2D, 'Text', 'Reset Z', ...
                'ButtonPushedFcn', @(~,~) app.resetZRange());
            uilabel(ctrl2D, 'Text', '');  % spacer
            
            % Axes container
            app.axCont2D = uigridlayout(pnl2DGrid, [1 1]);
            app.axCont2D.Layout.Row = 2;
            app.axCont2D.Padding = [0 0 0 0];
            
            % ---- Spectrum panel ----
            app.pnlSpec = uipanel(gl, 'Title', 'Spectrum  (at selected delay)');
            app.pnlSpec.Layout.Row = 4;  app.pnlSpec.Layout.Column = 2;
            glS = uigridlayout(app.pnlSpec, [2, 1]);
            glS.RowHeight = {'1x', 38};  glS.Padding = [5 5 5 5];
            app.axContSpec = uigridlayout(glS, [1 1]);
            app.axContSpec.Padding = [0 0 0 0];
            ctrlS = uigridlayout(glS, [1, 5]);
            ctrlS.ColumnWidth = {'fit', '1x', 'fit', 'fit', 'fit'};
            ctrlS.Padding = [0 0 0 0];  ctrlS.ColumnSpacing = 5;
            app.lblSelT = uilabel(ctrlS, 'Text', ...
                sprintf('Delay (%s):', app.tUnitTxt()));
            app.edSelT = uieditfield(ctrlS, 'numeric', 'Value', 0, ...
                'ValueChangedFcn', @(src,~) app.setSelT(src.Value));
            uibutton(ctrlS, 'Text', 'Pin', ...
                'ButtonPushedFcn', @(~,~) app.pinSpec());
            uibutton(ctrlS, 'Text', 'Clear pins', ...
                'ButtonPushedFcn', @(~,~) app.clearSpecOverlays());
            uilabel(ctrlS, 'Text', '');
            
            % ---- Kinetics panel ----
            app.pnlKin = uipanel(gl, 'Title', 'Kinetics  (at selected wavelength)');
            app.pnlKin.Layout.Row = 5;  app.pnlKin.Layout.Column = 2;
            glK = uigridlayout(app.pnlKin, [2, 1]);
            glK.RowHeight = {'1x', 38};  glK.Padding = [5 5 5 5];
            app.axContKin = uigridlayout(glK, [1 1]);
            app.axContKin.Padding = [0 0 0 0];
            ctrlK = uigridlayout(glK, [1, 6]);
            ctrlK.ColumnWidth = {'fit', '1x', 'fit', 'fit', 'fit', 'fit'};
            ctrlK.Padding = [0 0 0 0];  ctrlK.ColumnSpacing = 5;
            app.lblSelWL = uilabel(ctrlK, 'Text', ...
                sprintf('%s (%s):', app.xSymTxt(), app.xUnitTxt()));
            app.edSelWL = uieditfield(ctrlK, 'numeric', 'Value', 0, ...
                'ValueChangedFcn', @(src,~) app.setSelWL(src.Value));
            uibutton(ctrlK, 'Text', 'Pin', ...
                'ButtonPushedFcn', @(~,~) app.pinKin());
            uibutton(ctrlK, 'Text', 'Clear pins', ...
                'ButtonPushedFcn', @(~,~) app.clearKinOverlays());
            uibutton(ctrlK, 'Text', 'Fit trace...', ...
                'BackgroundColor', [0.55 0.30 0.70], 'FontColor', 'w', ...
                'Tooltip', ['Fit a sum of exponentials (with optional ' ...
                    'IRF convolution) to the kinetic trace at the ' ...
                    'currently selected wavelength.'], ...
                'ButtonPushedFcn', @(~,~) app.openKineticFitWindow());
            uilabel(ctrlK, 'Text', '');
            
            app.build2DAxes();
            app.buildSpecAxes();
            app.buildKinAxes();
        end
        
        % ----- Dynamic axes builders -----
        function build2DAxes(app)
            delete(app.axCont2D.Children);
            app.ax2D   = gobjects(0);
            app.hVLine = gobjects(0);
            app.hHLine = gobjects(0);
            app.cbarAx = [];
            
            if strcmp(app.delayScaleMode, 'split')
                app.axCont2D.RowHeight   = {'2x', '1x'};
                app.axCont2D.ColumnWidth = {'1x', 85};
                axTop = uiaxes(app.axCont2D);
                axTop.Layout.Row = 1;  axTop.Layout.Column = 1;
                axBot = uiaxes(app.axCont2D);
                axBot.Layout.Row = 2;  axBot.Layout.Column = 1;
                % Dedicated hidden axes for colorbar, spanning both rows
                cbAx = uiaxes(app.axCont2D);
                cbAx.Layout.Row = [1 2];  cbAx.Layout.Column = 2;
                cbAx.Visible = 'off';
                app.cbarAx = cbAx;
                app.ax2D = [axTop, axBot];
            else
                app.axCont2D.RowHeight   = {'1x'};
                app.axCont2D.ColumnWidth = {'1x'};
                app.ax2D = uiaxes(app.axCont2D);
            end
            for k = 1:numel(app.ax2D)
                ax = app.ax2D(k);
                ax.Color = [0.93 0.93 0.93];
                xlabel(ax, app.xAxLabel());
                ylabel(ax, sprintf('Delay time (%s)', app.tUnitAx()));
            end
        end
        
        function buildSpecAxes(app)
            delete(app.axContSpec.Children);
            app.axSpec = uiaxes(app.axContSpec);
            xlabel(app.axSpec, app.xAxLabel());
            ylabel(app.axSpec, '\DeltaA or Intensity');
        end
        
        function buildKinAxes(app)
            delete(app.axContKin.Children);
            app.axKin = gobjects(0);
            if strcmp(app.delayScaleMode, 'split')
                app.axContKin.RowHeight   = {'1x'};
                app.axContKin.ColumnWidth = {'1x', '2x'};
                axL = uiaxes(app.axContKin);  axL.Layout.Column = 1;
                axR = uiaxes(app.axContKin);  axR.Layout.Column = 2;
                app.axKin = [axL, axR];
            else
                app.axContKin.RowHeight   = {'1x'};
                app.axContKin.ColumnWidth = {'1x'};
                app.axKin = uiaxes(app.axContKin);
            end
            for k = 1:numel(app.axKin)
                xlabel(app.axKin(k), sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(app.axKin(k), '\DeltaA or Intensity');
            end
            
            % --- Pre-warm the uifigure web renderer (CEF) ---
            % The first uifigure created in a MATLAB session takes
            % noticeably longer to initialize than subsequent ones --
            % the Chromium-based renderer that drives uifigure has to
            % spin up its process. If a data load fires off a flood of
            % graphics commands before that startup is complete, the
            % SceneTree can race and the 2D map ends up frozen or
            % partially rendered. Force a redraw and brief pause here
            % so the renderer is fully ready before the user loads
            % anything.
            drawnow;
            pause(0.1);
            drawnow;
        end
        
        % =================================================================
        % Data I/O
        % =================================================================
        function loadData(app)
            [files, path] = uigetfile( ...
                {'*.csv;*.dat;*.txt', 'Data (*.csv, *.dat, *.txt)'; ...
                 '*.*', 'All Files'}, ...
                'Select TA data file(s)', app.startDir(), ...
                'MultiSelect', 'on');
            if isequal(files, 0), return; end
            app.rememberDir(path);
            if ~iscell(files), files = {files}; end
            
            if numel(files) == 1
                try
                    [wl, t, A] = app.parseDataFile(fullfile(path, files{1}));
                catch ME
                    uialert(app.fig, sprintf('Failed to read file:\n%s', ...
                        ME.message), 'Load Error');
                    return;
                end
                [wl2, t2, A2, cancel] = app.openSingleFileCheck( ...
                    wl, t, A, files{1});
                if cancel, return; end
                app.setLoadedData(wl2, t2, A2, ...
                    sprintf('Loaded: %s', files{1}));
            else
                app.openLoadAverageWindow(path, files);
            end
        end
        
        function swapLikely = detectAxisSwap(~, wl, t)
            % Heuristic check: does the parsed (wl, t) pair have the two
            % axes swapped? Looks for two independent red flags and
            % returns true if any one is strong enough.
            %
            %   "Wavelength" axis looks suspicious if it contains
            %   non-positive values (wl must be positive) or spans more
            %   than ~1000x (wl is typically 200-2500 nm, at most ~7x).
            %   "Delay" axis looks suspicious if every entry falls in a
            %   typical nm band (150-2500) AND its span is > 50, which is
            %   atypical for delay axes (delays are usually broader).
            wl = wl(:);  t = t(:);
            wlLooksLikeDelay = any(wl <= 0) || ...
                (min(wl) > 0 && max(wl)/min(wl) > 1000);
            tLooksLikeWL = all(t > 150) && all(t < 2500) && ...
                (max(t) - min(t) > 50);
            swapLikely = wlLooksLikeDelay || tLooksLikeWL;
        end
        
        function [wlOut, tOut, Aout, canceled] = openSingleFileCheck( ...
                app, wl, t, A, fname)
            % Modal preview/confirm dialog after reading a single file.
            % Shows the parsed 2D map, reports heuristic axis-swap
            % diagnosis, and lets the user flip the orientation before
            % committing to load.
            canceled = false;
            wlOut = wl;  tOut = t;  Aout = A;
            swapGuess = app.detectAxisSwap(wl, t);
            
            f = uifigure('Name', 'Check axes orientation', ...
                'Position', app.fitFigurePosition(820, 640), 'WindowStyle', 'modal');
            app.applyLightTheme(f);
            gl = uigridlayout(f, [4, 1]);
            gl.RowHeight   = {'fit', '1x', 'fit', 44};
            gl.Padding     = [12 12 12 12];
            gl.RowSpacing  = 8;
            
            info = uilabel(gl, 'Text', '', 'WordWrap', 'on', ...
                'FontColor', [0.2 0.2 0.2]);
            
            mp = uipanel(gl, 'BorderType', 'none');
            mgl = uigridlayout(mp, [1, 1]);  mgl.Padding = [0 0 0 0];
            ax = uiaxes(mgl);
            
            tog = uigridlayout(gl, [1, 3]);
            tog.ColumnWidth = {'fit', 180, '1x'};
            tog.Padding = [0 0 0 0];  tog.ColumnSpacing = 10;
            uilabel(tog, 'Text', ['If wavelength ends up on the y-axis ' ...
                'or delay on the x-axis, use:'], ...
                'FontAngle', 'italic', 'FontColor', [0.35 0.35 0.35]);
            cbT = uicheckbox(tog, 'Text', 'Transpose axes', ...
                'Value', swapGuess, ...
                'Tooltip', ['Swap the parsed wavelength and delay ' ...
                    'axes before loading. Useful when the source file ' ...
                    'was saved with rows = delay and columns = ' ...
                    'wavelength instead of the expected layout.']);
            uilabel(tog, 'Text', fname, ...
                'FontAngle', 'italic', 'FontColor', [0.4 0.4 0.4]);
            
            btnRow = uigridlayout(gl, [1, 3]);
            btnRow.ColumnWidth = {'1x', 140, 140};
            btnRow.Padding = [0 0 0 0];  btnRow.ColumnSpacing = 8;
            uilabel(btnRow, 'Text', '');   % spacer
            btnLoad = uibutton(btnRow, 'Text', 'Load', ...
                'BackgroundColor', [0.30 0.55 0.80], 'FontColor', 'w');
            btnCancel = uibutton(btnRow, 'Text', 'Cancel');
            
            drawPreview();
            cbT.ValueChangedFcn = @(~,~) drawPreview();
            btnLoad.ButtonPushedFcn = @(~,~) doLoad();
            btnCancel.ButtonPushedFcn = @(~,~) doCancel();
            f.CloseRequestFcn = @(~,~) doCancel();
            uiwait(f);
            
            % --- nested ---
            function [wlCur, tCur, Acur] = currentOrientation()
                if cbT.Value
                    wlCur = t(:);  tCur = wl(:).';  Acur = A.';
                    % Ensure monotonicity after swapping axes.
                    [wlCur, iw] = sort(wlCur);  Acur = Acur(iw, :);
                    [tCur,  it] = sort(tCur);   Acur = Acur(:,  it);
                else
                    wlCur = wl(:);  tCur = t(:).';  Acur = A;
                end
            end
            
            function drawPreview()
                [wlCur, tCur, Acur] = currentOrientation();
                msg = sprintf([ ...
                    'File: %s\n\n' ...
                    'With current orientation:\n' ...
                    '  x-axis = %s: %d points, %.3g ... %.3g\n' ...
                    '  y-axis = Delay time: %d points, %.3g ... %.3g'], ...
                    fname, app.pick(strcmp(app.specUnit, 'cm-1'), ...
                        'Wavenumber', 'Wavelength'), ...
                    numel(wlCur), min(wlCur), max(wlCur), ...
                    numel(tCur), min(tCur), max(tCur));
                if swapGuess && ~cbT.Value
                    msg = sprintf(['%s\n\n' ...
                        'WARNING: Heuristic suggests these axes may be ' ...
                        'swapped. Tick "Transpose axes" if the map ' ...
                        'below looks wrong.'], msg);
                elseif swapGuess && cbT.Value
                    msg = sprintf(['%s\n\n' ...
                        'Transpose applied (heuristic had suggested a ' ...
                        'swap). Verify the map below looks right.'], msg);
                end
                info.Text = msg;
                
                cla(ax);
                [W, T] = meshgrid(wlCur, tCur);
                cl = max(abs(Acur(:)), [], 'omitnan');
                if isempty(cl) || cl == 0 || isnan(cl), cl = 1; end
                surface(ax, W, T, zeros(size(Acur.')), Acur.', ...
                    'EdgeColor', 'none', 'FaceColor', 'flat');
                view(ax, 2);
                ax.YDir = 'normal';
                ax.CLim = [-cl cl];
                colormap(ax, app.getColormap());
                cb = colorbar(ax);  cb.Label.String = '\DeltaA or Intensity';
                xlabel(ax, app.xAxLabel());
                ylabel(ax, sprintf('Delay time (%s)', app.tUnitAx()));
                xlim(ax, [min(wlCur), max(wlCur)]);
                if min(tCur) < max(tCur)
                    ylim(ax, [min(tCur), max(tCur)]);
                end
                grid(ax, 'on');
                title(ax, 'Preview of loaded data');
            end
            
            function doLoad()
                [wlOut, tOut, Aout] = currentOrientation();
                uiresume(f);  delete(f);
            end
            
            function doCancel()
                canceled = true;
                uiresume(f);
                if isvalid(f), delete(f); end
            end
        end
        
        function [wl, t, A] = parseDataFile(~, fullpath)
            M = readmatrix(fullpath);
            if size(M,1) < 3 || size(M,2) < 3
                error('TRSpecAnalyzer:parseDataFile:tooSmall', ...
                    'Data matrix is too small (must be at least 3x3).');
            end
            t  = M(1,     2:end);   t  = t(:).';
            wl = M(2:end, 1    );   wl = wl(:);
            A  = M(2:end, 2:end);
            % Rows without a wavelength and columns without a delay are
            % not data: text lines that an export program appended to
            % the file (sample, date, comments) arrive as rows of NaN.
            okW = isfinite(wl);  okT = isfinite(t);
            if nnz(okW) < 2 || nnz(okT) < 2
                error('TRSpecAnalyzer:parseDataFile:noAxes', ...
                    ['The first column (wavelengths) or the first row ' ...
                     '(delays) of the file holds no numbers.']);
            end
            wl = wl(okW);  t = t(okT);  A = A(okW, okT);
            [wl, iw] = sort(wl);  A = A(iw, :);
            [t,  it] = sort(t );  A = A(:,  it);
        end
        
        function setLoadedData(app, wl, t, A, sourceDesc)
            % Store as both original and working data
            app.originalWavelength = wl;
            app.originalDelay      = t;
            app.originalDeltaA     = A;
            app.dataSourceDesc     = sourceDesc;
            
            app.wavelength  = wl;
            app.delay       = t;
            app.deltaA_raw  = A;
            app.deltaA      = A;
            
            % Reset all state
            app.bgApplied = false;  app.bgSpectrum = [];  app.bgN = 0;
            app.chirpApplied = false;
            app.chirpPts = zeros(0,2);  app.chirpParams = [];  app.chirpFitRMS = NaN;
            app.chirpViewTMin = t(1);
            app.chirpViewTMax = min(app.defaultChirpTMax(), t(end));
            app.ridgeWL = [];  app.ridgeT = [];
            app.tZeroShift = 0;
            app.maskedRegions = cell(0, 3);
            app.specOverlays = zeros(0,1);
            app.kinOverlays  = zeros(0,1);
            
            % Default selection
            app.selWL = wl(round(end/2));
            idx0 = find(t >= 0, 1);
            if isempty(idx0), idx0 = round(numel(t)/2); end
            app.selT = t(idx0);
            
            app.edSelWL.Limits = [min(wl), max(wl)];
            app.edSelT.Limits  = [min(t),  max(t) ];
            app.edSelWL.Value  = app.selWL;
            app.edSelT.Value   = app.selT;
            
            % Reset main display range to full
            app.mainViewTMin = t(1);
            app.mainViewTMax = t(end);
            if ~isempty(app.edViewTMin) && isvalid(app.edViewTMin)
                app.edViewTMin.Value = app.mainViewTMin;
                app.edViewTMax.Value = app.mainViewTMax;
            end
            
            % Reset Z range to auto based on new data
            app.mapZMin = [];  app.mapZMax = [];
            if ~isempty(app.edZMin) && isvalid(app.edZMin)
                cl = max(abs(A(:)), [], 'omitnan');
                if isempty(cl) || cl == 0 || isnan(cl), cl = 1; end
                app.edZMin.Value = -cl;
                app.edZMax.Value =  cl;
            end
            
            app.updateStatusLabel();
            app.dataChanged();
            app.updateAll();
            % Yield to the renderer once the first heatmap + spectrum
            % + kinetics have been built. Without this drawnow the
            % whole stack of graphics commands queued by updateAll()
            % can pile up before CEF gets a chance to flush, which is
            % the proximate cause of the "first-load freeze".
            drawnow;
        end
        
        function updateStatusLabel(app)
            if isempty(app.wavelength)
                app.lblStatus.Text = 'No data loaded';
                return;
            end
            base = app.dataSourceDesc;
            if isempty(base), base = 'Data'; end
            if app.isCropped()
                app.lblStatus.Text = sprintf( ...
                    ['%s   (%d ' app.uiStr('\lambda') ...
                     ' x %d t)   [cropped]'], ...
                    base, numel(app.wavelength), numel(app.delay));
            else
                app.lblStatus.Text = sprintf( ...
                    ['%s   (%d ' app.uiStr('\lambda') ' x %d t)'], ...
                    base, numel(app.wavelength), numel(app.delay));
            end
        end
        
        % =================================================================
        % Correction pipeline  (raw -> BG -> chirp)
        % =================================================================
        function recompute(app)
            if isempty(app.deltaA_raw), return; end
            data = app.deltaA_raw;
            if app.bgApplied && ~isempty(app.bgSpectrum)
                data = data - app.bgSpectrum;
            end
            if app.chirpApplied && ~isempty(app.chirpParams)
                t  = app.delay(:).';
                wl = app.wavelength(:);
                t0 = app.chirpModel(app.chirpParams, wl);
                newData = NaN(size(data));
                for i = 1:numel(wl)
                    newData(i,:) = interp1(t, data(i,:), t + t0(i), ...
                        'linear', NaN);
                end
                data = newData;
            end
            % Wavelength region masks (applied last so they survive
            % chirp interpolation).
            if ~isempty(app.maskedRegions)
                for k = 1:size(app.maskedRegions, 1)
                    wl1  = app.maskedRegions{k, 1};
                    wl2  = app.maskedRegions{k, 2};
                    mode = app.maskedRegions{k, 3};
                    mRow = app.wavelength >= wl1 & app.wavelength <= wl2;
                    if ~any(mRow), continue; end
                    if strcmp(mode, 'nan')
                        data(mRow, :) = NaN;
                    else  % 'zero'
                        data(mRow, :) = 0;
                    end
                end
            end
            % Single choke point for every correction (background, chirp,
            % masks): if the rebuilt working data differs from what the
            % analyses were run on, their windows and results are stale.
            changed = ~isempty(app.deltaA) && ~isequaln(data, app.deltaA);
            app.deltaA = data;
            if changed && ~app.holdDataChange
                app.dataChanged();
            end
        end
        
        function base = getChirpBaseData(app)
            base = app.deltaA_raw;
            if app.bgApplied && ~isempty(app.bgSpectrum)
                base = base - app.bgSpectrum;
            end
        end
        
        function t0 = chirpModel(~, p, w)
            % t0(w) = a * sqrt( (b*w^2 - 1) / (c*w^2 - 1) ) + d
            % p = [a, b, c, d].  w may be vector.
            w  = w(:);
            num = p(2).*w.^2 - 1;
            den = p(3).*w.^2 - 1;
            ratio = num ./ den;
            ratio(ratio < 0 | ~isfinite(ratio)) = NaN;
            t0 = p(1) .* sqrt(ratio) + p(4);
        end
        
        function [p, rms] = fitChirpParams(app, pts)
            % Fit the Sellmeier-like t0(w) formula to (w, t) points.
            % Multi-start fminsearch; wavelengths are internally scaled so
            % that fitted b, c are O(1) and then unscaled at the end.
            w = pts(:,1);  tObs = pts(:,2);
            if numel(w) < 4
                error('At least 4 points are required.');
            end
            w0 = median(w);
            ws = w / w0;                      % scaled wavelengths
            
            model = @(p, x) app.chirpModel(p, x);
            obj   = @(p) sum( (tObs - model(p, ws)).^2, 'omitnan' ) ...
                        + 1e6 * any(~isfinite(model(p, ws)));
            
            tRange = max(tObs) - min(tObs);
            tMid   = (max(tObs) + min(tObs)) / 2;
            if tRange == 0, tRange = 1; end
            
            starts = {
                [ tRange,      1.2, 1.5, tMid - tRange]
                [-tRange,      1.2, 1.5, tMid + tRange]
                [ tRange/2,    2.0, 3.0, tMid]
                [ tRange,      1.1, 2.0, tMid - tRange/2]
                [ 2*tRange,    0.8, 1.2, tMid - tRange]
                [ tRange,      5.0, 10.0, tMid - tRange]
            };
            
            opts = optimset('Display', 'off', ...
                'TolX', 1e-9, 'TolFun', 1e-11, ...
                'MaxIter', 20000, 'MaxFunEvals', 20000);
            
            bestP = [];  bestR = inf;
            for k = 1:numel(starts)
                try
                    [pk, rk] = fminsearch(obj, starts{k}, opts);
                    if all(isfinite(pk)) && rk < bestR
                        bestR = rk;  bestP = pk;
                    end
                catch
                end
            end
            if isempty(bestP)
                error('Chirp fit did not converge. Try adjusting points.');
            end
            % Unscale b, c so that the model works in original wavelength units
            bestP(2) = bestP(2) / w0^2;
            bestP(3) = bestP(3) / w0^2;
            p   = bestP;
            rms = sqrt(bestR / numel(w));
        end
        
        % =================================================================
        % Global analysis  (exponential decays convolved with Gaussian IRF)
        % =================================================================
        function c = expIRFConv(~, t, tau, t0, fwhm)
            % Exp decay convolved with a Gaussian IRF (FWHM).
            % Piecewise evaluation to avoid numerical overflow:
            %   far past t0 (dt > threshold)  -> asymptotic exp form
            %   near t0 (|dt| <= threshold)   -> full formula via erfcx
            %   far before t0 (dt < -threshold) -> ~0  (left as initialized)
            sigma = fwhm / (2 * sqrt(2 * log(2)));
            dt = t(:) - t0;
            c  = zeros(size(dt));
            
            if ~isfinite(sigma) || sigma <= 0
                % Degenerate IRF: fall back to ideal step*exp behavior.
                if isinf(tau)
                    c = double(dt >= 0);
                else
                    pos = dt >= 0;
                    c(pos) = exp(-dt(pos)/tau);
                end
                c(~isfinite(c)) = 0;
                return;
            end
            
            if isinf(tau)
                c = 0.5 * erfc(-dt / (sigma * sqrt(2)));
                c(~isfinite(c)) = 0;
                return;
            end
            
            if ~isfinite(tau) || tau <= 0
                % Invalid tau -> return zeros; the fit penalty elsewhere
                % will push the optimizer away from this region.
                return;
            end
            
            % Threshold chosen so that erfc(B) ~ 2 for dt > thresh, regardless
            % of how small tau is (sigma^2/tau term handles tau << sigma case).
            thresh = 5*sigma + sigma^2/tau;
            
            near   = abs(dt) <= thresh;
            farPos = dt > thresh;
            
            % Asymptotic past the IRF: pure decay with sigma^2/(2*tau^2) offset
            if any(farPos)
                dtf = dt(farPos);
                c(farPos) = exp(sigma^2/(2*tau^2) - dtf/tau);
            end
            
            % Full convolution near t0, evaluated with erfcx for stability
            if any(near)
                dtn = dt(near);
                B   = (sigma/tau - dtn/sigma) / sqrt(2);
                c(near) = 0.5 * exp(-dtn.^2 / (2*sigma^2)) .* erfcx(B);
            end
            
            c(~isfinite(c)) = 0;
        end
        
        function c = stretchedIRFConv(app, t, tau, beta, t0, fwhm)
            %STRETCHEDIRFCONV  Stretched-exp f(t) = exp(-((t-t0)/tau)^beta)
            %   (zero before t0) convolved with a Gaussian IRF by
            %   Gauss-Legendre quadrature at each delay. Without a valid
            %   IRF width (FWHM not finite or <= 0) the bare function is
            %   returned.
            %
            %   For beta == 1 this reduces to the standard exp; the
            %   caller is free to use expIRFConv instead in that case
            %   (which has a closed-form IRF treatment via erfcx).
            t  = t(:);
            dt = t - t0;
            c  = zeros(size(dt));
            
            if ~isfinite(tau) || tau <= 0 || ~isfinite(beta) || beta <= 0
                return;
            end
            
            sigma = fwhm / (2 * sqrt(2 * log(2)));
            if ~isfinite(sigma) || sigma <= 0
                % Bare stretched exponential (no IRF).
                c = exp(-(max(dt, 0) ./ tau) .^ beta) .* double(dt >= 0);
                return;
            end
            % c(x) = int_0^inf f(u) g(x - u) du,  x = t - t0,
            % with g the unit-area Gaussian IRF. g limits the
            % integrand to u in [x - 6 sigma, x + 6 sigma] and
            % f = 0 for u < 0, so every requested delay is
            % integrated over [max(0, x - 6 sigma), x + 6 sigma]
            % with a fixed Gauss-Legendre rule (two panels, so
            % that a decay much faster than the IRF is still
            % resolved next to u = 0).
            % Starting the interval AT u = 0 keeps the step of f
            % out of the integrand: a uniform-grid convolution
            % straddles that step and is only first-order
            % accurate around t0 (errors of a few % of the
            % amplitude). The cost is also independent of the
            % length of the delay axis.
            [gx, gw] = app.gaussLegendre();
            lo  = max(dt - 6 * sigma, 0);
            hi  = dt + 6 * sigma;
            act = hi > 0;
            if any(act)
                xa  = dt(act);
                loA = lo(act);
                hiA = hi(act);
                % Panel split a few tau above the lower limit
                % (never beyond the middle of the interval).
                spl = loA + min(4 * tau, (hiA - loA) / 2);
                acc = zeros(size(xa));
                for pnl = 1:2
                    if pnl == 1
                        a = loA;  b = spl;
                    else
                        a = spl;  b = hiA;
                    end
                    half = (b - a) / 2;
                    mid  = (b + a) / 2;
                    U = max(mid + half * gx.', 0);
                    F = exp(-(U ./ tau) .^ beta);
                    G = exp(-(xa - U).^2 / (2 * sigma^2));
                    acc = acc + half .* ((F .* G) * gw);
                end
                c(act) = acc / (sigma * sqrt(2 * pi));
            end
            c(~isfinite(c)) = 0;
        end
        
        function [x, w] = gaussLegendre(app)
            % Nodes x (in (-1, 1)) and weights w of a 48-point
            % Gauss-Legendre rule (Golub-Welsch), computed once and
            % cached. int_a^b f ~ (b-a)/2 * sum(w .* f(mid + half*x)).
            if isempty(app.quadGL)
                K  = 48;
                k  = 1:(K - 1);
                b  = k ./ sqrt(4 * k.^2 - 1);
                [V, E] = eig(diag(b, 1) + diag(b, -1));
                [xs, ix] = sort(diag(E));
                ws = 2 * (V(1, ix).') .^ 2;
                app.quadGL = struct('x', xs(:), 'w', ws(:));
            end
            x = app.quadGL.x;
            w = app.quadGL.w;
        end

        function ca = caOptions(~, ca)
            %caOptions  Settings of the coherent-artifact terms of a fit,
            %   with defaults.
            %   order : 0 none | 1 Gaussian | 2 + its 1st time derivative
            %           | 3 + its 2nd time derivative
            %   mode  : 'irf'   the terms have the width of the IRF
            %           'free'  width of their own, fitted
            %           'fixed' width of their own, not fitted
            %   fwhm  : that own width (starting or fixed value)
            d = struct('order', 0, 'mode', 'irf', 'fwhm', NaN);
            if ~isstruct(ca) || isempty(ca), ca = d;  return; end
            fn = fieldnames(d);
            for k = 1:numel(fn)
                if ~isfield(ca, fn{k}) || isempty(ca.(fn{k}))
                    ca.(fn{k}) = d.(fn{k});
                end
            end
            ca = ca(1);
            ca.order = max(0, min(3, round(double(ca.order))));
            if ~ischar(ca.mode) || ~any(strcmp(ca.mode, {'irf', 'free', 'fixed'}))
                ca.mode = 'irf';
            end
        end
        
        function B = caBasis(~, t, t0, fwhm, order)
            %caBasis  Time profiles of the coherent-artifact terms: a
            %   Gaussian of the given FWHM centred at t0 and its first and
            %   second derivative with respect to time, each scaled to a
            %   size of order 1 (numel(t) x order).
            %   A signal that exists only while pump and probe overlap
            %   (cross-phase modulation, two-photon absorption, stimulated
            %   Raman of the solvent) is a combination of these three.
            t = t(:);
            sg = fwhm / (2 * sqrt(2 * log(2)));
            x = (t - t0) / sg;
            G = exp(-x.^2 / 2);
            B = [G, -x .* G, (x.^2 - 1) .* G];
            B = B(:, 1:order);
        end
        
        function D = gaDataNoCA(app)
            %gaDataNoCA  Working data minus the coherent-artifact part of
            %   the stored global fit (the data themselves when that fit
            %   has no artifact terms). Used where spectra are solved
            %   again from the data with the kinetic model alone.
            D = app.deltaA;
            ca = app.gaResultCA;
            if isstruct(ca) && isfield(ca, 'order') && ca.order > 0 && ...
                    size(ca.spec, 1) == size(D, 1) && ~isempty(app.gaResultT0)
                B = app.caBasis(app.delay, app.gaResultT0, ca.fwhm, ca.order);
                sp = ca.spec;  sp(~isfinite(sp)) = 0;
                D = D - sp * B.';
            end
        end
        
        function C = buildGaBasis(app, t, tauVec, t0, fwhm, hasInf)
            k = numel(tauVec) + double(hasInf);
            C = zeros(numel(t), k);
            for j = 1:numel(tauVec)
                C(:, j) = app.expIRFConv(t, tauVec(j), t0, fwhm);
            end
            if hasInf
                C(:, end) = app.expIRFConv(t, Inf, t0, fwhm);
            end
        end
        
        function [tauOut, betaOut, t0Out, fwhmOut, A, fitMat, info] = ...
                fitGlobalAnalysis(app, D, t, tauInit, t0Init, fwhmInit, ...
                tauFixed, t0Fixed, fwhmFixed, hasInf, ...
                betaInit, stretchOn, irfMode, tRange, ca, fwMode)
            % Global fit via VARPRO: nonlinear params (tau, beta, t0, FWHM)
            % via fminsearch, amplitudes A solved linearly per wavelength.
            %
            % Stretched-exponential support:
            %   stretchOn(j)=true  -> component j uses
            %                         exp(-((t-t0)/tau_j)^beta_j)  conv IRF
            %                         (Gauss-Legendre quadrature)
            %   stretchOn(j)=false -> component j uses standard exp
            %                         (closed-form IRF via erfcx)
            %   beta_j is globally shared across all wavelengths.
            %   irfMode = 'none' (no IRF at all) | 'numerical' (Gaussian
            %             IRF; any other value is treated like this one)
            %   tRange  = optional [tmin tmax] fit window ([] = all delays)
            %   ca      = optional coherent-artifact terms (see caOptions):
            %             extra columns of the basis at time zero with free
            %             spectra. A (the DADS) holds the kinetic components
            %             only; the artifact spectra are in info.caSpec.
            %   fwMode  = what is done about a fitted IRF FWHM that the
            %             delays do not resolve (see below): 'auto'
            %             (default) | 'plain' | 'limit'
            %
            % A fitted FWHM and the delay step. A response narrower than
            % about half the delay step around time zero looks the same
            % at the measured delays as any narrower one: below that
            % width (info.fwhmLo) the loss is flat, and a search that
            % gets there -- typically early on, while the time constants
            % are still far off -- does not come back.
            %   'plain' : the fit as it is, nothing else.
            %   'limit' : the search is kept above that width (or above
            %             the start value of the FWHM, when that is
            %             smaller).
            %   'auto'  : the plain fit. When it lost the width -- the
            %             FWHM is below that width or within 20 % of it
            %             (irfUnresolved), or its uncertainty exceeds its
            %             value -- the fit is repeated in two steps: time
            %             constants with t0 and FWHM held at the start
            %             values ('plain'), then everything free from
            %             there ('limit'). A start value of the FWHM
            %             that is itself not resolved is replaced by two
            %             delay steps in the repeat. The repeat is
            %             returned when it fits better, or as well with
            %             the first FWHM not resolved.
            % info.irfRetry: 0 no repeat, 1 result of the repeat,
            % 2 repeated, first result kept. info.fwhmAtLimit: the
            % returned FWHM is not resolved (irfUnresolved).
            %
            % Backwards compatibility: if betaInit/stretchOn are not given
            % (or stretchOn is all-false), behavior matches the previous
            % standard-exp-only fit exactly.
            if nargin < 11 || isempty(betaInit)
                betaInit = ones(1, numel(tauInit));
            end
            if nargin < 12 || isempty(stretchOn)
                stretchOn = false(1, numel(tauInit));
            end
            if nargin < 13 || isempty(irfMode)
                irfMode = 'numerical';
            end
            if nargin < 14, tRange = []; end
            if nargin < 15, ca = []; end
            if nargin < 16 || isempty(fwMode), fwMode = 'auto'; end
            ca = app.caOptions(ca);
            if strcmpi(irfMode, 'none'), ca.order = 0; end
            caK    = ca.order;                           % number of terms
            caOwn  = caK > 0 && ~strcmp(ca.mode, 'irf');  % width of its own
            caFree = caK > 0 && strcmp(ca.mode, 'free');  % ... and fitted
            caW    = ca.fwhm;
            if caOwn && ~(caW > 0), caW = fwhmInit; end
            stretchOn = logical(stretchOn(:).');
            betaInit  = betaInit(:).';
            % Row vectors throughout, whatever orientation was passed;
            % a scalar option applies to every component.
            tauInit  = tauInit(:).';
            tauFixed = logical(tauFixed(:).');
            nCmpGA   = numel(tauInit);
            if isscalar(tauFixed),  tauFixed  = repmat(tauFixed,  1, nCmpGA); end
            if isscalar(betaInit),  betaInit  = repmat(betaInit,  1, nCmpGA); end
            if isscalar(stretchOn), stretchOn = repmat(stretchOn, 1, nCmpGA); end
            if numel(tauFixed) ~= nCmpGA || numel(betaInit) ~= nCmpGA || ...
                    numel(stretchOn) ~= nCmpGA
                error(['fitGlobalAnalysis: tauFixed, betaInit and ' ...
                       'stretchOn need one entry per time constant ' ...
                       '(%d).'], nCmpGA);
            end
            
            % Optional fit-time-range restriction: delay points outside
            % [tRange(1), tRange(2)] are excluded from BOTH the linear
            % amplitude solve and the loss / residuals. Default = use
            % the full delay axis.
            if nargin < 14 || isempty(tRange) || ~all(isfinite(tRange))
                tSel = true(1, numel(t));
            else
                tSel = (t(:).' >= min(tRange)) & (t(:).' <= max(tRange));
                if ~any(tSel), tSel = true(1, numel(t)); end
            end
            
            tauCur  = tauInit(:).';
            betaCur = betaInit;
            t0Cur   = t0Init;
            fwhmCur = fwhmInit;
            
            % In 'none' (no-IRF) mode t0 is exactly degenerate with the
            % linear amplitudes for pure exponentials
            % (exp(-(t-t0)/tau) = e^{t0/tau} exp(-t/tau)), and FWHM does
            % not enter the model at all. Force both fixed so the
            % optimizer never wanders along flat directions.
            if strcmpi(irfMode, 'none')
                t0Fixed   = true;
                fwhmFixed = true;
            end
            
            idxFreeTau  = find(~tauFixed);
            idxFreeBeta = find(stretchOn);  % all stretched betas are free
            nFreeTau    = numel(idxFreeTau);
            nFreeBeta   = numel(idxFreeBeta);
            
            % Pack free params (log-transform for positivity)
            x0 = zeros(0, 1);
            if nFreeTau > 0
                x0 = [x0; log(tauInit(idxFreeTau(:))).'];
            end
            if nFreeBeta > 0
                x0 = [x0; log(max(betaInit(idxFreeBeta(:)), 1e-3)).'];
            end
            if ~t0Fixed,   x0 = [x0; t0Init];         end
            if ~fwhmFixed, x0 = [x0; log(fwhmInit)];  end
            if caFree,     x0 = [x0; log(caW)];       end
            
            % Delay mask of the fit (the user's fit range), fixed for
            % the whole optimisation.
            fitMask = tSel;
            if ~any(fitMask)
                error(['No delay points left to fit (check the fit ' ...
                       't range and t0).']);
            end
            % Smallest FWHM the delays resolve, and the lower limit of
            % the search in 'limit' mode (see the header). The limit is
            % in force during the search only, not while the
            % uncertainties are evaluated at the result.
            fwRes = 0;  fwLo = 0;  fwWall = false;
            if ~fwhmFixed
                dzIrf = app.irfDelayStep(t(fitMask), t0Init);
                if isfinite(dzIrf) && dzIrf > 0
                    fwRes  = 0.5 * dzIrf;
                    fwLo   = min(fwRes, fwhmInit) * (1 - 1e-9);
                    fwWall = strcmp(fwMode, 'limit');
                end
            end
            if caK > 0
                % The artifact terms live within about one width of time
                % zero. Without enough delays there they cannot be told
                % apart.
                if caOwn, w0 = caW; else, w0 = fwhmInit; end
                nNear = nnz(fitMask & abs(t(:).' - t0Init) <= 1.5 * w0);
                if nNear < caK + 2
                    error(['Coherent artifact: %d fitted delay(s) lie ' ...
                        'within 1.5 FWHM of time zero, %d are needed for ' ...
                        '%d term(s). Use fewer terms or switch them off.'], ...
                        nNear, caK + 2, caK);
                end
            end
            nLoss = 0;     % residuals entering the last loss evaluation
            
            function [loss, As, fitM] = objective(x)
                idx = 1;
                for j = 1:nFreeTau
                    tauCur(idxFreeTau(j)) = exp(x(idx));
                    idx = idx + 1;
                end
                for j = 1:nFreeBeta
                    betaCur(idxFreeBeta(j)) = exp(x(idx));
                    idx = idx + 1;
                end
                if ~t0Fixed,   t0Cur   = x(idx);        idx = idx + 1; end
                if ~fwhmFixed, fwhmCur = exp(x(idx));   idx = idx + 1; end
                if caFree,     caW     = exp(x(idx));                  end
                
                % Reject unphysical regions early. A fitted artifact
                % width stays within a factor of 5 of the IRF width.
                sigmaCur = fwhmCur / (2*sqrt(2*log(2)));
                if strcmpi(irfMode, 'none'), sigmaCur = 0; end
                tSpan = max(t) - min(t);
                if any(tauCur < sigmaCur/4) || any(tauCur > 10 * tSpan) || ...
                        any(tauCur <= 0) || ...
                        (~strcmpi(irfMode, 'none') && ...
                         (fwhmCur <= 0 || fwhmCur > tSpan)) || ...
                        (fwWall && fwhmCur < fwLo) || ...
                        (caFree && (caW < fwhmCur / 5 || caW > 5 * fwhmCur)) || ...
                        any(betaCur(stretchOn) <= 0) || ...
                        any(betaCur(stretchOn) > 2)
                    As   = zeros(size(D, 1), numel(tauCur) + double(hasInf) + caK);
                    fitM = zeros(size(D));
                    loss = 1e30;
                    return;
                end
                
                % Build basis with mixed exp / stretched columns, then
                % the coherent-artifact terms.
                C = app.buildGaBasisExt(t, tauCur, betaCur, ...
                    stretchOn, t0Cur, fwhmCur, hasInf, irfMode);
                if caK > 0
                    if caOwn, wNow = caW; else, wNow = fwhmCur; end
                    C = [C, app.caBasis(t, t0Cur, wNow, caK)];
                end
                colNorms = sqrt(sum(C.^2, 1, 'omitnan'));
                if any(colNorms < 1e-12) || any(~isfinite(C(:)))
                    As   = zeros(size(D, 1), size(C, 2));
                    fitM = zeros(size(D));
                    loss = 1e30;
                    return;
                end
                
                % Column mask for the loss (fixed, see above).
                cMask = fitMask;
                
                % Solve amplitudes. If D has NaN, fall back to per-wavelength;
                % otherwise use vectorized min-norm LS.
                [M, N] = size(D);
                k = size(C, 2);
                As   = zeros(M, k);
                fitM = zeros(M, N);
                if all(isfinite(D(:))) && all(cMask)
                    try
                        At = lsqminnorm(C, D.');
                    catch
                        At = C \ D.';
                    end
                    if any(~isfinite(At(:))) || max(abs(At(:))) > 1e10
                        loss = 1e30;  return;
                    end
                    As   = At.';
                    fitM = As * C.';
                elseif all(isfinite(D(:)))
                    % Restricted columns but no NaN: one shared solve on
                    % the selected columns, reconstruction on all.
                    try
                        At = lsqminnorm(C(cMask, :), D(:, cMask).');
                    catch
                        At = C(cMask, :) \ D(:, cMask).';
                    end
                    if any(~isfinite(At(:))) || max(abs(At(:))) > 1e10
                        loss = 1e30;  return;
                    end
                    As   = At.';
                    fitM = As * C.';
                else
                    for ii = 1:M
                        row = D(ii, :);
                        mm  = isfinite(row) & cMask;
                        if nnz(mm) <= k, continue; end
                        Ci = C(mm, :);
                        di = row(mm).';
                        try
                            ai = lsqminnorm(Ci, di);
                        catch
                            ai = Ci \ di;
                        end
                        if any(~isfinite(ai)) || max(abs(ai)) > 1e8
                            continue;
                        end
                        As(ii, :)    = ai.';
                        fitM(ii, :)  = (ai.') * C.';
                    end
                end
                
                R = D - fitM;
                mfin = isfinite(R) & repmat(cMask, M, 1);
                nLoss = nnz(mfin);
                if ~any(mfin(:))
                    loss = 1e30;
                else
                    loss = sum(R(mfin).^2);
                end
            end
            
            if isempty(x0)
                % No free nonlinear params -- just compute amplitudes
                [loss, A, fitMat] = objective([]);
                iters = 0;
                initLoss = loss;
                nInit = nLoss;  nOpt = nLoss;
            else
                % Evaluate initial loss before optimization (diagnostic)
                [initLoss, ~, ~] = objective(x0);
                nInit = nLoss;
                opts = optimset('Display', 'off', ...
                    'TolX', 1e-8, 'TolFun', 1e-10, ...
                    'MaxIter', 5000, 'MaxFunEvals', 20000);
                [xOpt, loss, ~, out] = fminsearch(@objective, x0, opts);
                [~, A, fitMat] = objective(xOpt);
                nOpt = nLoss;
                iters = out.iterations;
            end
            tauOut  = tauCur;
            betaOut = betaCur;
            t0Out   = t0Cur;
            fwhmOut = fwhmCur;
            
            fwWall = false;      % the search is over

            % ---- 1-sigma uncertainties (numerical Jacobian) ----
            % Computed on the PHYSICAL parameter scale (tau, beta, t0,
            % fwhm), not the internal log scale, so the reported sigma
            % can be quoted directly as tau +/- sigma.
            tauStd  = NaN(1, numel(tauCur));
            betaStd = NaN(1, numel(betaCur));
            t0Std   = NaN;   fwhmStd = NaN;   caStd = NaN;
            if ~isempty(x0)
                pOpt = zeros(numel(x0), 1);
                idx = 1;
                for j = 1:nFreeTau
                    pOpt(idx) = tauCur(idxFreeTau(j));  idx = idx + 1;
                end
                for j = 1:nFreeBeta
                    pOpt(idx) = betaCur(idxFreeBeta(j)); idx = idx + 1;
                end
                if ~t0Fixed,   pOpt(idx) = t0Cur;   idx = idx + 1; end
                if ~fwhmFixed, pOpt(idx) = fwhmCur; idx = idx + 1; end
                if caFree,     pOpt(idx) = caW;                    end
                
                % Degrees of freedom: residuals in the loss minus ALL
                % fitted parameters -- the linear amplitudes here, the
                % nonlinear ones inside fitParamUncertainty.
                nLin  = size(A, 2) * nnz(any(A ~= 0, 2));
                nUsed = nOpt - nLin;
                sig = app.fitParamUncertainty( ...
                    @(p) gaResidualAtPhys(p), pOpt, loss, nUsed);
                
                idx = 1;
                for j = 1:nFreeTau
                    tauStd(idxFreeTau(j)) = sig(idx);  idx = idx + 1;
                end
                for j = 1:nFreeBeta
                    betaStd(idxFreeBeta(j)) = sig(idx); idx = idx + 1;
                end
                if ~t0Fixed,   t0Std   = sig(idx);  idx = idx + 1; end
                if ~fwhmFixed, fwhmStd = sig(idx);  idx = idx + 1; end
                if caFree,     caStd   = sig(idx);                 end
                
                % The Jacobian probing perturbed the nested state
                % (tauCur etc.); restore the optimum.
                [~, A, fitMat] = objective(xOpt);
            end
            
            % Kinetic components and artifact terms are returned apart.
            kKin   = numel(tauCur) + double(hasInf);
            caSpec = A(:, kKin+1:end);
            A      = A(:, 1:kKin);
            caWout = NaN;
            if caK > 0
                if caOwn, caWout = caW; else, caWout = fwhmCur; end
            end
            
            function r = gaResidualAtPhys(p)
                % Map physical params -> internal log-packed x, run the
                % objective, and return the masked residual vector.
                xx = zeros(numel(p), 1);
                ii2 = 1;
                for jj = 1:nFreeTau
                    xx(ii2) = log(max(p(ii2), 1e-12));  ii2 = ii2 + 1;
                end
                for jj = 1:nFreeBeta
                    xx(ii2) = log(max(p(ii2), 1e-12));  ii2 = ii2 + 1;
                end
                if ~t0Fixed,   xx(ii2) = p(ii2);        ii2 = ii2 + 1; end
                if ~fwhmFixed
                    xx(ii2) = log(max(p(ii2), 1e-12));  ii2 = ii2 + 1;
                end
                if caFree,     xx(ii2) = log(max(p(ii2), 1e-12));      end
                [~, ~, fitM_] = objective(xx);
                R_ = D - fitM_;
                mf = isfinite(D) & repmat(fitMask, size(R_, 1), 1);
                r = R_(mf);
                r(~isfinite(r)) = 0;
            end
            
            info = struct('rss', loss, 'iters', iters, ...
                'rms', sqrt(loss / max(1, nOpt)), ...
                'initialLoss', initLoss, ...
                'initialRMS', sqrt(initLoss / max(1, nInit)), ...
                'nUsed', nOpt, ...
                'tauStd', tauStd, 'betaStd', betaStd, ...
                't0Std', t0Std, 'fwhmStd', fwhmStd, ...
                'caOrder', caK, 'caSpec', caSpec, 'caFwhm', caWout, ...
                'caFwhmStd', caStd, 'caMode', ca.mode, ...
                'fwhmLo', NaN, 'fwhmAtLimit', false, 'irfRetry', 0);
            if ~fwhmFixed && fwRes > 0
                info.fwhmLo      = fwRes;
                info.fwhmAtLimit = app.irfUnresolved(fwhmOut, fwRes);
            end

            % ---- two-step repeat when the IRF width was lost ----
            if ~fwhmFixed && strcmp(fwMode, 'auto') && ~isempty(x0) && ...
                    app.irfWidthLost(fwhmOut, fwhmStd, fwRes)
                % A start value that the delays do not resolve gives
                % the search nothing to go by: two delay steps instead.
                fwS0 = fwhmInit;
                if app.irfUnresolved(fwhmInit, fwRes), fwS0 = 4 * fwRes; end
                try
                    [tauS1, betaS1, ~, ~, ~, ~, infoS1] = ...
                        app.fitGlobalAnalysis(D, t, tauInit, t0Init, ...
                        fwS0, tauFixed, true, true, hasInf, betaInit, ...
                        stretchOn, irfMode, tRange, ca, 'plain');
                    [tauS2, betaS2, t0S2, fwS2, AS2, fitS2, infoS2] = ...
                        app.fitGlobalAnalysis(D, t, tauS1, t0Init, ...
                        fwS0, tauFixed, t0Fixed, fwhmFixed, hasInf, ...
                        betaS1, stretchOn, irfMode, tRange, ca, 'limit');
                    itAll = info.iters + infoS1.iters + infoS2.iters;
                    if infoS2.rss < info.rss * (1 - 1e-6) || ...
                            (info.fwhmAtLimit && ...
                             infoS2.rss <= info.rss * (1 + 1e-6))
                        infoS2.initialLoss = info.initialLoss;
                        infoS2.initialRMS  = info.initialRMS;
                        tauOut = tauS2;  betaOut = betaS2;
                        t0Out  = t0S2;   fwhmOut = fwS2;
                        A = AS2;  fitMat = fitS2;
                        info = infoS2;
                        info.irfRetry = 1;
                    else
                        info.irfRetry = 2;
                    end
                    info.iters = itAll;
                catch
                    % the repeat failed: the first result stands
                end
            end
        end
        
        function dz = irfDelayStep(~, t, t0)
            %irfDelayStep  Delay step around time zero: the median of the
            %   (up to four) intervals between the delays next to t0. NaN
            %   when there are fewer than two different delays.
            dz = NaN;
            d = unique(t(isfinite(t)));
            if numel(d) < 2 || ~isfinite(t0), return; end
            [~, j] = min(abs(d - t0));
            k = max(1, j - 2):min(numel(d), j + 2);
            dz = median(diff(d(k)));
        end

        function tf = irfUnresolved(~, fw, fwRes)
            %irfUnresolved  True when a fitted IRF FWHM is not resolved by
            %   the delays: it is below the smallest width they resolve
            %   (fwRes, half the delay step around time zero; 0 =
            %   unknown) or within 20 % of it. Up to there the loss
            %   hardly changes with the width.
            tf = fwRes > 0 && fw <= 1.2 * fwRes;
        end
        
        function tf = irfWidthLost(app, fw, fwStd, fwRes)
            %irfWidthLost  True when a fitted IRF FWHM is not a usable
            %   result: the delays do not resolve it (irfUnresolved), its
            %   uncertainty is not a finite number, or the uncertainty
            %   exceeds the value.
            tf = app.irfUnresolved(fw, fwRes) || ...
                ~isfinite(fwStd) || fwStd > fw;
        end

        function L = irfLimitNote(app, info)
            %irfLimitNote  Lines for the result text of a fit (1 x n
            %   cell) when the fitted IRF FWHM is not resolved by the
            %   delays, and when the result comes from the two-step
            %   repeat (see fitGlobalAnalysis). Empty when there is
            %   nothing to say. info: struct with fwhmAtLimit, fwhmLo,
            %   irfRetry (missing fields are taken as "no").
            L = {};
            if ~isstruct(info) || ~isscalar(info), return; end
            if isfield(info, 'fwhmAtLimit') && isequal(info.fwhmAtLimit, true)
                lim = NaN;
                if isfield(info, 'fwhmLo') && isscalar(info.fwhmLo)
                    lim = info.fwhmLo;
                end
                L{end+1} = '  FWHM is close to, or below, half the';
                L{end+1} = sprintf('  delay step at time zero (%.3g %s):', ...
                    lim, app.tUnitTxt());
                L{end+1} = '  the delays do not resolve such a';
                L{end+1} = '  width. Fix it at a known value.';
            end
            if isempty(L) && isfield(info, 'irfRetry') && ...
                    isequal(info.irfRetry, 1)
                L{end+1} = '  The direct fit lost the IRF width;';
                L{end+1} = '  this is the result of a fit in two';
                L{end+1} = '  steps (IRF fixed, then free).';
            end
        end

        function sig = fitParamUncertainty(~, resFun, p, rss, nData)
            %fitParamUncertainty  1-sigma parameter uncertainties from a
            %   numerical Jacobian at the optimum.
            %
            %   resFun : @(p) -> residual column vector at parameters p
            %            (physical scale, NOT log-transformed)
            %   p      : optimal parameter vector (1xM or Mx1)
            %   rss    : residual sum of squares at optimum
            %   nData  : number of data points used in the fit
            %
            %   Returns sig (Mx1): sqrt(diag(s^2 * inv(J'J))) with
            %   s^2 = rss / (nData - M); callers pass nData already
            %   reduced by the number of linear amplitudes, so that s^2
            %   uses the full parameter count. Linear-approximation CIs --
            %   valid near the optimum for well-conditioned fits.
            %   Inf marks a parameter the data do not determine: it has
            %   no effect on the residual, or it acts only in combination
            %   with another parameter (e.g. two rate constants that
            %   enter through their sum). NaN = could not be evaluated.
            p = p(:);
            M = numel(p);
            sig = NaN(M, 1);
            if M == 0 || nData <= M, return; end
            try
                r0 = resFun(p);
                n  = numel(r0);
                J  = zeros(n, M);
                for j = 1:M
                    % Relative central-difference step
                    h = max(abs(p(j)), 1e-12) * 1e-4;
                    pp = p;  pp(j) = p(j) + h;
                    pm = p;  pm(j) = p(j) - h;
                    J(:, j) = (resFun(pp) - resFun(pm)) / (2 * h);
                end
                s2  = rss / max(1, (nData - M));
                % Unit-norm columns: the rank test must not depend on
                % the units of the parameters (a tau of 900 next to a
                % t0 of 0.1).
                cn   = sqrt(sum(J.^2, 1));
                dead = ~(cn > 0) | ~isfinite(cn);   % no effect at all
                und  = dead(:);
                live = find(~dead);
                if ~isempty(live)
                    Jn = J(:, live) ./ cn(live);
                    [~, S, V] = svd(Jn, 'econ');
                    sv = diag(S);
                    % Directions the data do not constrain. A pseudo-
                    % inverse would silently drop them and report a
                    % SMALL error for exactly the parameters that are
                    % not determined, so flag those parameters instead.
                    % The threshold reflects the accuracy of the
                    % finite-difference Jacobian.
                    flat = sv < 1e-5 * sv(1);
                    if any(flat)
                        wNull = sum(V(:, flat).^2, 2);
                        und(live(wNull > 1e-6)) = true;
                    end
                    Vg = V(:, ~flat);
                    dg = sum((Vg ./ sv(~flat).').^2, 2);
                    sg = sqrt(s2 * dg) ./ cn(live).';
                    okS = isfinite(sg);
                    sig(live(okS)) = sg(okS);
                end
                sig(und) = Inf;
            catch
                % leave NaN
            end
        end
        
        function C = buildGaBasisExt(app, t, tauVec, betaVec, ...
                stretchOn, t0, fwhm, hasInf, irfMode)
            %buildGaBasisExt  Like buildGaBasis but with optional
            %   stretched columns. Columns are aligned with tauVec order;
            %   the optional infinite-tau column (when hasInf=true)
            %   always uses standard exp. irfMode 'none' gives the bare
            %   kinetics; any other value includes the Gaussian IRF.
            n = numel(tauVec) + double(hasInf);
            C = zeros(numel(t), n);
            if strcmpi(irfMode, 'none')
                % Pure kinetics without any IRF convolution: each column
                % is exp(-((t-t0)/tau)^beta) . H(t-t0) (beta = 1 for
                % standard components). Intended for data cropped to
                % start after the IRF.
                tt = t(:) - t0;
                step = double(tt >= 0);
                for j = 1:numel(tauVec)
                    if stretchOn(j)
                        C(:, j) = step .* ...
                            exp(-(max(tt, 0) ./ tauVec(j)).^betaVec(j));
                    else
                        C(:, j) = step .* exp(-max(tt, 0) ./ tauVec(j));
                    end
                end
                if hasInf
                    C(:, end) = step;
                end
                return;
            end
            for j = 1:numel(tauVec)
                if stretchOn(j)
                    C(:, j) = app.stretchedIRFConv(t, tauVec(j), ...
                        betaVec(j), t0, fwhm);
                else
                    C(:, j) = app.expIRFConv(t, tauVec(j), t0, fwhm);
                end
            end
            if hasInf
                C(:, end) = app.expIRFConv(t, Inf, t0, fwhm);
            end
        end
        
        function [EADS, tauSorted, Bmat] = computeEADSFromDADS( ...
                ~, DADS, tauVec, hasInf)
            %computeEADSFromDADS  Convert parallel DADS to sequential EADS.
            %   Sequential model: species 1 -> 2 -> ... -> N (+ optional
            %   terminal species when hasInf=true). The concentration of
            %   species i under irreversible first-order chain kinetics with
            %   rates k_j = 1/tau_j (and k = 0 for the terminal species) is
            %       c_i(t) = sum_{j<=i} B(i,j) * exp(-k_j t)
            %   with B(i,j) = (prod_{m<i} k_m) / prod_{m<=i, m~=j}(k_m - k_j).
            %   From D = DADS * exp(-K t) = EADS * C(t) and C = B * exp(-K t),
            %   we get DADS = EADS * B, so EADS = DADS / B.
            %   Species are ordered by ascending tau; the terminal (tau=inf)
            %   species sits in the last column.
            nDecay = numel(tauVec);
            [tauSorted, sortIdx] = sort(tauVec(:), 'ascend');
            if hasInf
                kVec = [1 ./ tauSorted; 0];
                perm = [sortIdx(:); nDecay + 1];
            else
                kVec = 1 ./ tauSorted;
                perm = sortIdx(:);
            end
            n = numel(kVec);
            
            Bmat = zeros(n, n);
            for ii = 1:n
                if ii == 1, prodK = 1; else, prodK = prod(kVec(1:ii-1)); end
                for jj = 1:ii
                    denom = 1;
                    for mm = 1:ii
                        if mm ~= jj
                            denom = denom * (kVec(mm) - kVec(jj));
                        end
                    end
                    Bmat(ii, jj) = prodK / denom;
                end
            end
            
            DADSperm = DADS(:, perm);
            % Prefer left-divide-style mldivide on the transpose for numerical
            % stability: DADS = EADS * B  =>  EADS' = B' \ DADS'
            EADS = (Bmat.' \ DADSperm.').';
            if any(~isfinite(EADS(:)))
                EADS = DADSperm * pinv(Bmat);
            end
        end
        
        function [SAS, speciesConc, tauSorted] = computeSASsequential( ...
                app, D, t, tauVec, betaVec, stretchOn, hasInf, ...
                t0, fwhm, irfMode, tRange, seqModel)
            %COMPUTESASSEQUENTIAL  Species-associated spectra for the
            %   sequential model  S1 -> S2 -> ... -> Sn (-> ground) in
            %   which any step may be a stretched exponential.
            %
            %   A stretched step can be read in two ways (seqModel):
            %
            %   'age' (default): the clock of a species starts when the
            %   population ENTERS it. Species i has the survival function
            %       S_i(a) = exp(-(a / tau_i)^beta_i)     (beta = 1: exp)
            %   of the time a since entry. With I_i(t) the cumulative
            %   population that has entered species i up to time t,
            %       c_i(t)     = int S_i(t - s) dI_i(s)
            %       I_{i+1}(t) = I_i(t) - c_i(t)
            %   (whatever has left species i has entered species i+1),
            %   and I_1(t) is the integrated IRF (a step at t0 when no
            %   IRF is used).
            %
            %   't0': one clock for all steps, started by the excitation.
            %   Step i has the time-dependent rate coefficient
            %       k_i(x) = (beta_i / tau_i) (x / tau_i)^(beta_i - 1)
            %   at the time x since excitation,
            %       dc_i/dx = k_(i-1)(x) c_(i-1) - k_i(x) c_i,
            %   so population that entered species i at x = s is still
            %   there at x with probability
            %       exp(-[(x / tau_i)^beta_i - (s / tau_i)^beta_i]).
            %   With an IRF this delta-pulse response is convolved with
            %   the Gaussian (see seqClockConc).
            %
            %   The two readings are identical for a first-order step and
            %   for the first species, so they differ only when a species
            %   other than the first is stretched. In both, population is
            %   conserved at every step and the model reduces exactly to
            %   the Bateman (EADS) solution when no step is stretched.
            %   The optional terminal species (hasInf) does not decay: it
            %   collects everything that left the last step.
            %
            %   irfMode : 'none'  no IRF, the kinetics start at t0;
            %             any other value: Gaussian IRF included
            %   tRange  : optional [tmin tmax] -- the amplitude solve uses
            %             the same delay window as the fit ([] = all).
            %   seqModel: 'age' (default) | 't0', see above
            %
            %   D: M x N data matrix (wavelength x delay).
            %   Returns:
            %     SAS         : M x n_species (NaN where a spectrum cannot
            %                   be determined, e.g. a masked wavelength)
            %     speciesConc : N x n_species concentration profiles
            %     tauSorted   : tau sorted ascending (= species order)

            t  = t(:);
            nT = numel(t);
            [tauSorted, sortIdx] = sort(tauVec(:), 'ascend');
            betaVec   = betaVec(:);
            stretchOn = logical(stretchOn(:));
            betaSorted    = betaVec(sortIdx);
            stretchSorted = stretchOn(sortIdx);
            nDecay   = numel(tauSorted);
            nSpecies = nDecay + double(hasInf);
            if nargin < 11, tRange = []; end
            if nargin < 12 || isempty(seqModel), seqModel = 'age'; end

            sigma  = fwhm / (2 * sqrt(2 * log(2)));
            hasSig = isfinite(sigma) && sigma > 0;
            useIRF = hasSig && ~strcmpi(irfMode, 'none');
            % The two readings differ only if a later species is stretched.
            clockModel = strcmpi(seqModel, 't0') && nDecay > 1 && ...
                any(stretchSorted(2:end));

            SAS = nan(size(D, 1), nSpecies);
            speciesConc = zeros(nT, nSpecies);
            if nDecay < 1 || any(~isfinite(tauSorted)) || ...
                    any(tauSorted <= 0)
                return;
            end

            if clockModel
                % One clock for all steps, started by the excitation.
                speciesConc = app.seqClockConc(t, tauSorted, betaSorted, ...
                    stretchSorted, hasInf, t0, fwhm, useIRF);
            else
                % ---- Time grid: fine around t0 (where the kinetics are
                % fast), geometrically coarser later, and containing every
                % data delay so that no interpolation is needed. ----
                if useIRF
                    tStart = t0 - 6 * sigma;
                    h0 = min(0.02 * tauSorted(1), 0.1 * sigma);
                    h0 = max(h0, 6 * sigma / 1e4);
                else
                    tStart = t0;
                    h0 = 0.02 * tauSorted(1);
                end
                tEnd = max(t);
                if ~(tEnd > tStart), return; end    % all data before the pump
                growth = 0.01;                      % 1 % step growth
                xEnd = tEnd - t0;
                xUniEnd = h0 / growth;              % uniform part ends here
                xUni = ((tStart - t0) : h0 : min(xUniEnd, xEnd)).';
                if xEnd > xUniEnd
                    nGeo = ceil(log(xEnd / xUniEnd) / log(1 + growth));
                    xGeo = xUniEnd * (1 + growth) .^ (1:nGeo).';
                else
                    xGeo = zeros(0, 1);
                end
                inRange = t >= tStart;
                tAll = unique([t0 + xUni; t0 + xGeo; t(inRange)]);
                tAll = tAll(tAll >= tStart);
                [~, loc] = ismember(t(inRange), tAll);
                nG  = numel(tAll);
                x   = tAll - t0;
                dtG = diff(tAll);
                tMid = tAll(1:end-1) + dtG / 2;

                % ---- Chain of cumulative arrivals ----
                if useIRF
                    Icum = app.expIRFConv(tAll, Inf, t0, fwhm);  % int. IRF
                    Icum = Icum(:);
                else
                    Icum = ones(nG, 1);                          % step at t0
                end
                c = zeros(nG, nSpecies);
                for ii = 1:nDecay
                    tauI  = tauSorted(ii);
                    betaI = betaSorted(ii);
                    if ii == 1
                        % First species: the model function itself (the
                        % same expression as the fit basis).
                        if useIRF && stretchSorted(1)
                            ci = app.stretchedIRFConv(tAll, tauI, betaI, ...
                                t0, fwhm);
                        elseif useIRF
                            ci = app.expIRFConv(tAll, tauI, t0, fwhm);
                        elseif stretchSorted(1)
                            ci = exp(-(max(x, 0) ./ tauI) .^ betaI);
                        else
                            ci = exp(-max(x, 0) ./ tauI);
                        end
                        ci = ci(:);
                    else
                        dI = diff(Icum);     % arrivals in each grid cell
                        ci = zeros(nG, 1);
                        if stretchSorted(ii)
                            % Arrivals of each cell times the survival at
                            % the age of the cell mid-point.
                            for nn = 2:nG
                                age = tAll(nn) - tMid(1:nn-1);
                                ci(nn) = sum(dI(1:nn-1) .* ...
                                    exp(-(age ./ tauI) .^ betaI));
                            end
                        else
                            % First-order decay: exact update for arrivals
                            % spread evenly over the cell.
                            z   = dtG ./ tauI;
                            dec = exp(-z);
                            phi = ones(size(z));
                            big = z > 1e-8;
                            phi(big) = -expm1(-z(big)) ./ z(big);
                            for nn = 2:nG
                                ci(nn) = ci(nn-1) * dec(nn-1) + ...
                                    dI(nn-1) * phi(nn-1);
                            end
                        end
                    end
                    ci(~isfinite(ci)) = 0;
                    c(:, ii) = ci;
                    Icum = Icum - ci;        % what has left species ii
                end
                if hasInf
                    c(:, nSpecies) = Icum;
                end
                speciesConc(inRange, :) = c(loc, :);
            end
            speciesConc(~isfinite(speciesConc)) = 0;

            % ---- Amplitude solve ----
            % Drop columns that ended up identically zero (numerical
            % degeneracy) -- they would make the LS solve singular.
            colNorms = sqrt(sum(speciesConc.^2, 1));
            keepCols = colNorms > 1e-12 * max(colNorms);
            nK = nnz(keepCols);
            % Delays used: the same window as the fit.
            tSel = true(nT, 1);
            if ~isempty(tRange) && all(isfinite(tRange))
                tSel = (t >= min(tRange)) & (t <= max(tRange));
                if ~any(tSel), tSel = true(nT, 1); end
            end
            if nK == 0 || nnz(tSel) <= nK, return; end

            Cuse  = speciesConc(:, keepCols);
            Dsel  = D(:, tSel);
            Skept = nan(size(D, 1), nK);
            if all(isfinite(Dsel(:)))
                Skept = solveLS(Cuse(tSel, :), Dsel.').';
            else
                % Gaps (chirp-corrected edges) or masked wavelengths:
                % solve every wavelength on its own finite delays, so
                % one masked row does not take the others with it.
                for ii = 1:size(D, 1)
                    m = tSel & isfinite(D(ii, :)).';
                    if nnz(m) <= nK, continue; end
                    Skept(ii, :) = solveLS(Cuse(m, :), D(ii, m).').';
                end
            end
            % Degenerate species columns stay NaN (the legend still
            % lists them).
            SAS(:, keepCols) = Skept;

            function X = solveLS(Cm, Y)
                try
                    X = lsqminnorm(Cm, Y);
                catch
                    X = pinv(Cm) * Y;
                end
                if any(~isfinite(X(:)))
                    X = pinv(Cm) * Y;
                end
            end
        end
        
        function conc = seqClockConc(app, t, tauS, betaS, strS, hasInf, ...
                t0, fwhm, useIRF)
            %SEQCLOCKCONC  Populations of the sequential model in which
            %   all steps run on ONE clock started by the excitation
            %   (seqModel 't0' of computeSASsequential): step i has the
            %   rate coefficient k_i(x) = d/dx (x / tau_i)^beta_i at the
            %   time x since excitation.
            %
            %   t                 : delays
            %   tauS, betaS, strS : per step, in species order (tau
            %                       ascending); strS = stretched flags
            %   useIRF            : convolve with the Gaussian IRF
            %   conc              : numel(t) x n_species
            %
            %   The response to a delta pulse is integrated on a
            %   geometric grid in x (1 % steps, starting far below the
            %   fastest time constant so that an x^beta onset is
            %   resolved). With G_i(x) = (x / tau_i)^beta_i and dI the
            %   population arriving in a grid cell,
            %       h_i(x_n) = h_i(x_n-1) exp(-[G_i(x_n) - G_i(x_n-1)])
            %                + dI exp(-[G_i(x_n) - G_i(x_mid)]).
            %   With an IRF that response (piecewise linear on the grid)
            %   is convolved with the Gaussian; the first species uses
            %   the same expression as the fit basis and the terminal
            %   species follows from conservation.
            t      = t(:);
            nT     = numel(t);
            nDecay = numel(tauS);
            nSp    = nDecay + double(hasInf);
            conc   = zeros(nT, nSp);
            xd     = t - t0;                 % delays on the pump clock
            if useIRF
                sigma = fwhm / (2 * sqrt(2 * log(2)));
                xEnd  = max(xd) + 6 * sigma;
            else
                sigma = 0;
                xEnd  = max(xd);
            end
            if ~(xEnd > 0), return; end      % all data before the pump
            
            % ---- Grid: 0, geometric steps up to xEnd, and the data
            % delays (so that nothing is interpolated without an IRF).
            growth = 0.01;
            xLo  = 1e-40 * tauS(1);
            nGeo = ceil(log(xEnd / xLo) / log(1 + growth));
            xg = xLo * (1 + growth) .^ (0:nGeo).';
            xg = unique([0; xg; xd(xd > 0)]);
            nG = numel(xg);
            dx = diff(xg);
            xm = xg(1:end-1) + dx / 2;
            
            % ---- Response to a delta pulse at x = 0 ----
            H = zeros(nG, nSp);
            Icum = ones(nG, 1);              % arrived in species 1
            for ii = 1:nDecay
                tauI  = tauS(ii);
                betaI = betaS(ii);
                if ii == 1
                    if strS(1)
                        hi = exp(-(xg ./ tauI) .^ betaI);
                    else
                        hi = exp(-xg ./ tauI);
                    end
                else
                    dI = diff(Icum);         % arrivals in each grid cell
                    if strS(ii)
                        G   = (xg ./ tauI) .^ betaI;
                        dec = exp(-diff(G));
                        arr = dI .* exp(-(G(2:end) - ...
                            (xm ./ tauI) .^ betaI));
                    else
                        % First-order step: exact for arrivals spread
                        % evenly over the cell.
                        z   = dx ./ tauI;
                        dec = exp(-z);
                        phi = ones(size(z));
                        big = z > 1e-8;
                        phi(big) = -expm1(-z(big)) ./ z(big);
                        arr = dI .* phi;
                    end
                    hi = zeros(nG, 1);
                    for nn = 2:nG
                        hi(nn) = hi(nn-1) * dec(nn-1) + arr(nn-1);
                    end
                end
                hi(~isfinite(hi)) = 0;
                H(:, ii) = hi;
                Icum = Icum - hi;            % what has left species ii
            end
            if hasInf
                H(:, nSp) = Icum;
            end
            
            if ~useIRF
                pos = xd >= 0;
                [~, loc] = ismember(xd(pos), xg);
                conc(pos, :) = H(loc, :);
                return;
            end
            
            % ---- Convolution with the Gaussian IRF ----
            if strS(1)
                conc(:, 1) = app.stretchedIRFConv(t, tauS(1), ...
                    betaS(1), t0, fwhm);
            else
                conc(:, 1) = app.expIRFConv(t, tauS(1), t0, fwhm);
            end
            a  = xg(1:end-1);
            b  = xg(2:end);
            Ha = H(1:end-1, 2:nDecay);
            Hb = H(2:end,   2:nDecay);
            gN = 1 / (sigma * sqrt(2 * pi));
            s2 = sigma * sqrt(2);
            wide = dx >= 1e-3 * sigma;
            for n = 1:nT
                X = xd(n);
                if ~(X + 6 * sigma > 0), continue; end
                k1 = find(b > X - 6 * sigma, 1, 'first');
                k2 = find(a < X + 6 * sigma, 1, 'last');
                if isempty(k1) || isempty(k2) || k2 < k1, continue; end
                kk = (k1:k2).';
                va = X - a(kk);
                vb = X - b(kk);
                wA = zeros(numel(kk), 1);    % weight of h at the left node
                wB = zeros(numel(kk), 1);    % ... and at the right node
                % Cells not much narrower than the IRF: exact integral
                % of (linear h) x Gaussian.
                L = wide(kk);
                P = 0.5 * (erf(va(L) / s2) - erf(vb(L) / s2));
                Q = (va(L) .* P + sigma^2 * gN * ...
                    (exp(-va(L).^2 / (2 * sigma^2)) - ...
                     exp(-vb(L).^2 / (2 * sigma^2)))) ./ dx(kk(L));
                wB(L) = Q;
                wA(L) = P - Q;
                % Much narrower cells: mid-point value of the Gaussian
                % (the exact form would lose digits by cancellation).
                N = ~L;
                w = 0.5 * gN * dx(kk(N)) .* ...
                    exp(-(X - xm(kk(N))).^2 / (2 * sigma^2));
                wA(N) = w;
                wB(N) = w;
                conc(n, 2:nDecay) = wA.' * Ha(kk, :) + wB.' * Hb(kk, :);
            end
            if hasInf
                Phi = app.expIRFConv(t, Inf, t0, fwhm);
                conc(:, nSp) = Phi(:) - sum(conc(:, 1:nDecay), 2);
            end
            conc(~isfinite(conc)) = 0;
        end
        
        function sch = parseTargetScheme(~, txt)
            %PARSETARGETSCHEME  Text description -> kinetic scheme struct.
            %
            %   One statement per line (or separated by ';'). Compartment
            %   and parameter names are free.
            %     A -> B              first-order step, own rate parameter
            %     A -> B : k1         ... with a named rate parameter
            %     B -> C : k1         the same name = the same parameter
            %     B -> C : 2*k1       rate tied to a parameter (also k1/2)
            %     A -> B -> C         a chain is split into single steps
            %     A <-> B : kf, kb    reversible: two steps
            %     B -> 2 C            stoichiometric factor on the product
            %     A -> GS             loss to the ground state. GS, S0 and
            %                         0 are the dark sink: no spectrum, no
            %                         population is tracked
            %     species: A, B, C    declare compartments (fixes their
            %                         order; allows one without any step)
            %     init: A = 1         initial populations (default: the
            %                         first compartment = 1, the others 0)
            %     zero: B 700 900     spectrum of B = 0 in a spectral range
            %     equal: C = B        spectrum of C equals that of B
            %     equal: C = 0.5*B 600 700   ... proportional, in a range
            %     nonneg              all spectra >= 0 (nonneg: A, B for
            %                         some of them)
            %   Text after % or # is a comment.
            %
            %   Returns sch with fields
            %     species  cellstr, order of first appearance (no sink)
            %     nS       number of compartments
            %     params   cellstr, names of the rate parameters; nP
            %     steps    struct array: from, to (0 = sink), coef,
            %              par (index into params), mult (rate = mult *
            %              k_par), label ('A -> B'), expr ('2*k1')
            %     c0       nS x 1 initial populations
            %     zero     struct array: sp, lo, hi
            %     rel      struct array: dep, master, ratio, lo, hi
            %     nonneg   1 x nS logical
            %     sinkName name used for the ground state in the text
            if isstring(txt), txt = cellstr(txt); end
            if iscell(txt)
                txt = strjoin(cellfun(@(c) char(c), txt(:).', ...
                    'UniformOutput', false), newline);
            end
            txt = char(txt);
            stm = regexp(txt, '[\n\r;]+', 'split');

            sinkNames = {'gs', 's0', '0', 'ground'};
            sinkName  = '';
            species = {};
            params  = {};
            steps = struct('from', {}, 'to', {}, 'coef', {}, 'par', {}, ...
                'mult', {}, 'label', {}, 'expr', {});
            initTxt = {};  zeroTxt = {};  equalTxt = {};  nonnegTxt = {};
            nonnegAll = false;
            nAuto = 0;

            for ii = 1:numel(stm)
                s = strtrim(regexprep(stm{ii}, '[%#].*$', ''));
                if isempty(s), continue; end
                % Keyword lines. The colon is required, so a compartment
                % may be called e.g. "initial".
                tkKey = regexp(s, ...
                    '^(?i)(init|zero|equal|nonneg|species)\s*:\s*(.*)$', ...
                    'tokens', 'once');
                if ~isempty(tkKey)
                    arg = strtrim(tkKey{2});
                    switch lower(tkKey{1})
                        case 'init'
                            initTxt{end+1} = arg; %#ok<AGROW>
                        case 'zero'
                            zeroTxt{end+1} = arg; %#ok<AGROW>
                        case 'equal'
                            equalTxt{end+1} = arg; %#ok<AGROW>
                        case 'nonneg'
                            if isempty(arg)
                                nonnegAll = true;
                            else
                                nonnegTxt{end+1} = arg; %#ok<AGROW>
                            end
                        case 'species'
                            nmL = regexp(arg, '\s*,\s*', 'split');
                            for pp = 1:numel(nmL)
                                if isempty(nmL{pp}), continue; end
                                checkName(nmL{pp}, s);
                                if isSink(nmL{pp})
                                    error(['Scheme line "%s": the ground ' ...
                                        'state (%s) is not a compartment.'], ...
                                        s, nmL{pp});
                                end
                                speciesIndex(nmL{pp});
                            end
                    end
                    continue;
                end
                if strcmpi(s, 'nonneg')
                    nonnegAll = true;
                    continue;
                end

                % Reaction line, optionally followed by ": rate labels".
                rateLab = {};
                rx = s;
                ic = find(s == ':', 1);
                if ~isempty(ic)
                    rx = strtrim(s(1:ic-1));
                    rateLab = regexp(strtrim(s(ic+1:end)), '\s*,\s*', 'split');
                    rateLab = rateLab(~cellfun(@isempty, rateLab));
                end
                [parts, arrows] = regexp(rx, '<->|->', 'split', 'match');
                if isempty(arrows)
                    error('Scheme line "%s": expected "A -> B".', s);
                end
                nHere = numel(arrows) + nnz(strcmp(arrows, '<->'));
                if ~isempty(rateLab) && numel(rateLab) ~= nHere
                    error(['Scheme line "%s": %d step(s) but %d rate ' ...
                        'label(s) after ":".'], s, nHere, numel(rateLab));
                end
                nm = cell(1, numel(parts));
                cf = ones(1, numel(parts));
                for pp = 1:numel(parts)
                    [nm{pp}, cf(pp)] = parseTerm(parts{pp}, s);
                end
                iLab = 0;
                for aa = 1:numel(arrows)
                    src = nm{aa};  dst = nm{aa + 1};
                    if isSink(src)
                        error(['Scheme line "%s": the ground state ' ...
                            '(%s) cannot be the source of a step.'], s, src);
                    end
                    if strcmp(src, dst)
                        error('Scheme line "%s": a step needs two different compartments.', s);
                    end
                    iLab = iLab + 1;
                    addStep(src, dst, cf(aa + 1), pickLab(rateLab, iLab), s);
                    if strcmp(arrows{aa}, '<->')
                        if isSink(dst)
                            error(['Scheme line "%s": a reversible step ' ...
                                'cannot involve the ground state.'], s);
                        end
                        if cf(aa + 1) ~= 1
                            error(['Scheme line "%s": no stoichiometric ' ...
                                'factor on a reversible step.'], s);
                        end
                        iLab = iLab + 1;
                        addStep(dst, src, 1, pickLab(rateLab, iLab), s);
                    end
                end
            end
            nS = numel(species);
            if nS < 1
                error('The scheme has no compartments. Example:  A -> B');
            end

            % ---- name the parameters of unlabelled steps: k1, k2, ... ----
            nxt = 1;
            for pp = 1:numel(params)
                if params{pp}(1) ~= '#', continue; end
                while any(strcmp(params, sprintf('k%d', nxt)))
                    nxt = nxt + 1;
                end
                params{pp} = sprintf('k%d', nxt);
            end
            for qq = 1:numel(steps)
                if steps(qq).mult == 1
                    steps(qq).expr = params{steps(qq).par};
                else
                    steps(qq).expr = sprintf('%.10g*%s', steps(qq).mult, ...
                        params{steps(qq).par});
                end
            end

            % ---- initial populations ----
            c0 = zeros(nS, 1);
            if isempty(initTxt)
                c0(1) = 1;
            else
                items = regexp(strjoin(initTxt, ','), '\s*,\s*', 'split');
                for ii = 1:numel(items)
                    it = strtrim(items{ii});
                    if isempty(it), continue; end
                    tk = regexp(it, '^(\S+?)\s*=\s*(\S+)$', 'tokens', 'once');
                    if isempty(tk)
                        error('init: cannot read "%s" (use  init: A = 1).', it);
                    end
                    idx = find(strcmp(species, tk{1}), 1);
                    val = str2double(tk{2});
                    if isempty(idx)
                        error('init: "%s" is not a compartment of the scheme.', tk{1});
                    end
                    if ~isfinite(val) || val < 0
                        error('init: the population of %s must be >= 0.', tk{1});
                    end
                    c0(idx) = val;
                end
                if ~any(c0 > 0)
                    error('init: at least one compartment needs a population > 0.');
                end
            end

            % ---- zero constraints ----
            zr = struct('sp', {}, 'lo', {}, 'hi', {});
            for ii = 1:numel(zeroTxt)
                z = zeroTxt{ii};
                tk = regexp(z, '^(\S+)\s+(.*)$', 'tokens', 'once');
                if isempty(tk)
                    error('zero: cannot read "%s" (use  zero: B 700 900).', z);
                end
                idx = find(strcmp(species, tk{1}), 1);
                if isempty(idx)
                    error('zero: "%s" is not a compartment of the scheme.', tk{1});
                end
                v = readRange(tk{2});
                if numel(v) ~= 2
                    error('zero: "%s" needs two numbers (from, to).', z);
                end
                zr(end+1) = struct('sp', idx, 'lo', min(v), 'hi', max(v)); %#ok<AGROW>
            end

            % ---- equal / proportional spectra ----
            rl = struct('dep', {}, 'master', {}, 'ratio', {}, 'lo', {}, 'hi', {});
            for ii = 1:numel(equalTxt)
                q = equalTxt{ii};
                ie = find(q == '=', 1);
                tk = {};
                if ~isempty(ie)
                    rest = regexprep(strtrim(q(ie+1:end)), '\s*\*\s*', '*');
                    tk = regexp(rest, '^(\S+)\s*(.*)$', 'tokens', 'once');
                end
                if isempty(tk)
                    error(['equal: cannot read "%s" (use  equal: C = B  ' ...
                        'or  equal: C = 0.5*B 600 700).'], q);
                end
                depN = strtrim(q(1:ie-1));
                masN = tk{1};
                rt = 1;
                tm = regexp(masN, ['^([-+]?[0-9]*\.?[0-9]+' ...
                    '(?:[eE][-+]?[0-9]+)?)\*(.+)$'], 'tokens', 'once');
                if ~isempty(tm)
                    rt = str2double(tm{1});
                    masN = tm{2};
                end
                iD = find(strcmp(species, depN), 1);
                iM = find(strcmp(species, masN), 1);
                if isempty(iD) || isempty(iM)
                    error('equal: "%s" names a compartment that is not in the scheme.', q);
                end
                if iD == iM
                    error('equal: "%s" relates a compartment to itself.', q);
                end
                if ~isfinite(rt) || rt == 0
                    error('equal: the factor in "%s" must be a non-zero number.', q);
                end
                lo = -Inf;  hi = Inf;
                if ~isempty(strtrim(tk{2}))
                    v = readRange(tk{2});
                    if numel(v) ~= 2
                        error('equal: the range in "%s" needs two numbers.', q);
                    end
                    lo = min(v);  hi = max(v);
                end
                for kk = 1:numel(rl)
                    if rl(kk).dep == iD && max(rl(kk).lo, lo) <= min(rl(kk).hi, hi)
                        error(['equal: "%s" is tied twice in overlapping ' ...
                            'ranges.'], species{iD});
                    end
                end
                rl(end+1) = struct('dep', iD, 'master', iM, 'ratio', rt, ...
                    'lo', lo, 'hi', hi); %#ok<AGROW>
            end
            % No circular definitions (C = B together with B = C).
            for kk = 1:numel(rl)
                cur = rl(kk).master;
                for hop = 1:nS
                    nx = find([rl.dep] == cur, 1);
                    if isempty(nx), break; end
                    cur = rl(nx).master;
                    if cur == rl(kk).dep
                        error('equal: the constraints on "%s" are circular.', ...
                            species{rl(kk).dep});
                    end
                end
            end

            % ---- non-negative spectra ----
            nn = false(1, nS);
            if nonnegAll, nn(:) = true; end
            for ii = 1:numel(nonnegTxt)
                nmL = regexp(nonnegTxt{ii}, '\s*,\s*', 'split');
                for pp = 1:numel(nmL)
                    if isempty(nmL{pp}), continue; end
                    idx = find(strcmp(species, nmL{pp}), 1);
                    if isempty(idx)
                        error('nonneg: "%s" is not a compartment of the scheme.', nmL{pp});
                    end
                    nn(idx) = true;
                end
            end

            if isempty(sinkName), sinkName = 'GS'; end
            sch = struct('species', {species}, 'nS', nS, ...
                'params', {params}, 'nP', numel(params), 'steps', steps, ...
                'c0', c0, 'zero', zr, 'rel', rl, 'nonneg', nn, ...
                'sinkName', sinkName);

            % ---------------- nested helpers ----------------
            function tf = isSink(name)
                tf = any(strcmpi(name, sinkNames));
            end

            function r = pickLab(labs, k)
                if isempty(labs), r = ''; else, r = labs{k}; end
            end

            function v = readRange(str)
                r2 = strrep(str, '..', ' ');
                r2 = regexprep(r2, '(?<=[0-9.])\s*-\s*(?=[0-9.])', ' ');
                r2 = strrep(r2, ',', ' ');
                v = sscanf(r2, '%f');
                if any(~isfinite(v)), v = []; end
            end

            function checkName(name, lineTxt)
                if isempty(name) || ...
                        isempty(regexp(name, '^[A-Za-z0-9_()''^{}]+$', 'once'))
                    error(['Scheme line "%s": "%s" is not a valid ' ...
                        'compartment name (letters, digits, _ ( ) '' only).'], ...
                        lineTxt, name);
                end
                if ~isSink(name) && ...
                        ~isempty(regexp(name, '^[0-9.]+$', 'once'))
                    error(['Scheme line "%s": a compartment name cannot ' ...
                        'be a number ("%s").'], lineTxt, name);
                end
            end

            function [name, coef] = parseTerm(term, lineTxt)
                term = strtrim(term);
                coef = 1;
                tk2 = regexp(term, '^([0-9]*\.?[0-9]+)\s*\*\s*(\S+)$', ...
                    'tokens', 'once');
                if isempty(tk2)
                    tk2 = regexp(term, '^([0-9]*\.?[0-9]+)\s+(\S+)$', ...
                        'tokens', 'once');
                end
                if ~isempty(tk2)
                    coef = str2double(tk2{1});
                    name = tk2{2};
                else
                    name = term;
                end
                checkName(name, lineTxt);
                if ~isfinite(coef) || coef <= 0
                    error('Scheme line "%s": bad factor in "%s".', ...
                        lineTxt, term);
                end
            end

            function [pname, mult] = parseRate(str, lineTxt)
                str = strtrim(str);
                num = '([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)';
                nam = '([A-Za-z][A-Za-z0-9_]*)';
                mult = 1;
                t1 = regexp(str, ['^' nam '$'], 'tokens', 'once');
                t2 = regexp(str, ['^' num '\s*\*\s*' nam '$'], 'tokens', 'once');
                t3 = regexp(str, ['^' nam '\s*\*\s*' num '$'], 'tokens', 'once');
                t4 = regexp(str, ['^' nam '\s*/\s*' num '$'], 'tokens', 'once');
                if ~isempty(t1)
                    pname = t1{1};
                elseif ~isempty(t2)
                    pname = t2{2};  mult = str2double(t2{1});
                elseif ~isempty(t3)
                    pname = t3{1};  mult = str2double(t3{2});
                elseif ~isempty(t4)
                    pname = t4{1};  mult = 1 / str2double(t4{2});
                else
                    error(['Scheme line "%s": cannot read the rate "%s" ' ...
                        '(use  k1,  2*k1  or  k1/2).'], lineTxt, str);
                end
                if ~isfinite(mult) || mult <= 0
                    error(['Scheme line "%s": the factor in "%s" must ' ...
                        'be positive.'], lineTxt, str);
                end
            end

            function k = speciesIndex(name)
                k = find(strcmp(species, name), 1);
                if isempty(k)
                    species{end+1} = name;
                    k = numel(species);
                end
            end

            function addStep(src, dst, coef, rateTxt, lineTxt)
                iFrom = speciesIndex(src);
                if isSink(dst)
                    iTo = 0;
                    if isempty(sinkName), sinkName = dst; end
                else
                    iTo = speciesIndex(dst);
                end
                for qs = 1:numel(steps)
                    if steps(qs).from == iFrom && steps(qs).to == iTo
                        error('The step "%s -> %s" appears twice.', src, dst);
                    end
                end
                if isempty(rateTxt)
                    nAuto = nAuto + 1;
                    pname = sprintf('#%d', nAuto);
                    mlt = 1;
                else
                    [pname, mlt] = parseRate(rateTxt, lineTxt);
                end
                ip = find(strcmp(params, pname), 1);
                if isempty(ip)
                    params{end+1} = pname;
                    ip = numel(params);
                end
                if coef == 1
                    lab = sprintf('%s -> %s', src, dst);
                else
                    lab = sprintf('%s -> %g %s', src, coef, dst);
                end
                steps(end+1) = struct('from', iFrom, 'to', iTo, ...
                    'coef', coef, 'par', ip, 'mult', mlt, ...
                    'label', lab, 'expr', '');
            end
        end

        function txtL = targetSchemeText(~, sch)
            %TARGETSCHEMETEXT  Scheme struct -> canonical text (cellstr).
            %   The inverse of parseTargetScheme, used when the scheme is
            %   edited through the K-matrix table. Every step is written
            %   with its rate expression, so parameter names survive.
            %   Only species, steps (from, to, coef, expr), c0, zero, rel,
            %   nonneg and sinkName are read.
            sp = sch.species;
            txtL = {['species: ' strjoin(sp, ', ')]};
            for jj = 1:numel(sch.steps)
                st = sch.steps(jj);
                if st.to == 0
                    dst = sch.sinkName;
                else
                    dst = sp{st.to};
                end
                if st.coef ~= 1
                    dst = sprintf('%g %s', st.coef, dst);
                end
                txtL{end+1} = sprintf('%s -> %s : %s', sp{st.from}, dst, st.expr); %#ok<AGROW>
            end
            c0 = sch.c0(:);
            def = zeros(numel(sp), 1);  def(1) = 1;
            if ~isequal(c0, def)
                it = {};
                for jj = 1:numel(sp)
                    if c0(jj) ~= 0
                        it{end+1} = sprintf('%s = %.10g', sp{jj}, c0(jj)); %#ok<AGROW>
                    end
                end
                txtL{end+1} = ['init: ' strjoin(it, ', ')];
            end
            for jj = 1:numel(sch.zero)
                z = sch.zero(jj);
                txtL{end+1} = sprintf('zero: %s %.10g %.10g', sp{z.sp}, z.lo, z.hi); %#ok<AGROW>
            end
            for jj = 1:numel(sch.rel)
                r = sch.rel(jj);
                if r.ratio == 1
                    rhs = sp{r.master};
                else
                    rhs = sprintf('%.10g*%s', r.ratio, sp{r.master});
                end
                ln = sprintf('equal: %s = %s', sp{r.dep}, rhs);
                if isfinite(r.lo) && isfinite(r.hi)
                    ln = sprintf('%s %.10g %.10g', ln, r.lo, r.hi);
                end
                txtL{end+1} = ln; %#ok<AGROW>
            end
            if any(sch.nonneg)
                if all(sch.nonneg)
                    txtL{end+1} = 'nonneg';
                else
                    txtL{end+1} = ['nonneg: ' strjoin(sp(sch.nonneg), ', ')];
                end
            end
            txtL = txtL(:);
        end

        function D = targetDiagramLayout(app, sch, W, H, fs0, tauStep, undet, unitTxt)
            %TARGETDIAGRAMLAYOUT  Box-and-arrow picture of a kinetic scheme.
            %   Geometry only: returns the rectangles, lines, arrow heads
            %   and texts of the picture in the pixel coordinates of a
            %   W x H canvas (origin bottom-left), so that drawing it is
            %   a plain loop and the layout can be tested without a
            %   figure.
            %
            %   sch       scheme struct of parseTargetScheme
            %   fs0       largest font size (pixels); it is reduced until
            %             the picture fits the canvas (not below 6)
            %   tauStep   time constant of every step ([] = none shown)
            %   undet     true for steps whose time constant the data do
            %             not determine (marked "?")
            %   unitTxt   time unit printed after the time constants
            %
            %   Layout: compartments stand in columns along the flow
            %   (longest path from the populated compartments; a step
            %   that closes a cycle points back), separate groups of
            %   compartments side by side. A step that skips a column
            %   runs through a free row of the columns in between. A loss
            %   to the ground state is a short arrow below its box.
            %
            %   D.W, D.H   canvas used (larger than W x H only when the
            %              picture does not fit at font size 6)
            %   D.fs       font size used;  D.fits  true if it fits
            %   D.rect     pos [x y w h], sp (compartment), init
            %   D.line     x, y (polyline)
            %   D.head     x, y (filled triangle)
            %   D.text     x, y, str, ha, va, size, kind, sp, bg, bb
            %              kind: 'name' | 'annot' | 'rate' | 'value' |
            %              'sink' | 'note';  bg = draw on a white patch;
            %              bb = estimated extent [x0 y0 x1 y1]
            %   D.col, D.row   column and row coordinate of every
            %              compartment;  D.noteY, D.noteSize  baseline
            %              and font size of the note line at the bottom
            if nargin < 5 || isempty(fs0), fs0 = 11; end
            if nargin < 6, tauStep = []; end
            if nargin < 7, undet = []; end
            if nargin < 8, unitTxt = ''; end
            nS  = sch.nS;
            st  = sch.steps;
            nSt = numel(st);
            hasVal = nSt > 0 && numel(tauStep) == nSt;
            if numel(undet) ~= nSt, undet = false(1, nSt); end
            eRect = struct('pos', {}, 'sp', {}, 'init', {});
            eLine = struct('x', {}, 'y', {});
            eText = struct('x', {}, 'y', {}, 'str', {}, 'ha', {}, 'va', {}, ...
                'size', {}, 'kind', {}, 'sp', {}, 'bg', {}, 'bb', {});
            D = struct('W', W, 'H', H, 'fs', fs0, 'fits', true, ...
                'rect', eRect, 'line', eLine, 'head', eLine, 'text', eText, ...
                'col', zeros(1, 0), 'row', zeros(1, 0), ...
                'noteY', 0, 'noteSize', fs0);
            if nS < 1, return; end

            fromV = zeros(1, nSt);  toV = zeros(1, nSt);
            for q = 1:nSt
                fromV(q) = st(q).from;  toV(q) = st(q).to;
            end
            c0 = sch.c0(:).';
            sinkName = 'GS';
            if isfield(sch, 'sinkName') && ~isempty(sch.sinkName)
                sinkName = sch.sinkName;
            end

            % ---- Text of every step: rate; time constant, branching ----
            lab = cell(1, nSt);
            for q = 1:nSt
                s1 = st(q).expr;
                if st(q).coef ~= 1
                    s1 = sprintf('%s [x%g]', s1, st(q).coef);
                end
                lab{q} = {s1};
                if hasVal
                    tv = tauStep(q);
                    if ~(isfinite(tv) && tv > 0)
                        s2 = '?';
                    elseif tv >= 1e3 && tv < 1e5
                        s2 = strtrim(sprintf('%.0f %s', tv, unitTxt));
                    else
                        s2 = strtrim(sprintf('%.3g %s', tv, unitTxt));
                    end
                    if undet(q), s2 = [s2 ' ?']; end %#ok<AGROW>
                    outI = find(fromV == fromV(q));
                    kOut = 1 ./ tauStep(outI);
                    if numel(outI) > 1 && all(isfinite(kOut)) && sum(kOut) > 0
                        s2 = sprintf('%s (%.0f%%)', s2, ...
                            100 / tauStep(q) / sum(kOut));
                    end
                    lab{q}{2} = s2;
                end
            end

            % ---- Note under the name of a compartment ----
            ann = cell(1, nS);
            defC0 = zeros(1, nS);  defC0(1) = 1;
            showInit = ~isequal(c0, defC0);
            for q = 1:nS
                it = {};
                if showInit && c0(q) > 0
                    it{end+1} = sprintf('init %g', c0(q)); %#ok<AGROW>
                end
                nz = 0;
                for w = 1:numel(sch.zero)
                    if sch.zero(w).sp ~= q, continue; end
                    nz = nz + 1;
                    if nz == 1
                        it{end+1} = sprintf('zero %g-%g', ...
                            sch.zero(w).lo, sch.zero(w).hi); %#ok<AGROW>
                    elseif nz == 2
                        it{end} = [it{end} ' ..'];
                    end
                end
                for w = 1:numel(sch.rel)
                    if sch.rel(w).dep ~= q, continue; end
                    if sch.rel(w).ratio == 1
                        it{end+1} = ['= ' sch.species{sch.rel(w).master}]; %#ok<AGROW>
                    else
                        it{end+1} = sprintf('= %g %s', sch.rel(w).ratio, ...
                            sch.species{sch.rel(w).master}); %#ok<AGROW>
                    end
                end
                if sch.nonneg(q), it{end+1} = 'nonneg'; end %#ok<AGROW>
                a = strjoin(it, ', ');
                if numel(a) > 18, a = [a(1:16) '..']; end
                ann{q} = a;
            end

            % ---- Graph of the compartments (the sink is not a node) ----
            A = false(nS);
            for q = 1:nSt
                if toV(q) > 0, A(fromV(q), toV(q)) = true; end
            end
            U = A | A.';
            comp = zeros(1, nS);  nC = 0;
            for s0 = 1:nS
                if comp(s0) > 0, continue; end
                nC = nC + 1;
                comp(s0) = nC;
                stk = s0;
                while ~isempty(stk)
                    u = stk(end);  stk(end) = [];
                    nb = find(U(u, :) & comp == 0);
                    comp(nb) = nC;
                    stk = [stk, nb]; %#ok<AGROW>
                end
            end

            % Steps that close a cycle (depth-first search starting at
            % the populated compartments) are left out of the ordering.
            Adag = A;
            state = zeros(1, nS);            % 0 new, 1 open, 2 done
            order = [find(c0 > 0), find(~(c0 > 0))];
            for r = order
                if state(r) ~= 0, continue; end
                stk = r;  nxt = 1;
                state(r) = 1;
                while ~isempty(stk)
                    u = stk(end);
                    nb = find(A(u, :));
                    k = nxt(end);
                    if k > numel(nb)
                        state(u) = 2;
                        stk(end) = [];  nxt(end) = [];
                    else
                        nxt(end) = k + 1;
                        v = nb(k);
                        if state(v) == 1
                            Adag(u, v) = false;
                        elseif state(v) == 0
                            state(v) = 1;
                            stk(end+1) = v; %#ok<AGROW>
                            nxt(end+1) = 1; %#ok<AGROW>
                        end
                    end
                end
            end

            % Column = longest path from a compartment nothing leads to.
            layer = zeros(1, nS);
            indeg = sum(Adag, 1);
            queue = find(indeg == 0);
            while ~isempty(queue)
                u = queue(1);  queue(1) = [];
                for v = find(Adag(u, :))
                    layer(v) = max(layer(v), layer(u) + 1);
                    indeg(v) = indeg(v) - 1;
                    if indeg(v) == 0, queue(end+1) = v; end %#ok<AGROW>
                end
            end
            % A start that only feeds a later column stands next to it.
            for u = find(sum(Adag, 1) == 0)
                nb = find(Adag(u, :));
                if ~isempty(nb), layer(u) = min(layer(nb)) - 1; end
            end

            % ---- Connections: one per pair of compartments ----
            E = struct('lo', {}, 'hi', {}, 'fwd', {}, 'bwd', {}, 'dum', {});
            for u = 1:nS
                for v = u+1:nS
                    if ~U(u, v), continue; end
                    if layer(u) <= layer(v)
                        e = struct('lo', u, 'hi', v);
                    else
                        e = struct('lo', v, 'hi', u);
                    end
                    e.fwd = find(fromV == e.lo & toV == e.hi);
                    e.bwd = find(fromV == e.hi & toV == e.lo);
                    e.dum = zeros(1, 0);
                    E(end+1) = e; %#ok<AGROW>
                end
            end
            nE = numel(E);
            eComp = zeros(1, nE);
            for k = 1:nE, eComp(k) = comp(E(k).lo); end

            % ---- Rows: every column is ordered by the mean row of what
            %      it is connected to in the column on its left ----
            ysp   = zeros(1, nS);            % row coordinate, up > 0
            colOf = zeros(1, nS);            % column (0-based), all groups
            nRow  = 1;
            colBase = 0;
            colComp = zeros(1, 0);           % group of every column
            colN    = zeros(1, 0);           % rows used in every column
            for c = 1:nC
                mem = find(comp == c);
                nL  = max(layer(mem)) + 1;
                ec  = find(eComp == c);
                for l = 0:nL-1
                    it = zeros(0, 4);        % [-key, isDummy, id, 0]
                    for s = mem(layer(mem) == l)
                        yp = zeros(1, 0);
                        for k = ec
                            if E(k).hi ~= s, continue; end
                            if layer(E(k).lo) == l - 1
                                yp(end+1) = ysp(E(k).lo); %#ok<AGROW>
                            elseif ~isempty(E(k).dum)
                                yp(end+1) = E(k).dum(end); %#ok<AGROW>
                            end
                        end
                        key = -Inf;
                        if ~isempty(yp), key = mean(yp); end
                        it(end+1, :) = [-key, 0, s, 0]; %#ok<AGROW>
                    end
                    for k = ec
                        if ~(layer(E(k).lo) < l && l < layer(E(k).hi)), continue; end
                        if layer(E(k).lo) == l - 1
                            key = ysp(E(k).lo);
                        else
                            key = E(k).dum(end);
                        end
                        it(end+1, :) = [-key, 1, k, 0]; %#ok<AGROW>
                    end
                    it = sortrows(it, [1 2 3]);
                    n  = size(it, 1);
                    for w = 1:n
                        yv = (n - 1) / 2 - (w - 1);
                        if it(w, 2) == 0
                            ysp(it(w, 3)) = yv;
                        else
                            E(it(w, 3)).dum(end+1) = yv;
                        end
                    end
                    nRow = max(nRow, n);
                    colN(colBase + l + 1) = n; %#ok<AGROW>
                end
                colOf(mem) = colBase + layer(mem);
                colComp(colBase + (1:nL)) = c; %#ok<AGROW>
                colBase = colBase + nL;
            end
            nCol = colBase;

            % A label that hangs below a line next to a box (the returning
            % step of a pair, or a line that runs downwards) keeps the
            % label of a loss on the other side of that box; with both
            % sides taken, the loss is labelled under the sink name.
            sideR = false(1, nS);  sideL = false(1, nS);
            for k = 1:nE
                yr = [ysp(E(k).lo), E(k).dum, ysp(E(k).hi)];
                isP = ~isempty(E(k).fwd) && ~isempty(E(k).bwd);
                if isP || yr(2) < yr(1), sideR(E(k).lo) = true; end
                if isP || yr(end-1) < yr(end), sideL(E(k).hi) = true; end
            end

            % Label lines of every connection and of every loss.
            eF = cell(1, nE);  eB = cell(1, nE);
            for k = 1:nE
                eF{k} = cell(1, 0);  eB{k} = cell(1, 0);
                if ~isempty(E(k).fwd), eF{k} = [lab{E(k).fwd}]; end
                if ~isempty(E(k).bwd), eB{k} = [lab{E(k).bwd}]; end
            end
            sinkOf = cell(1, nS);            % loss steps of a compartment
            sLab   = cell(1, nS);
            for q = 1:nS
                sinkOf{q} = find(fromV == q & toV == 0);
                sLab{q}   = cell(1, 0);
                if ~isempty(sinkOf{q}), sLab{q} = [lab{sinkOf{q}}]; end
            end
            anySink = any(~cellfun(@isempty, sinkOf));
            under   = false(1, nS);          % loss labelled under the sink name
            for q = 1:nS
                under(q) = ~isempty(sinkOf{q}) && sideR(q) && sideL(q);
            end
            anyAnn  = any(~cellfun(@isempty, ann));
            maxCh = 0;  maxLn = 0;  maxSCh = 0;  maxSLn = 0;  anyPair = false;
            for k = 1:nE
                both = [eF{k}, eB{k}];
                maxCh = max([maxCh, cellfun(@numel, both)]);
                maxLn = max([maxLn, numel(eF{k}), numel(eB{k})]);
                anyPair = anyPair || (~isempty(eF{k}) && ~isempty(eB{k}));
            end
            for q = 1:nS
                maxSCh = max([maxSCh, cellfun(@numel, sLab{q})]);
                maxSLn = max(maxSLn, numel(sLab{q}));
            end

            % ---- Geometry at font size fs; shrink until it fits ----
            fs = fs0;
            for pass = 1:6
                fsS = max(6, round(1.68 * fs) / 2);  % labels and notes
                cwN = 0.66 * fs;                     % bold names
                cwS = 0.56 * fsS;
                lh  = 1.22 * fsS;                    % line height of labels
                bh  = 1.9 * fs;
                if anyAnn, bh = bh + 1.15 * fsS; end
                colW = 3.0 * fs * ones(1, nCol);
                for q = 1:nS
                    wq = max([3.0 * fs, ...
                        numel(sch.species{q}) * cwN + 1.3 * fs, ...
                        numel(ann{q}) * cwS + 1.0 * fs]);
                    colW(colOf(q) + 1) = max(colW(colOf(q) + 1), wq);
                end
                gx  = max(2.6 * fs, max(maxCh, maxSCh) * cwS + 1.3 * fs);
                sOver = 0;                 % loss label beyond its box
                if anySink
                    sOver = max(0, 0.45 * fs + maxSCh * cwS - min(colW) / 2);
                end
                gxC = 1.6 * fs + sOver;    % between separate groups
                mS  = 0.8 * fs + sOver;
                stubLen = 0;  stubH = 0;
                if anySink
                    stubLen = max(1.7 * fs, maxSLn * lh + 0.6 * fs);
                    stubH   = stubLen + 1.3 * fsS;
                    if any(under), stubH = stubH + maxSLn * lh; end
                end
                labH = 0;
                if nE > 0, labH = maxLn * lh + 0.35 * fs; end
                vGap = 1.0 * fs;
                if nRow > 1, vGap = max(vGap, 0.5 * labH + 0.5 * fs); end
                ry  = bh + stubH + vGap;
                dOff = 0.34 * fs;          % half distance of a <-> pair
                noteH = lh + 0.5 * fs;
                mT  = 0.6 * fs + max(0, labH + dOff - bh / 2);
                mB  = 0.4 * fs + noteH + max(stubH, ...
                    anyPair * max(0, labH + dOff - bh / 2));
                xc  = zeros(1, nCol);
                x   = mS;
                for j = 1:nCol
                    if j > 1
                        if colComp(j) == colComp(j - 1), x = x + gx; else, x = x + gxC; end
                    end
                    xc(j) = x + colW(j) / 2;
                    x = x + colW(j);
                end
                Wn = x + mS;
                Hn = mB + nRow * bh + (nRow - 1) * (ry - bh) + mT;
                if (Wn <= W && Hn <= H) || fs <= 6, break; end
                sc = min([W / Wn, H / Hn, 0.97]);
                fs = max(6, floor(2 * fs * sc) / 2);
            end
            ox = max(0, (W - Wn) / 2);
            oy = max(0, (H - Hn) / 2);
            yOf = @(r) oy + mB + bh / 2 + ((nRow - 1) / 2 + r) * ry;
            xc = xc + ox;

            % ---- Boxes ----
            R = zeros(nS, 4);
            T = eText;
            for q = 1:nS
                wq = colW(colOf(q) + 1);
                R(q, :) = [xc(colOf(q) + 1) - wq / 2, yOf(ysp(q)) - bh / 2, wq, bh];
                D.rect(end+1) = struct('pos', R(q, :), 'sp', q, ...
                    'init', c0(q) > 0);
                cx = R(q, 1) + wq / 2;
                if isempty(ann{q})
                    T(end+1) = struct('x', cx, 'y', R(q, 2) + bh / 2, ...
                        'str', sch.species{q}, 'ha', 'center', 'va', 'middle', ...
                        'size', fs, 'kind', 'name', 'sp', q, 'bg', false, 'bb', []); %#ok<AGROW>
                else
                    T(end+1) = struct('x', cx, 'y', R(q, 2) + bh - 0.95 * fs, ...
                        'str', sch.species{q}, 'ha', 'center', 'va', 'middle', ...
                        'size', fs, 'kind', 'name', 'sp', q, 'bg', false, 'bb', []); %#ok<AGROW>
                    T(end+1) = struct('x', cx, 'y', R(q, 2) + 0.3 * fs + 0.5 * fsS, ...
                        'str', ann{q}, 'ha', 'center', 'va', 'middle', ...
                        'size', fsS, 'kind', 'annot', 'sp', q, 'bg', false, 'bb', []); %#ok<AGROW>
                end
            end

            % ---- Steps between compartments ----
            ah = 0.8 * fs;  aw = 0.3 * fs;         % arrow head
            for k = 1:nE
                lo = E(k).lo;  hi = E(k).hi;
                P = [R(lo, 1) + R(lo, 3) / 2, R(lo, 2) + bh / 2];
                for w = 1:numel(E(k).dum)
                    P(end+1, :) = [xc(colOf(lo) + 1 + w), yOf(E(k).dum(w))]; %#ok<AGROW>
                end
                P(end+1, :) = [R(hi, 1) + R(hi, 3) / 2, R(hi, 2) + bh / 2]; %#ok<AGROW>
                isPair = ~isempty(eF{k}) && ~isempty(eB{k});
                for dirn = 1:2
                    if dirn == 1, ln = eF{k}; else, ln = eB{k}; end
                    if isempty(ln), continue; end
                    sh = 0;
                    if isPair, sh = dOff * (3 - 2 * dirn); end   % +dOff, -dOff
                    Q = P;  Q(:, 2) = Q(:, 2) + sh;
                    rs = R(lo, :);  re = R(hi, :);
                    if dirn == 2
                        Q = flipud(Q);
                        rs = R(hi, :);  re = R(lo, :);
                    end
                    Q(1, :)   = app.tdBoxExit(Q(1, :), Q(2, :), rs, 0);
                    Q(end, :) = app.tdBoxExit(Q(end, :), Q(end-1, :), re, 0.12 * fs);
                    uv = Q(end, :) - Q(end-1, :);
                    nu = sqrt(sum(uv .^ 2));
                    if nu > 0, uv = uv / nu; else, uv = [1 0]; end
                    tip = Q(end, :);
                    bas = tip - ah * uv;
                    nv  = [-uv(2), uv(1)];
                    D.head(end+1) = struct( ...
                        'x', [tip(1), bas(1) + aw * nv(1), bas(1) - aw * nv(1)], ...
                        'y', [tip(2), bas(2) + aw * nv(2), bas(2) - aw * nv(2)]);
                    Q(end, :) = tip - 0.6 * ah * uv;
                    D.line(end+1) = struct('x', Q(:, 1).', 'y', Q(:, 2).');

                    % Label, centred at C: above a level line (below it
                    % for the returning step of a pair). Next to a
                    % slanted line it stands on the outer side of the
                    % fan the line belongs to, clear of the line. At the
                    % bend of a step that skips a column it stands below
                    % or above the bend, or on it (white patch).
                    nL1 = numel(ln);
                    hw  = max(cellfun(@numel, ln)) * cwS / 2;
                    hh  = nL1 * lh / 2;
                    onPatch = false;
                    side = 1;
                    if isPair && dirn == 2, side = -1; end
                    if isempty(E(k).dum)
                        M  = (Q(1, :) + tip) / 2;
                        dv = tip - Q(1, :);
                        if abs(dv(2)) < 0.2 * abs(dv(1))
                            C = M + side * [0, 0.3 * fs + hh];
                        else
                            if ~isPair
                                yRef = 0;
                                if colN(colOf(lo) + 1) < colN(colOf(hi) + 1)
                                    yRef = ysp(lo);
                                elseif colN(colOf(hi) + 1) < colN(colOf(lo) + 1)
                                    yRef = ysp(hi);
                                end
                                if (ysp(lo) + ysp(hi)) / 2 < yRef, side = -1; end
                            end
                            nv = [-dv(2), dv(1)] / sqrt(sum(dv .^ 2));
                            if nv(2) < 0, nv = -nv; end
                            C = M + side * (hw * abs(nv(1)) + hh * abs(nv(2)) ...
                                + 0.3 * fs) * nv;
                        end
                    else
                        iv = 1 + ceil(numel(E(k).dum) / 2);
                        V  = P(iv, :) + [0, sh];
                        isLow  = P(iv, 2) < P(iv - 1, 2) && P(iv, 2) < P(iv + 1, 2);
                        isHigh = P(iv, 2) > P(iv - 1, 2) && P(iv, 2) > P(iv + 1, 2);
                        if isPair
                            C = V + side * [0, 0.35 * fs + hh];
                            onPatch = true;
                        elseif isLow
                            C = V - [0, 0.35 * fs + hh];
                        elseif isHigh
                            C = V + [0, 0.35 * fs + hh];
                        else
                            C = V;
                            onPatch = true;
                        end
                    end
                    for w = 1:nL1
                        knd = 'rate';
                        if hasVal && mod(w, 2) == 0, knd = 'value'; end
                        T(end+1) = struct('x', C(1), ...
                            'y', C(2) + ((nL1 - 1) / 2 - (w - 1)) * lh, ...
                            'str', ln{w}, 'ha', 'center', 'va', 'middle', ...
                            'size', fsS, 'kind', knd, 'sp', 0, ...
                            'bg', onPatch, 'bb', []); %#ok<AGROW>
                    end
                end
            end

            % ---- Losses to the ground state ----
            for q = 1:nS
                if isempty(sinkOf{q}), continue; end
                cx = R(q, 1) + R(q, 3) / 2;
                yb = R(q, 2);
                tip = [cx, yb - stubLen];
                D.line(end+1) = struct('x', [cx, cx], 'y', [yb, tip(2) + 0.6 * ah]);
                D.head(end+1) = struct('x', [cx, cx + aw, cx - aw], ...
                    'y', [tip(2), tip(2) + ah, tip(2) + ah]);
                T(end+1) = struct('x', cx, 'y', tip(2) - 0.15 * fs, ...
                    'str', sinkName, 'ha', 'center', 'va', 'top', ...
                    'size', fsS, 'kind', 'sink', 'sp', 0, 'bg', false, 'bb', []); %#ok<AGROW>
                ln = sLab{q};
                for w = 1:numel(ln)
                    knd = 'rate';
                    if hasVal && mod(w, 2) == 0, knd = 'value'; end
                    ty = yb - 0.25 * fs - (w - 1) * lh;
                    if under(q)
                        tx = cx;  tha = 'center';
                        ty = tip(2) - 0.15 * fs - 1.3 * fsS - (w - 1) * lh;
                    elseif sideR(q)
                        tx = cx - 0.45 * fs;  tha = 'right';
                    else
                        tx = cx + 0.45 * fs;  tha = 'left';
                    end
                    T(end+1) = struct('x', tx, 'y', ty, ...
                        'str', ln{w}, 'ha', tha, 'va', 'top', 'size', fsS, ...
                        'kind', knd, 'sp', 0, 'bg', false, 'bb', []); %#ok<AGROW>
                end
            end

            % ---- Note line ----
            note = 'thick frame = initial population';
            if hasVal && any(undet)
                note = [note ';  ? = not determined'];
            end
            D.noteY = 0.4 * fs;
            D.noteSize = fsS;
            T(end+1) = struct('x', 0.6 * fs, 'y', D.noteY, 'str', note, ...
                'ha', 'left', 'va', 'bottom', 'size', fsS, 'kind', 'note', ...
                'sp', 0, 'bg', false, 'bb', []);

            % Estimated extent of every text (for tests of the layout).
            for w = 1:numel(T)
                if strcmp(T(w).kind, 'name'), cw = 0.66 * T(w).size; else, cw = 0.56 * T(w).size; end
                tw = numel(T(w).str) * cw;
                th = 1.15 * T(w).size;
                switch T(w).ha
                    case 'center', x0 = T(w).x - tw / 2;
                    case 'right',  x0 = T(w).x - tw;
                    otherwise,     x0 = T(w).x;
                end
                switch T(w).va
                    case 'middle', y0 = T(w).y - th / 2;
                    case 'top',    y0 = T(w).y - th;
                    otherwise,     y0 = T(w).y;
                end
                T(w).bb = [x0, y0, x0 + tw, y0 + th];
            end
            D.text = T;
            D.fs   = fs;
            D.fits = Wn <= W + 1e-6 && Hn <= H + 1e-6;
            D.W    = max(W, Wn);
            D.H    = max(H, Hn);
            D.col  = colOf;
            D.row  = ysp;
        end

        function B = tdBoxExit(~, P, Q, R, pad)
            %TDBOXEXIT  Point where the segment P -> Q leaves the box R
            %   = [x y w h], grown by pad on every side. P lies inside.
            x0 = R(1) - pad;  x1 = R(1) + R(3) + pad;
            y0 = R(2) - pad;  y1 = R(2) + R(4) + pad;
            d = Q - P;
            t = 1;
            if d(1) > 0
                t = min(t, (x1 - P(1)) / d(1));
            elseif d(1) < 0
                t = min(t, (x0 - P(1)) / d(1));
            end
            if d(2) > 0
                t = min(t, (y1 - P(2)) / d(2));
            elseif d(2) < 0
                t = min(t, (y0 - P(2)) / d(2));
            end
            B = P + max(t, 0) * d;
        end

        function K = targetRateMatrix(~, sch, tauPar)
            %TARGETRATEMATRIX  Rate matrix K of dc/dt = K*c.
            %   tauPar(p) = 1 / k_p is the time constant of rate parameter
            %   p. Step j (compartment a -> coef * compartment b) has the
            %   rate mult_j * k_par(j): it removes rate * c_a from a and
            %   adds coef * rate * c_a to b (nothing for the sink).
            K = zeros(sch.nS);
            for jj = 1:numel(sch.steps)
                kk = sch.steps(jj).mult / tauPar(sch.steps(jj).par);
                a  = sch.steps(jj).from;
                b  = sch.steps(jj).to;
                K(a, a) = K(a, a) - kk;
                if b > 0
                    K(b, a) = K(b, a) + sch.steps(jj).coef * kk;
                end
            end
        end

        function [C, obsTau, modeW] = targetConc(app, t, K, c0, t0, fwhm, irfMode)
            %TARGETCONC  Concentration profiles of the compartment model
            %   dc/dt = K*c,  c(t0) = c0,  optionally convolved with a
            %   Gaussian IRF.
            %
            %   K = V*diag(lam)/V, so c(t) = sum_m V(:,m) * a_m * E_m(t)
            %   with a = V \ c0 and E_m the m-th decay mode: the closed-
            %   form exp (x) Gaussian of expIRFConv (irfMode 'irf'), or
            %   exp(lam_m (t - t0)) for t >= t0 (irfMode 'none'). A zero
            %   eigenvalue (non-decaying compartment) gives the
            %   integrated IRF / a step.
            %
            %   C      : numel(t) x nCompartments
            %   obsTau : lifetimes -1/lam of the decay modes, ascending
            %            (Inf for a non-decaying mode)
            %   modeW  : nCompartments x nModes weights, C = E * modeW.'
            %            (so the decay-associated spectra are SAS * modeW)
            t  = t(:);
            n  = size(K, 1);
            c0 = c0(:);
            sc = max(abs(diag(K)));
            if ~(sc > 0)
                sc = 1;
            end
            % A defective or nearly degenerate K (e.g. two equal rates in
            % a chain) has no usable eigenvector basis: nudge the rates
            % out of each compartment apart by a relative 1e-6 ... 1e-3
            % until it has one. Scaling a whole COLUMN keeps the balance
            % between what leaves a compartment and what arrives
            % elsewhere.
            Kp = K;
            for attempt = 0:4
                if attempt > 0
                    Kp = K .* (1 + 10^(attempt - 7) * (1:n));
                end
                [V, E] = eig(Kp);
                lam = diag(E);
                if rcond(V) > 1e-9, break; end
            end
            if any(abs(imag(lam)) > 1e-9 * sc)
                error(['The scheme gives oscillating populations ' ...
                    '(complex decay rates, e.g. a closed loop of three ' ...
                    'or more compartments). This is not supported.']);
            end
            lam = real(lam);
            V   = real(V);
            if any(lam > 1e-9 * sc)
                error(['The scheme makes the populations grow without ' ...
                    'bound (a loop that contains a multiplying step?).']);
            end
            if rcond(V) < 1e-12
                error('The rate matrix is degenerate; change the time constants slightly.');
            end
            a = V \ c0;

            noIRF = strcmpi(irfMode, 'none');
            B = zeros(numel(t), n);
            tauM = inf(n, 1);
            for m = 1:n
                if lam(m) < -1e-12 * sc
                    tauM(m) = -1 / lam(m);
                end
                if noIRF
                    x = t - t0;
                    b = zeros(size(x));
                    on = x >= 0;
                    if isfinite(tauM(m))
                        b(on) = exp(-x(on) ./ tauM(m));
                    else
                        b(on) = 1;
                    end
                else
                    b = app.expIRFConv(t, tauM(m), t0, fwhm);
                end
                B(:, m) = b(:);
            end
            [tauS, ord] = sort(tauM, 'ascend');
            modeW = V .* a.';
            modeW = modeW(:, ord);
            C = B(:, ord) * modeW.';
            C(~isfinite(C)) = 0;
            obsTau = tauS(:).';
        end

        function [tauOut, t0Out, fwhmOut, SAS, fitMat, info] = ...
                fitTargetAnalysis(app, D, wl, t, sch, tauInit, tauFixed, ...
                t0Init, fwhmInit, t0Fixed, fwhmFixed, irfMode, tRange, ...
                tauLo, tauHi, ca, fwMode)
            %FITTARGETANALYSIS  Target analysis: fit a kinetic scheme.
            %
            %   Model:  D(lambda, t) = sum_s SAS_s(lambda) * c_s(t), with
            %   c(t) the populations of the scheme `sch` (see
            %   parseTargetScheme / targetConc). Variable projection, as
            %   in fitGlobalAnalysis: the nonlinear parameters (the time
            %   constant 1/k of every rate parameter, t0, FWHM) are
            %   optimised with fminsearch and the spectra are solved by
            %   linear least squares at every step, subject to the
            %   spectral constraints of the scheme (zero ranges, equal /
            %   proportional spectra, non-negativity).
            %
            %   D        : M x N data (wavelength x delay)
            %   wl, t    : spectral axis (M) and delays (N)
            %   tauInit  : time constant of every rate parameter (order
            %              of sch.params);  tauFixed: logical, same size
            %   irfMode  : 'irf' (Gaussian IRF) | 'none' (no IRF: t0 and
            %              FWHM are not fitted)
            %   tRange   : optional [tmin tmax] fit window ([] = all)
            %   tauLo/Hi : optional limits per parameter (NaN = none)
            %   ca       : optional coherent-artifact terms (see caOptions):
            %              extra terms at time zero with free spectra, not
            %              subject to the spectral constraints. Their
            %              spectra are returned in info.caSpec
            %   fwMode   : 'auto' (default) | 'plain' | 'limit': what is
            %              done about a fitted FWHM that the delays do
            %              not resolve -- lower limit of the search and
            %              two-step repeat, as in fitGlobalAnalysis
            %              (info.fwhmLo, fwhmAtLimit, irfRetry)
            %
            %   SAS      : M x nCompartments (NaN where not determinable)
            %   info     : rms, initialRMS, iters, converged, nUsed,
            %              tauStd / t0Std / fwhmStd (1 sigma; NaN where
            %              fixed, Inf where not determined), atLimit,
            %              conc (N x nCompartments), K, obsTau, nObs,
            %              DAS (M x nModes, lifetimes obsTau), reordered
            t  = t(:);
            wl = wl(:);
            [M, N] = size(D);
            nS = sch.nS;
            nP = sch.nP;
            tauInit  = tauInit(:).';
            tauFixed = logical(tauFixed(:).');
            if isscalar(tauFixed), tauFixed = repmat(tauFixed, 1, nP); end
            if numel(tauInit) ~= nP || numel(tauFixed) ~= nP
                error('fitTargetAnalysis: one time constant per rate parameter (%d) is needed.', nP);
            end
            if any(~isfinite(tauInit)) || any(tauInit <= 0)
                error('All time constants must be positive and finite.');
            end
            if numel(t) ~= N || numel(wl) ~= M
                error('fitTargetAnalysis: axes do not match the data size.');
            end
            if nargin < 13, tRange = []; end
            if nargin < 14 || isempty(tauLo), tauLo = nan(1, nP); end
            if nargin < 15 || isempty(tauHi), tauHi = nan(1, nP); end
            tauLo = tauLo(:).';  tauHi = tauHi(:).';
            if numel(tauLo) ~= nP || numel(tauHi) ~= nP
                error('fitTargetAnalysis: limits need one entry per rate parameter.');
            end
            noIRF = strcmpi(irfMode, 'none');
            if noIRF
                t0Fixed = true;  fwhmFixed = true;
            elseif ~(fwhmInit > 0)
                error('The IRF FWHM must be positive.');
            end
            if nargin < 16, ca = []; end
            if nargin < 17 || isempty(fwMode), fwMode = 'auto'; end
            ca = app.caOptions(ca);
            if noIRF, ca.order = 0; end
            caK    = ca.order;
            caOwn  = caK > 0 && ~strcmp(ca.mode, 'irf');
            caFree = caK > 0 && strcmp(ca.mode, 'free');
            caW    = ca.fwhm;
            if caOwn && ~(caW > 0), caW = fwhmInit; end

            % ---- delays used in the fit (fixed for the whole run) ----
            tSel = true(1, N);
            if ~isempty(tRange) && numel(tRange) == 2 && all(isfinite(tRange))
                tSel = (t.' >= min(tRange)) & (t.' <= max(tRange));
            end
            if ~any(tSel)
                error('No delay points inside the fit t range.');
            end
            tSpan = max(t(tSel)) - min(t(tSel));
            if ~(tSpan > 0), tSpan = max(abs(t)) + eps; end
            % Smallest FWHM the delays resolve, and the lower limit of
            % the search in 'limit' mode (see fitGlobalAnalysis); the
            % limit is in force during the search only.
            fwRes = 0;  fwLo = 0;  fwWall = false;
            if ~fwhmFixed
                dzIrf = app.irfDelayStep(t(tSel), t0Init);
                if isfinite(dzIrf) && dzIrf > 0
                    fwRes  = 0.5 * dzIrf;
                    fwLo   = min(fwRes, fwhmInit) * (1 - 1e-9);
                    fwWall = strcmp(fwMode, 'limit');
                end
            end
            if caK > 0
                if caOwn, w0 = caW; else, w0 = fwhmInit; end
                nNear = nnz(tSel & abs(t.' - t0Init) <= 1.5 * w0);
                if nNear < caK + 2
                    error(['Coherent artifact: %d fitted delay(s) lie ' ...
                        'within 1.5 FWHM of time zero, %d are needed for ' ...
                        '%d term(s). Use fewer terms or switch them off.'], ...
                        nNear, caK + 2, caK);
                end
            end

            % ---- spectral constraints: one linear map per wavelength ----
            % code(i, s): 0 free, -1 forced to zero, k > 0 tied by
            % relation k. Wavelengths with the same code row share the
            % map Tm (nS x nUnknown): SAS_row = a * Tm.', so the row is
            % fitted with the basis C * Tm.
            code = zeros(M, nS);
            for zz = 1:numel(sch.rel)
                rows = wl >= sch.rel(zz).lo & wl <= sch.rel(zz).hi;
                code(rows, sch.rel(zz).dep) = zz;
            end
            for zz = 1:numel(sch.zero)
                rows = wl >= sch.zero(zz).lo & wl <= sch.zero(zz).hi;
                code(rows, sch.zero(zz).sp) = -1;
            end
            [pat, ~, grp] = unique(code, 'rows');
            nG = size(pat, 1);
            Tg = cell(1, nG);  nnG = cell(1, nG);
            for g = 1:nG
                Tm = eye(nS);
                Tm(pat(g, :) == -1, :) = 0;
                for pass = 1:nS          % resolve chains of relations
                    for s = find(pat(g, :) > 0)
                        rr = sch.rel(pat(g, s));
                        Tm(s, :) = rr.ratio * Tm(rr.master, :);
                    end
                end
                keepC = any(Tm ~= 0, 1);
                Tg{g}  = Tm(:, keepC);
                nnG{g} = sch.nonneg(keepC);
            end
            rowFinite = all(isfinite(D(:, tSel)), 2);

            % ---- nonlinear parameters ----
            % Free time constants are searched in log space. One with
            % user limits is mapped through a sine so that the search
            % can never leave [lo, hi]; the others only have the wide
            % sanity range below.
            idxFree = find(~tauFixed);
            nFree   = numel(idxFree);
            gLo = 1e-7 * tSpan;  gHi = 1e3 * tSpan;
            hasLim = false(1, nFree);
            lLo = zeros(1, nFree);  lHi = zeros(1, nFree);
            for q = 1:nFree
                lo = tauLo(idxFree(q));  hi = tauHi(idxFree(q));
                hasLim(q) = (isfinite(lo) && lo > 0) || (isfinite(hi) && hi > 0);
                if ~(isfinite(lo) && lo > 0), lo = gLo; end
                if ~(isfinite(hi) && hi > 0), hi = gHi; end
                if ~(hi > lo)
                    error('Parameter %s: the upper limit must exceed the lower limit.', ...
                        sch.params{idxFree(q)});
                end
                lLo(q) = log(lo);  lHi(q) = log(hi);
                if hasLim(q) && (tauInit(idxFree(q)) < lo || tauInit(idxFree(q)) > hi)
                    error('Parameter %s: the starting value %g is outside its limits [%g, %g].', ...
                        sch.params{idxFree(q)}, tauInit(idxFree(q)), lo, hi);
                end
            end
            x0 = tau2x(tauInit(idxFree)).';
            if ~t0Fixed,   x0 = [x0; t0Init];        end
            if ~fwhmFixed, x0 = [x0; log(fwhmInit)]; end
            if caFree,     x0 = [x0; log(caW)];      end

            tauCur = tauInit;  t0Cur = t0Init;  fwhmCur = fwhmInit;
            Ccur = zeros(N, nS);
            CcaCur = zeros(N, caK);    % artifact terms and their spectra
            AcaCur = zeros(M, caK);    % at the last objective call
            nLoss = 0;  nLinCur = 0;
            rowDone = false(M, 1);   % rows whose spectra were solved
            reordered = false;

            function xv = tau2x(tv)
                % physical time constants of the free parameters -> x
                xv = zeros(1, nFree);
                for q1 = 1:nFree
                    lt = log(max(tv(q1), realmin));
                    if hasLim(q1)
                        fr = (lt - lLo(q1)) / (lHi(q1) - lLo(q1));
                        fr = min(max(fr, 0), 1);
                        xv(q1) = asin(2 * fr - 1);
                    else
                        xv(q1) = lt;
                    end
                end
            end

            function unpack(x)
                tauCur = tauInit;
                for q2 = 1:nFree
                    if hasLim(q2)
                        lt = lLo(q2) + (lHi(q2) - lLo(q2)) * (sin(x(q2)) + 1) / 2;
                    else
                        lt = x(q2);
                    end
                    tauCur(idxFree(q2)) = exp(lt);
                end
                p = nFree;
                if ~t0Fixed,   p = p + 1;  t0Cur   = x(p);       end
                if ~fwhmFixed, p = p + 1;  fwhmCur = exp(x(p));  end
                if caFree,     p = p + 1;  caW     = exp(x(p));  end
            end

            function [loss, A, fitM] = objective(x)
                unpack(x);
                A = zeros(M, nS);  fitM = zeros(M, N);  loss = 1e30;
                AcaCur = zeros(M, caK);
                nLinCur = 0;
                rowDone(:) = false;
                % Keep the search inside a sane region.
                if nFree > 0
                    tf = tauCur(idxFree);
                    if any(~isfinite(tf)) || any(tf > gHi * (1 + 1e-9)) || ...
                            any(tf < gLo * (1 - 1e-9))
                        return;
                    end
                end
                if ~noIRF
                    if ~isfinite(fwhmCur) || fwhmCur <= 0 || ...
                            fwhmCur > tSpan || ~isfinite(t0Cur) || ...
                            t0Cur < min(t) - tSpan || t0Cur > max(t) || ...
                            (fwWall && fwhmCur < fwLo)
                        return;
                    end
                    % a fitted artifact width stays within a factor of
                    % 5 of the IRF width
                    if caFree && (caW < fwhmCur / 5 || caW > 5 * fwhmCur)
                        return;
                    end
                end
                try
                    Kc = app.targetRateMatrix(sch, tauCur);
                    Ccur = app.targetConc(t, Kc, sch.c0, t0Cur, fwhmCur, irfMode);
                catch
                    return;
                end
                if caK > 0
                    if caOwn, wNow = caW; else, wNow = fwhmCur; end
                    CcaCur = app.caBasis(t, t0Cur, wNow, caK);
                end
                for gg = 1:nG
                    Tm2 = Tg{gg};
                    nU  = size(Tm2, 2);
                    if nU == 0 && caK == 0
                        % every spectrum is constrained to zero here
                        rowDone(grp == gg) = true;
                        continue;
                    end
                    % Species (with the spectral constraints) followed by
                    % the artifact terms, which are never constrained.
                    Bg  = [Ccur * Tm2, CcaCur];       % N x (nU + caK)
                    nnc = [reshape(nnG{gg}, 1, []), false(1, caK)];
                    rws = find(grp == gg);
                    fin = rowFinite(rws);
                    if any(fin)
                        r1 = rws(fin);
                        Xa = lsSolve(Bg(tSel, :), D(r1, tSel).').';   % rows x nU
                        if any(nnc)
                            % Only the rows whose free solution is
                            % negative need the constrained solve.
                            tolN = 1e-10 * max(abs(Xa(:)));
                            bad = find(any(Xa(:, nnc) < -tolN, 2)).';
                            for ib = bad
                                Xa(ib, :) = nnSolve(Bg(tSel, :), ...
                                    D(r1(ib), tSel).', nnc);
                            end
                            Xa(:, nnc) = max(Xa(:, nnc), 0);
                        end
                        A(r1, :) = Xa(:, 1:nU) * Tm2.';
                        AcaCur(r1, :) = Xa(:, nU+1:end);
                        rowDone(r1) = true;
                        nLinCur = nLinCur + nnz(Xa);
                    end
                    for r2 = rws(~fin).'
                        mm = tSel & isfinite(D(r2, :));
                        if nnz(mm) <= nU + caK, continue; end
                        xa = lsSolve(Bg(mm, :), D(r2, mm).').';
                        if any(nnc) && any(xa(nnc) < -1e-10 * max(abs(xa)))
                            xa = nnSolve(Bg(mm, :), D(r2, mm).', nnc);
                        end
                        if any(nnc), xa(nnc) = max(xa(nnc), 0); end
                        A(r2, :) = xa(1:nU) * Tm2.';
                        AcaCur(r2, :) = xa(nU+1:end);
                        rowDone(r2) = true;
                        nLinCur = nLinCur + nnz(xa);
                    end
                end
                if any(~isfinite(A(:))) || any(~isfinite(AcaCur(:)))
                    A(:) = 0;  AcaCur(:) = 0;  return;
                end
                fitM = A * Ccur.' + AcaCur * CcaCur.';
                R = D - fitM;
                mfin = isfinite(R) & repmat(tSel, M, 1);
                nLoss = nnz(mfin);
                if nLoss > 0
                    loss = sum(R(mfin).^2);
                end
            end

            function X = lsSolve(Cm, Y)
                try
                    X = lsqminnorm(Cm, Y);
                catch
                    X = pinv(Cm) * Y;
                end
            end

            function xa = nnSolve(Bm, y, nnc)
                % Least squares with xa(nnc) >= 0 and the other
                % amplitudes free (written as a difference of two
                % non-negative parts).
                if all(nnc)
                    xa = lsqnonneg(Bm, y).';
                else
                    n1 = nnz(nnc);  n2 = nnz(~nnc);
                    zz2 = lsqnonneg([Bm(:, nnc), Bm(:, ~nnc), -Bm(:, ~nnc)], y);
                    xa = zeros(1, numel(nnc));
                    xa(nnc)  = zz2(1:n1);
                    xa(~nnc) = zz2(n1+1:n1+n2) - zz2(n1+n2+1:end);
                end
            end

            % Evaluate the model once outside the objective, so that a
            % scheme that cannot be computed reports WHY (the objective
            % only returns a penalty).
            K0 = app.targetRateMatrix(sch, tauInit);
            app.targetConc(t, K0, sch.c0, t0Init, fwhmInit, irfMode);

            if isempty(x0)
                [loss, A, fitMat] = objective([]);
                iters = 0;  converged = true;
                initLoss = loss;  nInit = nLoss;  nOpt = nLoss;
                xOpt = [];
            else
                [initLoss, ~, ~] = objective(x0);
                nInit = nLoss;
                if initLoss >= 1e30
                    error(['The starting values are outside the allowed ' ...
                        'range (time constants between 1e-7 and 1e3 ' ...
                        'times the fitted delay span; t0 inside the ' ...
                        'data; FWHM below the span).']);
                end
                opts = optimset('Display', 'off', ...
                    'TolX', 1e-8, 'TolFun', 1e-10, ...
                    'MaxIter', 5000, 'MaxFunEvals', 20000);
                [xOpt, loss, flag, out] = fminsearch(@objective, x0, opts);
                iters = out.iterations;
                converged = flag > 0;
                % Equivalent solutions: without spectral constraints,
                % exchanging the time constants of two steps can fit
                % the data exactly as well (a chain A -> B -> C with
                % the fast and the slow step interchanged has the same
                % decay rates). The data cannot tell these apart, so
                % keep the ORDER of the starting values.
                swp = find(~hasLim);
                if numel(swp) > 1
                    [~, rkInit] = sort(tauInit(idxFree(swp)));
                    tFit = exp(xOpt(swp)).';
                    tAlt = tFit;
                    tAlt(rkInit) = sort(tFit);
                    if any(abs(tAlt - tFit) > 1e-9 * tFit)
                        xAlt = xOpt;
                        xAlt(swp) = log(tAlt(:));
                        if objective(xAlt) <= loss * (1 + 1e-3)
                            [xAlt, lossAlt] = fminsearch(@objective, xAlt, opts);
                            if lossAlt <= loss * (1 + 1e-6)
                                xOpt = xAlt;  loss = lossAlt;
                                reordered = true;
                            end
                        end
                    end
                end
                [~, A, fitMat] = objective(xOpt);
                nOpt  = nLoss;
            end
            nLinOpt = nLinCur;
            tauOut = tauCur;  t0Out = t0Cur;  fwhmOut = fwhmCur;
            Copt = Ccur;

            % A free parameter that ended on one of its limits is not a
            % fitted value: say so and leave it out of the error
            % estimate.
            atLimit = false(1, nP);
            for q = 1:nFree
                if hasLim(q)
                    fr0 = (log(tauOut(idxFree(q))) - lLo(q)) / (lHi(q) - lLo(q));
                    atLimit(idxFree(q)) = fr0 < 1e-4 || fr0 > 1 - 1e-4;
                end
            end

            fwWall = false;      % the search is over

            % ---- 1-sigma uncertainties (same machinery as GA) ----
            tauStd = nan(1, nP);  t0Std = NaN;  fwhmStd = NaN;  caStd = NaN;
            if ~isempty(x0)
                useU = ~atLimit(idxFree);          % free and not on a limit
                pOpt = tauOut(idxFree(useU));
                nTau = numel(pOpt);
                if ~t0Fixed,   pOpt = [pOpt, t0Out];   end
                if ~fwhmFixed, pOpt = [pOpt, fwhmOut]; end
                if caFree,     pOpt = [pOpt, caW];     end
                if ~isempty(pOpt)
                    sig = app.fitParamUncertainty(@(p) residualAtPhys(p), ...
                        pOpt(:), loss, nOpt - nLinOpt);
                    fU = idxFree(useU);
                    tauStd(fU) = sig(1:nTau);
                    p2 = nTau;
                    if ~t0Fixed,   p2 = p2 + 1;  t0Std   = sig(p2);  end
                    if ~fwhmFixed, p2 = p2 + 1;  fwhmStd = sig(p2);  end
                    if caFree,     p2 = p2 + 1;  caStd   = sig(p2);  end
                end
                % Leave the shared state at the optimum again.
                [~, A, fitMat] = objective(xOpt);
            end

            function r = residualAtPhys(pp)
                tv = tauOut(idxFree);
                tv(useU) = pp(1:nTau).';
                xx = zeros(numel(x0), 1);
                xx(1:nFree) = tau2x(tv);
                p3 = nTau;  p4 = nFree;
                if ~t0Fixed
                    p3 = p3 + 1;  p4 = p4 + 1;  xx(p4) = pp(p3);
                end
                if ~fwhmFixed
                    p3 = p3 + 1;  p4 = p4 + 1;
                    xx(p4) = log(max(pp(p3), realmin));
                end
                if caFree
                    p3 = p3 + 1;  p4 = p4 + 1;
                    xx(p4) = log(max(pp(p3), realmin));
                end
                [~, ~, fitM_] = objective(xx);
                R_ = D - fitM_;
                mf = isfinite(D) & repmat(tSel, M, 1);
                r = R_(mf);
                r(~isfinite(r)) = 0;
            end

            % ---- outputs ----
            SAS = A;
            % Rows that could not be solved (masked wavelengths). A row
            % solved as all zeros -- possible with non-negative spectra
            % -- is a result, not a gap.
            SAS(~rowDone, :) = NaN;
            caSpec = AcaCur;
            caSpec(~rowDone, :) = NaN;
            caWout = NaN;
            if caK > 0
                if caOwn, caWout = caW; else, caWout = fwhmOut; end
            end
            Kopt = app.targetRateMatrix(sch, tauOut);
            [~, obsTau, modeW] = app.targetConc(t, Kopt, sch.c0, ...
                t0Out, fwhmOut, irfMode);
            finTau = obsTau(isfinite(obsTau));
            nObs = 0;
            if ~isempty(finTau)
                nObs = 1 + nnz(diff(finTau) > 1e-6 * finTau(1:end-1));
            end
            info = struct('rss', loss, 'iters', iters, ...
                'converged', converged, ...
                'rms', sqrt(loss / max(1, nOpt)), ...
                'initialRMS', sqrt(initLoss / max(1, nInit)), ...
                'nUsed', nOpt, ...
                'tauStd', tauStd, 't0Std', t0Std, 'fwhmStd', fwhmStd, ...
                'atLimit', atLimit, ...
                'conc', Copt, 'K', Kopt, 'obsTau', obsTau, 'nObs', nObs, ...
                'DAS', SAS * modeW, ...
                'nFree', nFree, 'tSel', tSel, 'reordered', reordered, ...
                'caOrder', caK, 'caSpec', caSpec, 'caFwhm', caWout, ...
                'caFwhmStd', caStd, 'caMode', ca.mode, ...
                'fwhmLo', NaN, 'fwhmAtLimit', false, 'irfRetry', 0);
            if ~fwhmFixed && fwRes > 0
                info.fwhmLo      = fwRes;
                info.fwhmAtLimit = app.irfUnresolved(fwhmOut, fwRes);
            end

            % ---- two-step repeat when the IRF width was lost ----
            if ~fwhmFixed && strcmp(fwMode, 'auto') && ~isempty(x0) && ...
                    app.irfWidthLost(fwhmOut, fwhmStd, fwRes)
                % A start value that the delays do not resolve gives
                % the search nothing to go by: two delay steps instead.
                fwS0 = fwhmInit;
                if app.irfUnresolved(fwhmInit, fwRes), fwS0 = 4 * fwRes; end
                try
                    [tauS1, ~, ~, ~, ~, infoS1] = app.fitTargetAnalysis( ...
                        D, wl, t, sch, tauInit, tauFixed, t0Init, ...
                        fwS0, true, true, irfMode, tRange, tauLo, ...
                        tauHi, ca, 'plain');
                    % a start value must lie inside its limits
                    qLo = ~tauFixed & isfinite(tauLo) & tauLo > 0;
                    qHi = ~tauFixed & isfinite(tauHi) & tauHi > 0;
                    tauS1(qLo) = max(tauS1(qLo), tauLo(qLo));
                    tauS1(qHi) = min(tauS1(qHi), tauHi(qHi));
                    [tauS2, t0S2, fwS2, SASS2, fitS2, infoS2] = ...
                        app.fitTargetAnalysis(D, wl, t, sch, tauS1, ...
                        tauFixed, t0Init, fwS0, t0Fixed, fwhmFixed, ...
                        irfMode, tRange, tauLo, tauHi, ca, 'limit');
                    itAll = info.iters + infoS1.iters + infoS2.iters;
                    if infoS2.rss < info.rss * (1 - 1e-6) || ...
                            (info.fwhmAtLimit && ...
                             infoS2.rss <= info.rss * (1 + 1e-6))
                        infoS2.initialRMS = info.initialRMS;
                        tauOut = tauS2;  t0Out = t0S2;  fwhmOut = fwS2;
                        SAS = SASS2;  fitMat = fitS2;
                        info = infoS2;
                        info.irfRetry = 1;
                    else
                        info.irfRetry = 2;
                    end
                    info.iters = itAll;
                catch
                    % the repeat failed: the first result stands
                end
            end
        end
        
        function [Amap, Drec, info] = computeLDA(app, D, t, tauGrid, ...
                t0, fwhm, alpha, regType, irfMode, tRange)
            %COMPUTELDA  Lifetime Density Analysis.
            %   Decompose D(lambda, t) as
            %       D(:, j) = sum_k  A(:, k) * exp(-(t_j - t0)/tau_k)  (IRF-conv.)
            %   for a fixed log-spaced tau grid, with Tikhonov
            %   regularization to control ill-posedness.
            %
            %   Inputs:
            %     D        : M x N   (wavelength x delay)
            %     t        : 1 x N   delay axis
            %     tauGrid  : 1 x K   lifetime values (log-spaced typical)
            %     t0,fwhm  : Gaussian-IRF parameters
            %     alpha    : regularization strength (>= 0)
            %     regType  : 'l2'      -> ||A||^2  penalty
            %                'l2deriv' -> ||L*A||^2 (1st-derivative on tau)
            %
            %   Outputs:
            %     Amap : M x K  amplitude density (signed, real)
            %     Drec : M x N  reconstruction
            %     info : struct with rms, ratio, tauGrid (echo)
            
            t       = t(:).';
            tauGrid = tauGrid(:).';
            K       = numel(tauGrid);
            N       = numel(t);
            M       = size(D, 1);
            
            if nargin < 9 || isempty(irfMode), irfMode = 'irf'; end
            % Optional fit-time-range restriction: delay points outside
            % [tRange(1), tRange(2)] are excluded from the LS solve;
            % the reconstruction is still evaluated on the full axis.
            if nargin < 10 || isempty(tRange) || ~all(isfinite(tRange))
                tSel = true(1, N);
            else
                tSel = (t >= min(tRange)) & (t <= max(tRange));
                if ~any(tSel), tSel = true(1, N); end
            end
            
            % Build basis matrix E (N x K): each column is the decay for
            % one tau on the data delay grid -- IRF-convolved by
            % default, or a pure exponential exp(-(t-t0)/tau).H(t-t0)
            % when irfMode is 'none'.
            E = zeros(N, K);
            if strcmpi(irfMode, 'none')
                tt  = t(:) - t0;
                stp = double(tt >= 0);
                for k = 1:K
                    E(:, k) = stp .* exp(-max(tt, 0) ./ tauGrid(k));
                end
            else
                for k = 1:K
                    E(:, k) = app.expIRFConv(t, tauGrid(k), t0, fwhm);
                end
            end
            
            % Regularization matrix
            switch lower(regType)
                case 'l2'
                    L = eye(K);
                case 'l2deriv'
                    % First-difference operator on the tau grid (in
                    % log-tau spacing for log-spaced grids -- weight by
                    % spacing so the penalty is uniform in log space).
                    if K >= 2
                        d = diff(log(tauGrid(:)));
                        % Use mean spacing as weight (works for nearly
                        % log-uniform grids, which is the common case).
                        ws = mean(d);
                        L = (-eye(K) + diag(ones(K-1, 1), 1)) / ws;
                        L = L(1:end-1, :);
                    else
                        L = eye(K);
                    end
                otherwise
                    L = eye(K);
            end
            
            % Solve (E'E + alpha^2 L'L) A' = E' D'  ->  A = (D * E_pinv).
            % We do it via the augmented LS system for numerical stability:
            %     [ E      ]            [ D' ]
            %     [ alpha*L]  *  A'  =  [ 0  ]
            % which is what lsqminnorm naturally handles.
            Eaug = [E(tSel, :); alpha * L];
            % Replace NaNs in the data with zeros so the LS solve is
            % well-defined; downstream callers may want to mask before
            % calling LDA, but we don't crash on NaN.
            Dclean = D(:, tSel);
            Dclean(~isfinite(Dclean)) = 0;
            Daug = [Dclean.'; zeros(size(L, 1), M)];
            
            try
                Atrans = lsqminnorm(Eaug, Daug);
            catch
                Atrans = pinv(Eaug) * Daug;
            end
            Amap = Atrans.';   % M x K
            Amap(~isfinite(Amap)) = 0;
            
            % Reconstruction on the FULL axis; residual metrics on the
            % fitted window only (consistent with what was minimized).
            Drec = Amap * E.';
            R    = D(:, tSel) - Drec(:, tSel);
            rms  = sqrt(mean(R(isfinite(R)).^2, 'omitnan'));
            Dw   = D(:, tSel);
            denom = sqrt(mean(Dw(isfinite(Dw)).^2, 'omitnan'));
            if denom > 0, ratio = rms / denom; else, ratio = NaN; end
            
            info = struct('rms', rms, 'ratio', ratio, ...
                'tauGrid', tauGrid, 'alpha', alpha, ...
                'regType', regType, 'K', K, 'M', M, 'N', N);
        end
        
        function [alphaList, resNorm, solNorm, alphaCorner] = ...
                computeLcurve(app, D, t, tauGrid, t0, fwhm, regType, ...
                    alphaList, irfMode, tRange)
            %COMPUTELCURVE  Sweep regularization parameter alpha over a
            %   range and return (residual norm, solution norm) pairs
            %   suitable for plotting an L-curve. The "corner" is
            %   estimated by the maximum-curvature point in log-log
            %   space (a robust simple heuristic).
            
            if nargin < 8 || isempty(alphaList)
                alphaList = logspace(-6, 2, 25);
            end
            alphaList = alphaList(:).';
            
            if nargin < 9,  irfMode = 'irf'; end
            if nargin < 10, tRange  = [];    end
            % Residual norms must be evaluated on the same delay window
            % the LDA solve uses, otherwise the L-curve mixes fitted
            % and unfitted regions and its corner shifts.
            tR = t(:).';
            if isempty(tRange) || ~all(isfinite(tRange))
                tSelL = true(1, numel(tR));
            else
                tSelL = (tR >= min(tRange)) & (tR <= max(tRange));
                if ~any(tSelL), tSelL = true(1, numel(tR)); end
            end
            n = numel(alphaList);
            resNorm = zeros(1, n);
            solNorm = zeros(1, n);
            for ii = 1:n
                [Amap, Drec, ~] = app.computeLDA(D, t, tauGrid, ...
                    t0, fwhm, alphaList(ii), regType, irfMode, tRange);
                R = D(:, tSelL) - Drec(:, tSelL);
                resNorm(ii) = sqrt(sum(R(isfinite(R)).^2));
                solNorm(ii) = sqrt(sum(Amap(isfinite(Amap)).^2));
            end
            
            % Maximum-curvature corner in log-log space
            x = log(resNorm);  y = log(solNorm);
            ok = isfinite(x) & isfinite(y);
            x = x(ok);  y = y(ok);
            alphaListOk = alphaList(ok);
            if numel(x) < 4
                alphaCorner = alphaList(round(n/2));
                return;
            end
            dx = gradient(x);  dy = gradient(y);
            ddx = gradient(dx); ddy = gradient(dy);
            curv = abs(dx .* ddy - dy .* ddx) ./ (dx.^2 + dy.^2).^1.5;
            curv(~isfinite(curv)) = 0;
            [~, iMax] = max(curv);
            alphaCorner = alphaListOk(iMax);
        end
        
        function [C, S, info] = computeMCR(~, D, N, options)
            %COMPUTEMCR  Multivariate Curve Resolution by Alternating
            %   Least Squares (MCR-ALS).
            %
            %   D = C * S' + E
            %     D : (Ndelay x Nwl)  data matrix (caller orients it so
            %         that rows are time, cols are wavelength)
            %     C : (Ndelay x N)    concentration profiles
            %     S : (Nwl    x N)    pure spectra
            %
            %   Inputs:
            %     D       : data matrix (Nt x Nwl). NaNs are treated as
            %               missing in the LS solves.
            %     N       : number of components
            %     options : struct with fields
            %         .nnC        (logical) non-neg constraint on C
            %         .uniC       (logical) unimodality constraint on C
            %         .nnS        (logical) non-neg constraint on S
            %         .tol        (double)  rel. LOF change tolerance
            %         .maxIter    (int)     max ALS iterations
            %         .initMode   ('svd' | 'random' | 'custom')
            %         .seed       'random': seed of the starting profiles
            %                     (positive integer). Empty / 0 = a new
            %                     seed from the clock; info.seed tells
            %                     which one was used, so that the run can
            %                     be repeated.
            %         .Cinit      (Nt x N) custom initial C (if initMode='custom')
            %         .Sinit      (Nwl x k) custom initial S (used only if
            %                     Cinit is empty in custom mode). Fewer
            %                     than N columns are completed from the
            %                     data (see startSpectra); non-finite
            %                     entries count as 0.
            %
            %   Outputs:
            %     C, S  : final concentration profiles and pure spectra
            %     info  : struct with .lof, .iter, .lofHistory, .converged,
            %             .seed (NaN unless 'random') and .nSinit (number
            %             of supplied starting spectra used; NaN unless
            %             they were used)
            
            if nargin < 4 || isempty(options), options = struct(); end
            if ~isfield(options, 'nnC'),     options.nnC = true;   end
            if ~isfield(options, 'uniC'),    options.uniC = false; end
            if ~isfield(options, 'nnS'),     options.nnS = false;  end
            if ~isfield(options, 'tol'),     options.tol = 1e-4;   end
            if ~isfield(options, 'maxIter'), options.maxIter = 200; end
            if ~isfield(options, 'initMode'), options.initMode = 'svd'; end
            if ~isfield(options, 'Cinit'),   options.Cinit = [];   end
            if ~isfield(options, 'Sinit'),   options.Sinit = [];   end
            if ~isfield(options, 'seed'),    options.seed = [];    end
            
            [Nt, Nwl] = size(D);
            
            % Replace NaN with 0 (and track the mask). The clean matrix is
            % used in the LS solves; LOF is computed only over valid cells.
            validMask = isfinite(D);
            Dclean    = D;
            Dclean(~validMask) = 0;
            
            % --- Initial estimate ---
            S = [];
            seedUsed = NaN;     % seed of a 'random' start
            nSinit   = NaN;     % supplied starting spectra that were used
            switch lower(options.initMode)
                case 'custom'
                    if ~isempty(options.Cinit) && size(options.Cinit, 2) == N
                        C = options.Cinit;
                    else
                        [S, nUsed] = startSpectra(Dclean, options.Sinit, N);
                        if isempty(S)
                            C = initSVD(Dclean, N);  % nothing usable given
                        else
                            C = solveC(Dclean, S);   % seed C from the spectra
                            nSinit = nUsed;
                        end
                    end
                case 'random'
                    % Positive random numbers from a random stream of
                    % its own: the run can be repeated from its seed,
                    % and MATLAB's global stream (the user's other
                    % random numbers) is left alone. Without a seed, a
                    % new one is taken from the clock for every run.
                    seedUsed = options.seed;
                    if isempty(seedUsed) || ~isnumeric(seedUsed) || ...
                            ~isscalar(seedUsed) || ~isfinite(seedUsed) || ...
                            seedUsed < 1
                        seedUsed = 1 + mod(floor(now * 8.64e7), 999999); %#ok<TNOW1>
                    end
                    seedUsed = min(floor(double(seedUsed)), 2^32 - 1);
                    rs = RandStream('mt19937ar', 'Seed', seedUsed);
                    C  = abs(randn(rs, Nt, N));
                otherwise   % 'svd'
                    C = initSVD(Dclean, N);
            end
            if isempty(S)
                S = solveS(Dclean, C);
            end
            
            % --- ALS loop ---
            lofHistory = zeros(1, options.maxIter);
            prevLOF = inf;
            converged = false;
            for it = 1:options.maxIter
                % Solve S given C
                S = solveS(Dclean, C);
                S = applyConstraintsS(S, options);
                
                % Solve C given S
                C = solveC(Dclean, S);
                C = applyConstraintsC(C, options);
                
                % LOF
                Drec = C * S.';
                R    = D - Drec;
                Rv   = R(validMask);
                Dv   = D(validMask);
                lof  = 100 * sqrt(sum(Rv.^2) / sum(Dv.^2));
                lofHistory(it) = lof;
                
                if abs(prevLOF - lof) / max(prevLOF, eps) < options.tol
                    converged = true;
                    lofHistory = lofHistory(1:it);
                    break;
                end
                prevLOF = lof;
            end
            if ~converged, lofHistory = lofHistory(1:it); end
            
            info = struct('lof', lof, 'iter', it, ...
                'lofHistory', lofHistory, 'converged', converged, ...
                'options', options, 'seed', seedUsed, 'nSinit', nSinit);
            
            % ============= nested helpers =============
            function C0 = initSVD(D0, n)
                % Starting concentration profiles: the absolute values
                % of the first n temporal (left) singular vectors of
                % the data.
                [U0, ~, ~] = svds(D0, n);
                % Make the dominant entry of each column positive.
                for jj = 1:n
                    [~, idx] = max(abs(U0(:, jj)));
                    if U0(idx, jj) < 0
                        U0(:, jj) = -U0(:, jj);
                    end
                end
                C0 = abs(U0);
                if any(~isfinite(C0(:)))
                    C0 = abs(randn(size(D0, 1), n));
                end
            end
            
            function [S0, nUsed] = startSpectra(D0, Sin, n)
                % Starting spectra of the 'custom' mode: the supplied
                % columns (at most n; non-finite entries, e.g. a masked
                % wavelength band, count as 0; empty columns and
                % columns that repeat earlier ones are dropped). When
                % fewer than n remain, the set is completed with the
                % leading spectral singular vectors of the part of the
                % data that the supplied spectra do not span -- a
                % column of zeros would stay zero through every ALS
                % step. nUsed = supplied spectra kept; S0 = [] if none.
                S0 = [];  nUsed = 0;
                if isempty(Sin) || size(Sin, 1) ~= size(D0, 2), return; end
                Sin(~isfinite(Sin)) = 0;
                Sin = Sin(:, 1:min(n, size(Sin, 2)));
                keep = false(1, size(Sin, 2));
                Q = zeros(size(Sin, 1), 0);     % orthonormal basis of the kept ones
                for jj = 1:size(Sin, 2)
                    v  = Sin(:, jj);
                    nv = norm(v);
                    if nv == 0, continue; end
                    r = v - Q * (Q.' * v);
                    r = r - Q * (Q.' * r);
                    if norm(r) > 1e-6 * nv
                        keep(jj) = true;
                        Q = [Q, r / norm(r)]; %#ok<AGROW>
                    end
                end
                nUsed = nnz(keep);
                if nUsed == 0, return; end
                S0 = Sin(:, keep);
                if nUsed < n
                    Dr = D0 - (D0 * Q) * Q.';
                    [~, ~, V0] = svds(Dr, n - nUsed);
                    if size(V0, 2) < n - nUsed || any(~isfinite(V0(:)))
                        S0 = [];  nUsed = 0;  return;   % -> SVD start
                    end
                    for jj = 1:size(V0, 2)
                        % Sign: the time profile that goes with the added
                        % spectrum has its dominant lobe positive.
                        p = Dr * V0(:, jj);
                        [~, ip] = max(abs(p));
                        if p(ip) < 0, V0(:, jj) = -V0(:, jj); end
                    end
                    % ... and the size of the supplied spectra.
                    S0 = [S0, V0 * mean(sqrt(sum(S0.^2, 1)))];
                end
            end
            
            function Sloc = solveS(Dloc, Cloc)
                % S = (C \ D)'  (least squares)
                try
                    Sloc = (Cloc \ Dloc).';
                catch
                    Sloc = (pinv(Cloc) * Dloc).';
                end
                Sloc(~isfinite(Sloc)) = 0;
            end
            function Cloc = solveC(Dloc, Sloc)
                % C = D / (S')
                try
                    Cloc = Dloc / Sloc.';
                catch
                    Cloc = Dloc * pinv(Sloc.');
                end
                Cloc(~isfinite(Cloc)) = 0;
            end
            function Sloc = applyConstraintsS(Sloc, opt)
                if opt.nnS
                    Sloc(Sloc < 0) = 0;
                end
            end
            function Cloc = applyConstraintsC(Cloc, opt)
                if opt.nnC
                    Cloc(Cloc < 0) = 0;
                end
                if opt.uniC
                    % Force each column to have a single peak: scan from
                    % the global max outward and impose monotonic non-
                    % increase. Cheap "horizontal" unimodality.
                    for jj = 1:size(Cloc, 2)
                        col = Cloc(:, jj);
                        [~, ipk] = max(col);
                        % left of peak: monotonic non-decrease
                        for kk = ipk-1:-1:1
                            if col(kk) > col(kk+1)
                                col(kk) = col(kk+1);
                            end
                        end
                        % right of peak: monotonic non-increase
                        for kk = ipk+1:numel(col)
                            if col(kk) > col(kk-1)
                                col(kk) = col(kk-1);
                            end
                        end
                        Cloc(:, jj) = col;
                    end
                end
            end
        end
        
        function [Pmap, freq, info] = computeCoherence(app, R, t, options)
            %COMPUTECOHERENCE  2D oscillation map (|FFT|^2 of residual).
            %
            %   Inputs:
            %     R       : Nwl x Nt residual matrix (kinetic component
            %               already subtracted; same orientation as
            %               app.deltaA).
            %     t       : 1 x Nt delay vector (in app.timeUnit units).
            %     options : struct with
            %         .tMin, .tMax       (in app.timeUnit units)
            %         .apod              'rect' | 'hann' | 'hamming' | 'blackman'
            %         .zeroPad           integer factor >= 1 (1 = none, 2 = 2x, ...)
            %         .detrend           'none' | 'mean' | 'linear'
            %         .freqUnit          'cm-1' | 'THz' | 'Hz'
            %         .normMode          'peak' | 'none' | 'perWl'
            %
            %   Outputs:
            %     Pmap    : Nwl x Nfreq, |FFT|^2 amplitude map
            %     freq    : 1 x Nfreq frequency axis (in options.freqUnit)
            %     info    : struct (resampled grid, dt, etc.)
            
            if nargin < 4 || isempty(options), options = struct(); end
            if ~isfield(options, 'tMin'),   options.tMin = min(t); end
            if ~isfield(options, 'tMax'),   options.tMax = max(t); end
            if ~isfield(options, 'apod'),   options.apod = 'hann'; end
            if ~isfield(options, 'zeroPad'),options.zeroPad = 2;   end
            if ~isfield(options, 'detrend'),options.detrend = 'linear'; end
            if ~isfield(options, 'freqUnit'),options.freqUnit = 'cm-1'; end
            if ~isfield(options, 'normMode'),options.normMode = 'peak';end
            
            t = t(:).';
            
            % --- Restrict to time window ---
            mask = (t >= options.tMin) & (t <= options.tMax);
            tw   = t(mask);
            Rw   = R(:, mask);
            if numel(tw) < 8
                error('Time window is too short (need >= 8 points, got %d).', ...
                    numel(tw));
            end
            
            % --- Resample to a uniform grid (FFT requires uniform dt) ---
            % Use the median dt of the masked window as the uniform step.
            dts   = diff(tw);
            dtUni = median(dts);
            % The uniform grid uses the MEDIAN step of the window: a long
            % window over a delay axis that is dense early and sparse
            % late would explode into millions of points.
            nUni = abs(tw(end) - tw(1)) / max(abs(dtUni), realmin) + 1;
            if nUni * max(1, round(options.zeroPad)) * size(Rw, 1) > 4e7
                error(['The time window resamples to about %.3g ' ...
                    'uniform points (step %.3g) for %d rows -- too ' ...
                    'large for the FFT map. Narrow the window to the ' ...
                    'densely sampled part of the delay axis.'], ...
                    nUni, dtUni, size(Rw, 1));
            end
            tUni  = tw(1):dtUni:tw(end);
            if numel(tUni) < 8
                error('Resampled grid too short (n=%d).', numel(tUni));
            end
            Nt    = numel(tUni);
            
            % Resample each row.
            Ru = zeros(size(Rw, 1), Nt);
            for ii = 1:size(Rw, 1)
                row = Rw(ii, :);
                m   = isfinite(row);
                if nnz(m) < 4
                    Ru(ii, :) = 0;
                else
                    Ru(ii, :) = interp1(tw(m), row(m), tUni, ...
                        'linear', 0);
                end
            end
            
            % --- Detrend (per row) ---
            switch lower(options.detrend)
                case 'none'
                    % no-op
                case 'mean'
                    Ru = Ru - mean(Ru, 2, 'omitnan');
                case 'linear'
                    % Fit and remove a linear trend per row.
                    A = [tUni(:), ones(numel(tUni), 1)];
                    for ii = 1:size(Ru, 1)
                        y = Ru(ii, :).';
                        m = isfinite(y);
                        if nnz(m) < 2, continue; end
                        c = A(m, :) \ y(m);
                        Ru(ii, :) = y.' - (A * c).';
                    end
            end
            Ru(~isfinite(Ru)) = 0;
            
            % --- Apodization (window function) ---
            n = numel(tUni);
            switch lower(options.apod)
                case 'rect',     w = ones(1, n);
                case 'hann',     w = 0.5 * (1 - cos(2*pi*(0:n-1)/(n-1)));
                case 'hamming',  w = 0.54 - 0.46*cos(2*pi*(0:n-1)/(n-1));
                case 'blackman'
                    k = (0:n-1)/(n-1);
                    w = 0.42 - 0.5*cos(2*pi*k) + 0.08*cos(4*pi*k);
                otherwise,       w = ones(1, n);
            end
            Ru = Ru .* w;
            
            % --- Zero-pad ---
            zpad = max(1, round(options.zeroPad));
            Nfft = max(2^nextpow2(Nt * zpad), Nt);
            
            % --- FFT ---
            X = fft(Ru, Nfft, 2);
            % One-sided power spectrum and phase map.
            Nhalf = floor(Nfft/2) + 1;
            Xhalf = X(:, 1:Nhalf);
            P     = abs(Xhalf).^2;
            Phi   = atan2(imag(Xhalf), real(Xhalf));   % raw [-pi, pi]
            
            % Frequency axis. The FFT bin spacing in 1/timeUnit is
            % 1/(Nfft*dtUni). Convert to the requested unit.
            fNyq    = 1 / (2 * dtUni);          % cycles per timeUnit
            f_perTU = linspace(0, fNyq, Nhalf);  % Nhalf points
            
            % timeUnit -> seconds factor
            switch app.timeUnit
                case 'fs', tu2s = 1e-15;
                case 'ps', tu2s = 1e-12;
                case 'ns', tu2s = 1e-9;
                case 'us', tu2s = 1e-6;
                otherwise, tu2s = 1e-12;
            end
            f_Hz = f_perTU / tu2s;     % cycles per second (Hz)
            
            switch lower(options.freqUnit)
                case 'cm-1'
                    cLightCm = 2.99792458e10;     % cm/s
                    freq = f_Hz / cLightCm;       % cm^-1
                case 'thz'
                    freq = f_Hz / 1e12;
                otherwise
                    freq = f_Hz;
            end
            
            % --- Normalization ---
            switch lower(options.normMode)
                case 'peak'
                    pk = max(P(:), [], 'omitnan');
                    if isfinite(pk) && pk > 0, P = P / pk; end
                case 'perwl'
                    pkW = max(P, [], 2, 'omitnan');
                    pkW(pkW <= 0) = 1;
                    P = P ./ pkW;
                case 'none'
                    % keep raw
            end
            P(~isfinite(P)) = 0;
            
            Pmap = P;
            info = struct('dtUni', dtUni, 'Nfft', Nfft, ...
                'tUni', tUni, 'Nwin', Nt, ...
                'apod', options.apod, 'zeroPad', zpad, ...
                'detrend', options.detrend, ...
                'freqUnit', options.freqUnit, ...
                'normMode', options.normMode, ...
                'Phimap', Phi);    % raw phase, Nwl x Nfreq, [-pi, pi]
        end
        
        function [Ptf, freq, tc, info] = computeTimeFrequency(app, y, t, options)
            %COMPUTETIMEFREQUENCY  Time-frequency map of a single trace.
            %
            %   Two methods, both FFT-based (no Wavelet Toolbox needed):
            %     'stft' : sliding-window FFT. Fixed window -> uniform
            %              time / frequency resolution everywhere.
            %     'cwt'  : Morlet continuous wavelet transform. Window
            %              scales with 1/f -> good frequency resolution at
            %              low wavenumber, good time resolution at high.
            %
            %   Inputs:
            %     y       : 1 x Nt trace (residual at one lambda, or a
            %               lambda-averaged residual)
            %     t       : 1 x Nt delay vector (app.timeUnit units)
            %     options : struct
            %         .tMin, .tMax   analysis window (timeUnit)
            %         .detrend       'none' | 'mean' | 'linear'
            %         .method        'stft' | 'cwt'
            %         .freqUnit      'cm-1' | 'THz' | 'Hz'
            %         .normMode      'peak' | 'none'
            %         STFT: .winLen (timeUnit), .hop (timeUnit),
            %               .apod ('hann'|'hamming'|'blackman'|'rect'),
            %               .zeroPad (integer >= 1)
            %         CWT : .omega0 (Morlet parameter, typ. 5..12),
            %               .nFreq (frequency grid points),
            %               .fLo, .fHi (freqUnit; fHi <= 0 -> Nyquist)
            %
            %   Outputs:
            %     Ptf  : Nfreq x Ntc power map (|X|^2)
            %     freq : 1 x Nfreq frequency axis (freqUnit)
            %     tc   : 1 x Ntc time axis of the map (timeUnit)
            %     info : struct -- tUni, yu (the analysed trace on the
            %            uniform grid), dtUni, method parameters, and a
            %            cone-of-influence descriptor:
            %              .coiFreq : 1 x Ntc, lowest reliable frequency
            %                         at each tc (freqUnit); Inf marks a
            %                         column that is entirely edge-
            %                         affected (STFT half-window).

            if nargin < 4 || isempty(options), options = struct(); end
            def = struct('tMin', min(t), 'tMax', max(t), ...
                'detrend', 'linear', 'method', 'stft', ...
                'freqUnit', 'cm-1', 'normMode', 'peak', ...
                'winLen', NaN, 'hop', NaN, 'apod', 'hann', 'zeroPad', 4, ...
                'omega0', 6, 'nFreq', 200, 'fLo', 0, 'fHi', 0);
            fn = fieldnames(def);
            for k = 1:numel(fn)
                if ~isfield(options, fn{k}), options.(fn{k}) = def.(fn{k}); end
            end

            t = t(:).';  y = y(:).';

            % --- Restrict to the analysis window, resample uniformly ---
            mask = (t >= options.tMin) & (t <= options.tMax) & isfinite(y);
            tw = t(mask);  yw = y(mask);
            if numel(tw) < 16
                error('Time window is too short (need >= 16 points, got %d).', ...
                    numel(tw));
            end
            dtUni = median(diff(tw));
            nUni  = abs(tw(end) - tw(1)) / max(abs(dtUni), realmin) + 1;
            if nUni > 2e5
                error(['The time window resamples to about %.3g ' ...
                    'uniform points (step %.3g) -- too many for a ' ...
                    'time-frequency map. Narrow the window to the ' ...
                    'densely sampled part of the delay axis.'], ...
                    nUni, dtUni);
            end
            tUni  = tw(1):dtUni:tw(end);
            Nt    = numel(tUni);
            if Nt < 16
                error('Resampled grid too short (n=%d).', Nt);
            end
            yu = interp1(tw, yw, tUni, 'linear', 0);

            % --- Detrend ---
            switch lower(options.detrend)
                case 'mean'
                    yu = yu - mean(yu, 'omitnan');
                case 'linear'
                    A = [tUni(:), ones(Nt, 1)];
                    c = A \ yu(:);
                    yu = yu - (A * c).';
            end
            yu(~isfinite(yu)) = 0;

            % --- timeUnit -> seconds, and cycles/timeUnit -> freqUnit ---
            switch app.timeUnit
                case 'fs', tu2s = 1e-15;
                case 'ps', tu2s = 1e-12;
                case 'ns', tu2s = 1e-9;
                case 'us', tu2s = 1e-6;
                otherwise, tu2s = 1e-12;
            end
            switch lower(options.freqUnit)
                case 'cm-1', fConv = 1 / (tu2s * 2.99792458e10);
                case 'thz',  fConv = 1 / (tu2s * 1e12);
                otherwise,   fConv = 1 / tu2s;
            end
            fNyqTU = 1 / (2 * dtUni);            % cycles per timeUnit

            info = struct('dtUni', dtUni, 'tUni', tUni, 'yu', yu, ...
                'method', lower(options.method), ...
                'freqUnit', options.freqUnit, 'detrend', options.detrend);

            switch lower(options.method)
                % =========================================================
                case 'stft'
                    winLen = options.winLen;
                    if ~isfinite(winLen) || winLen <= 0
                        winLen = (tUni(end) - tUni(1)) / 4;
                    end
                    nWin = round(winLen / dtUni);
                    nWin = max(8, min(nWin, Nt));
                    if mod(nWin, 2) == 0, nWin = nWin - 1; end   % odd -> centred
                    halfW = (nWin - 1) / 2;
                    hop = options.hop;
                    if ~isfinite(hop) || hop <= 0, hopN = 1;
                    else, hopN = max(1, round(hop / dtUni)); end

                    switch lower(options.apod)
                        case 'rect',    w = ones(1, nWin);
                        case 'hamming', w = 0.54 - 0.46*cos(2*pi*(0:nWin-1)/(nWin-1));
                        case 'blackman'
                            kk = (0:nWin-1)/(nWin-1);
                            w = 0.42 - 0.5*cos(2*pi*kk) + 0.08*cos(4*pi*kk);
                        otherwise,      w = 0.5 * (1 - cos(2*pi*(0:nWin-1)/(nWin-1)));
                    end
                    zpad  = max(1, round(options.zeroPad));
                    Nfft  = max(2^nextpow2(nWin * zpad), nWin);
                    Nhalf = floor(Nfft/2) + 1;

                    centers = 1:hopN:Nt;
                    Ntc  = numel(centers);
                    ypad = [zeros(1, halfW), yu, zeros(1, halfW)];
                    Ptf  = zeros(Nhalf, Ntc);
                    for k = 1:Ntc
                        c   = centers(k);
                        seg = ypad(c : c + nWin - 1) .* w;
                        X   = fft(seg, Nfft);
                        Ptf(:, k) = abs(X(1:Nhalf)).^2;
                    end
                    tc   = tUni(centers);
                    freq = linspace(0, fNyqTU, Nhalf) * fConv;

                    % Cone of influence: columns whose window extends past
                    % either edge of the data are zero-padded -> unreliable.
                    edgeT = halfW * dtUni;
                    coiFreq = zeros(1, Ntc);
                    coiFreq(tc < tUni(1) + edgeT | tc > tUni(end) - edgeT) = Inf;

                    info.nWin = nWin;  info.winLen = nWin * dtUni;
                    info.hopN = hopN;  info.Nfft = Nfft;
                    info.apod = options.apod;  info.zeroPad = zpad;
                    info.coiFreq = coiFreq;
                    info.coiT = [tUni(1) + edgeT, tUni(end) - edgeT];

                % =========================================================
                case 'cwt'
                    om0 = options.omega0;
                    if ~isfinite(om0) || om0 < 2, om0 = 6; end
                    nF = max(16, round(options.nFreq));
                    % Frequency grid (linear in freqUnit, as on the FFT map)
                    fHi = options.fHi;
                    if ~isfinite(fHi) || fHi <= 0, fHi = fNyqTU * fConv; end
                    fHi = min(fHi, fNyqTU * fConv);
                    fLo = options.fLo;
                    fMinUseful = 2 / (tUni(end) - tUni(1)) * fConv;  % ~2 cycles in window
                    if ~isfinite(fLo) || fLo <= 0, fLo = fMinUseful; end
                    fLo = max(fLo, fMinUseful / 4);
                    if fLo >= fHi, fLo = fHi / 50; end
                    freq  = linspace(fLo, fHi, nF);
                    fTU   = freq / fConv;                 % cycles / timeUnit
                    % Scale <-> frequency via the Morlet CENTRE frequency
                    % f_c = om0 / (2*pi*s). (The Fourier-period relation
                    % used in geophysics biases peak positions by ~1% for
                    % om0 = 6; the centre-frequency mapping puts a pure
                    % sinusoid's response maximum exactly at its f.)
                    sVec  = om0 ./ (2*pi*fTU);              % scales
                    
                    % FFT of the (zero-padded) trace; analytic Morlet in the
                    % frequency domain. PEAK-RESPONSE normalisation: the
                    % wavelet is scaled so a unit-amplitude cosine at its
                    % centre frequency gives |W| = 1 at every scale, i.e.
                    % equal-amplitude modes have equal power, exactly as
                    % on the |FFT|^2 map (unlike the Torrence & Compo
                    % energy normalisation, where power grows as 1/f).
                    Npad  = 2^nextpow2(2 * Nt);
                    Yhat  = fft(yu, Npad);
                    omega = 2*pi * (0:Npad-1) / (Npad * dtUni);
                    omega(omega > pi/dtUni) = omega(omega > pi/dtUni) - 2*pi/dtUni;
                    pos   = omega > 0;
                    Ptf   = zeros(nF, Nt);
                    for k = 1:nF
                        s = sVec(k);
                        psiHat = zeros(1, Npad);
                        psiHat(pos) = 2 * exp(-0.5 * (s*omega(pos) - om0).^2);
                        W = ifft(Yhat .* conj(psiHat));
                        Ptf(k, :) = abs(W(1:Nt)).^2;
                    end
                    tc = tUni;
                    
                    % Cone of influence: the Gaussian envelope's e-folding
                    % time is sqrt(2)*s; at distance d from the nearest
                    % edge, scales with sqrt(2)*s > d are edge-affected
                    % ->  f < coiFreq(d) is unreliable there.
                    d = min(tc - tUni(1), tUni(end) - tc);
                    d(d <= 0) = dtUni / 2;
                    coiFTU  = om0 ./ (2*pi) * sqrt(2) ./ d;
                    coiFreq = coiFTU * fConv;

                    info.omega0 = om0;  info.nFreq = nF;
                    info.scales = sVec;
                    info.coiFreq = coiFreq;
                    info.coiT = [];
                otherwise
                    error('Unknown time-frequency method "%s".', options.method);
            end

            % --- Normalisation ---
            if strcmpi(options.normMode, 'peak')
                pk = max(Ptf(:), [], 'omitnan');
                if isfinite(pk) && pk > 0, Ptf = Ptf / pk; end
            end
            Ptf(~isfinite(Ptf)) = 0;
            info.normMode = options.normMode;
        end

        function [tauOut, betaOut, t0Out, fwhmOut, Aout, fitVec, ...
                resVec, info] = fitSingleTrace(app, t, y, ...
                tauInit, tauFixed, betaInit, betaFixed, stretchOn, ...
                t0Init, t0Fixed, fwhmInit, fwhmFixed, hasInf, irfMode, ...
                tRange, fwMode)
            %FITSINGLETRACE  Fit a single kinetic trace y(t) with a sum
            %   of (possibly stretched) exponentials convolved with a
            %   Gaussian IRF.  Linear amplitudes A_i are solved per
            %   nonlinear-parameter step (VARPRO style); nonlinear
            %   parameters (tau, beta, t0, FWHM) are optimized via
            %   fminsearch over a logged / unbounded parameterization.
            %
            %   Inputs:
            %     t,y         : column vectors, same length
            %     tauInit     : 1xN vector, time-constants in same units as t
            %     tauFixed    : 1xN logical
            %     betaInit    : 1xN vector, stretching exponents (1 = exp)
            %     betaFixed   : 1xN logical
            %     stretchOn   : 1xN logical, true => use stretched form for
            %                   that component (and beta is a free / fixed
            %                   parameter); false => standard exp.
            %     t0Init,Fixed,fwhmInit,Fixed,hasInf : scalars / logical
            %     irfMode     : 'none' (no IRF convolution) | 'numerical'
            %                   (Gaussian IRF; any other value is treated
            %                   like this one)
            %     tRange      : optional [tmin tmax] fit window ([] = all)
            %     fwMode      : 'auto' (default) | 'plain' | 'limit': what
            %                   is done about a fitted FWHM that the
            %                   delays do not resolve -- lower limit of
            %                   the search and two-step repeat, as in
            %                   fitGlobalAnalysis (info.fwhmLo,
            %                   fwhmAtLimit, irfRetry)
            %
            %   Outputs:
            %     tauOut,betaOut : final values (row vectors, length N)
            %     t0Out,fwhmOut  : scalars
            %     Aout           : amplitudes, length N + (hasInf ? 1 : 0)
            %     fitVec, resVec : column vectors, same length as t
            %     info           : struct with iters, exitFlag, RMS, mask
            %                      (logical of points used in the fit)
            
            if nargin < 15, tRange = []; end
            if nargin < 16 || isempty(fwMode), fwMode = 'auto'; end
            t = t(:);  y = y(:);
            if numel(t) ~= numel(y)
                error('fitSingleTrace: t and y must have the same length.');
            end
            tauInit  = tauInit(:).';
            tauFixed = logical(tauFixed(:).');
            betaInit  = betaInit(:).';
            betaFixed = logical(betaFixed(:).');
            stretchOn = logical(stretchOn(:).');
            % A scalar per-component option applies to every component
            % (e.g. stretchOn = false for a plain multi-exponential).
            % Without this a scalar would be indexed out of range for
            % the second component.
            nCmp = numel(tauInit);
            if isscalar(tauFixed),  tauFixed  = repmat(tauFixed,  1, nCmp); end
            if isscalar(betaInit),  betaInit  = repmat(betaInit,  1, nCmp); end
            if isscalar(betaFixed), betaFixed = repmat(betaFixed, 1, nCmp); end
            if isscalar(stretchOn), stretchOn = repmat(stretchOn, 1, nCmp); end
            if numel(tauFixed) ~= nCmp || numel(betaInit) ~= nCmp || ...
                    numel(betaFixed) ~= nCmp || numel(stretchOn) ~= nCmp
                error(['fitSingleTrace: tauFixed, betaInit, betaFixed ' ...
                       'and stretchOn need one entry per time ' ...
                       'constant (%d).'], nCmp);
            end
            
            % In 'none' (no-IRF) mode t0 is degenerate with the linear
            % amplitude and FWHM does not enter the model; force both
            % fixed (mirrors the GA engine).
            if strcmpi(irfMode, 'none')
                t0Fixed   = true;
                fwhmFixed = true;
            end
            
            % Build mask of valid data points.
            mask = isfinite(y);
            % Optional fit-time-range restriction: points outside
            % [tRange(1), tRange(2)] are excluded from the fit.
            if nargin >= 15 && ~isempty(tRange) && all(isfinite(tRange))
                mask = mask & (t >= min(tRange)) & (t <= max(tRange));
            end
            % Also need at least N+1 (or N+2 if hasInf) points to fit.
            nNL = numel(tauInit);
            if nnz(mask) < nNL + double(hasInf) + 1
                error(['Not enough data points for the fit ' ...
                       '(need at least %d, have %d after masking).'], ...
                       nNL + double(hasInf) + 1, nnz(mask));
            end
            
            % Pack/unpack of free nonlinear parameters.
            % We optimize log(tau), log(beta) (clamped), t0 directly, log(FWHM).
            freeIdxTau  = find(~tauFixed);
            freeIdxBeta = find(stretchOn & ~betaFixed);   % beta only free if stretched & not fixed
            t0Free   = ~t0Fixed;
            fwhmFree = ~fwhmFixed;
            
            x0 = [];
            if ~isempty(freeIdxTau)
                x0 = [x0, log(max(tauInit(freeIdxTau), eps))];
            end
            if ~isempty(freeIdxBeta)
                % beta in (0, 1.5]; map via log so optimizer stays positive.
                x0 = [x0, log(max(betaInit(freeIdxBeta), 1e-3))];
            end
            if t0Free,   x0 = [x0, t0Init];                 end
            if fwhmFree, x0 = [x0, log(max(fwhmInit, eps))]; end

            % Smallest FWHM the delays resolve, and the lower limit of
            % the search in 'limit' mode (see fitGlobalAnalysis); the
            % limit is in force during the search only.
            fwRes = 0;  fwLo = 0;  fwWall = false;
            if fwhmFree
                dzIrf = app.irfDelayStep(t(mask), t0Init);
                if isfinite(dzIrf) && dzIrf > 0
                    fwRes  = 0.5 * dzIrf;
                    fwLo   = min(fwRes, fwhmInit) * (1 - 1e-9);
                    fwWall = strcmp(fwMode, 'limit');
                end
            end
            
            tau0  = tauInit;
            beta0 = betaInit;
            t00   = t0Init;
            fw0   = fwhmInit;
            
            opts = optimset('Display', 'off', 'MaxIter', 1500, ...
                'MaxFunEvals', 5000, 'TolFun', 1e-8, 'TolX', 1e-7);
            
            [xHat, fval, exitFlag, output] = fminsearch( ...
                @(x) lossFn(x), x0, opts);
            fwWall = false;      % the search is over
            
            % Final params and fit vector
            [tauOut, betaOut, t0Out, fwhmOut, Aout, fitVec] = unpackFit(xHat);
            resVec = y - fitVec;
            % Residual is meaningful only where mask is true; outside,
            % set residual to NaN so plots clearly show the excluded
            % region instead of a misleading large spike.
            resVec(~mask) = NaN;
            rms = sqrt(mean(resVec(mask).^2, 'omitnan'));
            
            % ---- 1-sigma uncertainties (numerical Jacobian, physical
            %      scale). Same free-parameter ordering as x0. ----
            tauStd  = NaN(1, numel(tauInit));
            betaStd = NaN(1, numel(betaInit));
            t0Std   = NaN;   fwhmStd = NaN;
            if ~isempty(x0)
                pOpt = [];
                % tauOut / betaOut are row vectors: reshape (not .')
                % so that the pieces concatenate for any number of free
                % parameters.
                if ~isempty(freeIdxTau)
                    pOpt = [pOpt, reshape(tauOut(freeIdxTau), 1, [])];
                end
                if ~isempty(freeIdxBeta)
                    pOpt = [pOpt, reshape(betaOut(freeIdxBeta), 1, [])];
                end
                if t0Free,   pOpt = [pOpt, t0Out];   end
                if fwhmFree, pOpt = [pOpt, fwhmOut]; end
                
                % nData - (linear amplitudes); the nonlinear parameters
                % are subtracted inside fitParamUncertainty.
                sig = app.fitParamUncertainty(@(p) sfResidualAtPhys(p), ...
                    pOpt(:), fval, nnz(mask) - numel(Aout));
                
                ii2 = 1;
                for jj = 1:numel(freeIdxTau)
                    tauStd(freeIdxTau(jj)) = sig(ii2);  ii2 = ii2 + 1;
                end
                for jj = 1:numel(freeIdxBeta)
                    betaStd(freeIdxBeta(jj)) = sig(ii2); ii2 = ii2 + 1;
                end
                if t0Free,   t0Std   = sig(ii2);  ii2 = ii2 + 1; end
                if fwhmFree, fwhmStd = sig(ii2);                 end
            end
            
            info = struct('exitFlag', exitFlag, ...
                'iter', output.iterations, ...
                'fval', fval, ...
                'rms', rms, ...
                'mask', mask, ...
                'irfMode', irfMode, ...
                'tauStd', tauStd, 'betaStd', betaStd, ...
                't0Std', t0Std, 'fwhmStd', fwhmStd, ...
                'fwhmLo', NaN, 'fwhmAtLimit', false, 'irfRetry', 0);
            if fwhmFree && fwRes > 0
                info.fwhmLo      = fwRes;
                info.fwhmAtLimit = app.irfUnresolved(fwhmOut, fwRes);
            end

            % ---- two-step repeat when the IRF width was lost ----
            if fwhmFree && strcmp(fwMode, 'auto') && ...
                    app.irfWidthLost(fwhmOut, fwhmStd, fwRes)
                % A start value that the delays do not resolve gives
                % the search nothing to go by: two delay steps instead.
                fwS0 = fwhmInit;
                if app.irfUnresolved(fwhmInit, fwRes), fwS0 = 4 * fwRes; end
                try
                    [tauS1, betaS1, ~, ~, ~, ~, ~, infoS1] = ...
                        app.fitSingleTrace(t, y, tauInit, tauFixed, ...
                        betaInit, betaFixed, stretchOn, t0Init, true, ...
                        fwS0, true, hasInf, irfMode, tRange, 'plain');
                    [tauS2, betaS2, t0S2, fwS2, AS2, fitS2, resS2, infoS2] = ...
                        app.fitSingleTrace(t, y, tauS1, tauFixed, betaS1, ...
                        betaFixed, stretchOn, t0Init, t0Fixed, fwS0, ...
                        fwhmFixed, hasInf, irfMode, tRange, 'limit');
                    itAll = info.iter + infoS1.iter + infoS2.iter;
                    if infoS2.fval < info.fval * (1 - 1e-6) || ...
                            (info.fwhmAtLimit && ...
                             infoS2.fval <= info.fval * (1 + 1e-6))
                        tauOut = tauS2;  betaOut = betaS2;
                        t0Out  = t0S2;   fwhmOut = fwS2;
                        Aout = AS2;  fitVec = fitS2;  resVec = resS2;
                        info = infoS2;
                        info.irfRetry = 1;
                    else
                        info.irfRetry = 2;
                    end
                    info.iter = itAll;
                catch
                    % the repeat failed: the first result stands
                end
            end
            
            function r = sfResidualAtPhys(p)
                % Physical params -> internal packed x (log where the
                % optimizer used log), then masked residual.
                xx = zeros(1, numel(p));
                kk = 1;
                for q = 1:numel(freeIdxTau)
                    xx(kk) = log(max(p(kk), 1e-12));  kk = kk + 1;
                end
                for q = 1:numel(freeIdxBeta)
                    xx(kk) = log(max(p(kk), 1e-12));  kk = kk + 1;
                end
                if t0Free,   xx(kk) = p(kk);  kk = kk + 1;        end
                if fwhmFree, xx(kk) = log(max(p(kk), 1e-12));     end
                [~, ~, ~, ~, ~, fv_] = unpackFit(xx);
                r = y - fv_;
                r = r(mask);
            end
            
            % --------- nested helpers ---------
            function L = lossFn(x)
                [~, ~, ~, ~, ~, fv] = unpackFit(x);
                r = (y - fv);
                r = r(mask);
                L = sum(r.^2);
                if ~isfinite(L), L = 1e30; end
            end
            
            function [tauC, betaC, t0C, fwC, AC, fv] = unpackFit(x)
                tauC  = tau0;
                betaC = beta0;
                t0C   = t00;
                fwC   = fw0;
                idx = 1;
                if ~isempty(freeIdxTau)
                    tauC(freeIdxTau) = exp(x(idx : idx + numel(freeIdxTau) - 1));
                    idx = idx + numel(freeIdxTau);
                end
                if ~isempty(freeIdxBeta)
                    betaC(freeIdxBeta) = exp(x(idx : idx + numel(freeIdxBeta) - 1));
                    idx = idx + numel(freeIdxBeta);
                end
                if t0Free,   t0C = x(idx);              idx = idx + 1; end
                if fwhmFree, fwC = exp(x(idx)); end
                
                % Reasonable bounds: reject pathological regions cheaply.
                % (FWHM bounds are irrelevant in 'none' mode.)
                fwBad = ~strcmpi(irfMode, 'none') && ...
                    (~isfinite(fwC) || fwC <= 0 || ...
                     fwC > (max(t)-min(t)) || (fwWall && fwC < fwLo));
                if any(~isfinite(tauC)) || any(tauC <= 0) || ...
                        any(~isfinite(betaC)) || any(betaC <= 0) || ...
                        any(betaC > 2) || fwBad
                    AC = zeros(nNL + double(hasInf), 1);
                    fv = zeros(size(t));
                    return;
                end
                
                % Build basis matrix (only at masked points to keep loss
                % consistent with what we minimize).
                tm = t;
                M = zeros(numel(tm), nNL + double(hasInf));
                if strcmpi(irfMode, 'none')
                    % Pure kinetics without IRF convolution:
                    % exp(-((t-t0)/tau)^beta) . H(t-t0).
                    tt  = tm - t0C;
                    stp = double(tt >= 0);
                    for jj = 1:nNL
                        if stretchOn(jj)
                            M(:, jj) = stp .* ...
                                exp(-(max(tt,0) ./ tauC(jj)).^betaC(jj));
                        else
                            M(:, jj) = stp .* exp(-max(tt,0) ./ tauC(jj));
                        end
                    end
                    if hasInf
                        M(:, end) = stp;
                    end
                else
                    for jj = 1:nNL
                        if stretchOn(jj)
                            M(:, jj) = app.stretchedIRFConv(tm, tauC(jj), ...
                                betaC(jj), t0C, fwC);
                        else
                            M(:, jj) = app.expIRFConv(tm, tauC(jj), t0C, fwC);
                        end
                    end
                    if hasInf
                        M(:, end) = app.expIRFConv(tm, Inf, t0C, fwC);
                    end
                end
                
                Mm = M(mask, :);
                ym = y(mask);
                % lsqminnorm gives the minimum-norm least-squares
                % solution silently when the basis is rank-deficient
                % (two taus that have collapsed onto each other during
                % the search, or a tau far shorter than dt). The
                % backslash operator solves the same problem but prints
                % a "rank deficient" warning on every optimizer
                % iteration -- harmless, but it floods the command
                % window. Same choice as the GA engine.
                try
                    AC = lsqminnorm(Mm, ym);
                catch
                    AC = pinv(Mm) * ym;
                end
                if any(~isfinite(AC))
                    AC = pinv(Mm) * ym;
                end
                fv = M * AC;
            end
        end
        
        % =================================================================
        % Automatic analysis: proposal for the corrections + summary of
        % the kinetics without a kinetic model (calculation part)
        % =================================================================
        function o = autoOptions(~, o)
            %autoOptions  Options of the automatic analysis, with defaults.
            if nargin < 2 || ~isstruct(o), o = struct(); end
            d = struct( ...
                'doBg',      true, ...  % subtract the mean of the delays before the signal
                'doMask',    true, ...  % mask channels with excess noise
                'doChirp',   true, ...  % wavelength-dependent time zero
                'doCrop',    true, ...  % drop delays emptied by the chirp correction
                'maxComp',   5, ...     % largest number of time constants tried
                'minGain',   0.05, ...  % a further component must lower the RMS by this fraction
                'zOnset',    4, ...     % departure from the baseline, in noise levels
                'maskK',     3, ...     % noise excess over the neighbouring channels
                'maskSN',    10, ...    % ... masked only if the signal is below this many
                ...                     %   noise levels (Inf: every noisy channel)
                'maskSevere', 10, ...   % ... or the noise excess is above this
                'nChirpPts', 10, ...    % points used for the chirp curve
                'ca',        'auto', ...% coherent-artifact terms: 'auto' | 'off' | 'on'
                'caGain',    0.15, ...  % 'auto': RMS gain that calls for them
                'caTauMin',  1/3, ...   % shortest time constant told apart from an
                ...                     %   artifact, in units of the IRF FWHM
                'tryAll',    false);    % fit every number up to maxComp
            fn = fieldnames(d);
            for k = 1:numel(fn)
                if ~isfield(o, fn{k}) || isempty(o.(fn{k}))
                    o.(fn{k}) = d.(fn{k});
                end
            end
            o.maxComp = max(1, min(5, round(o.maxComp)));
            if ~ischar(o.ca) || ~any(strcmp(o.ca, {'auto', 'off', 'on'}))
                o.ca = 'auto';
            end
            o.tryAll = logical(o.tryAll(1));
        end

        function m = autoNanMean(~, X, dim)
            %autoNanMean  Mean along dim, ignoring NaN (NaN where none).
            ok = isfinite(X);
            X(~ok) = 0;
            n = sum(ok, dim);
            m = sum(X, dim) ./ n;
            m(n == 0) = NaN;
        end

        function y = autoRunMed(~, x, w)
            %autoRunMed  Running median over w points, ignoring NaN.
            x = x(:);  n = numel(x);  y = NaN(n, 1);
            h = floor(w / 2);
            for i = 1:n
                seg = x(max(1, i-h):min(n, i+h));
                seg = seg(isfinite(seg));
                if ~isempty(seg), y(i) = median(seg); end
            end
        end

        function r = autoAc1(~, x)
            %autoAc1  Lag-1 autocorrelation of a vector: near 0 for noise,
            %   near 1 for a smooth curve.
            x = x(:);  x = x(isfinite(x));
            r = NaN;
            if numel(x) < 3, return; end
            x = x - mean(x);
            den = sum(x.^2);
            if den > 0
                r = sum(x(1:end-1) .* x(2:end)) / den;
            end
        end

        function sg = autoNoise(~, D, t)
            %autoNoise  Noise level of every wavelength (row of D).
            %   Each point is compared with the straight line through its
            %   two neighbours on the delay axis. For a smooth signal the
            %   difference is noise only, and its median absolute
            %   deviation ignores the few points at the signal rise.
            %   Valid on a non-uniform delay axis.
            t = t(:).';
            M = size(D, 1);
            sg = NaN(M, 1);
            for i = 1:M
                m = isfinite(D(i, :));
                if nnz(m) < 8, continue; end
                y = D(i, m);  x = t(m);
                w1 = (x(3:end) - x(2:end-1)) ./ (x(3:end) - x(1:end-2));
                w2 = 1 - w1;
                r  = y(2:end-1) - (w1 .* y(1:end-2) + w2 .* y(3:end));
                r  = r ./ sqrt(1 + w1.^2 + w2.^2);
                r  = r(isfinite(r));
                if numel(r) < 5, continue; end
                sg(i) = 1.4826 * median(abs(r - median(r)));
            end
        end

        function R = autoRise(~, D, t, sg, nBase, zThr)
            %autoRise  Where the signal starts, wavelength by wavelength.
            %   Baseline = mean of the first nBase delays. The onset of a
            %   wavelength is the first delay at which the signal leaves
            %   the baseline by more than zThr noise levels at two
            %   delays in a row, on the same side.
            %   R.on        index of the onset (NaN: none)
            %   R.r30, r80  index at which 30 % / 80 % of the largest
            %               departure from the baseline is reached
            %   R.abrupt    true when 30 % is reached within 4 delays, or
            %               within 2 % of the time left after the onset
            [M, N] = size(D);
            t = t(:).';
            b = NaN(M, 1);
            for i = 1:M
                seg = D(i, 1:nBase);  seg = seg(isfinite(seg));
                if ~isempty(seg), b(i) = mean(seg); end
            end
            thr = zThr * sg * sqrt(1 + 1 / nBase);
            on  = NaN(M, 1);  r30 = NaN(M, 1);  r80 = NaN(M, 1);
            abrupt = false(M, 1);
            for i = 1:M
                if ~isfinite(sg(i)) || sg(i) <= 0 || ~isfinite(b(i)), continue; end
                d = D(i, :) - b(i);
                big = abs(d) > thr(i);          % NaN compares false
                j = nBase + 1;
                found = false;
                while j < N
                    if big(j) && big(j+1) && d(j) * d(j+1) > 0
                        found = true;  break;
                    end
                    j = j + 1;
                end
                if ~found, continue; end
                on(i) = j;
                ad = abs(d);  ad(~isfinite(ad)) = 0;
                mx = max(ad(j:end));
                k = find(ad(j:end) >= 0.3 * mx, 1) + j - 1;
                r30(i) = k;
                r80(i) = find(ad(j:end) >= 0.8 * mx, 1) + j - 1;
                abrupt(i) = (k - j) <= 4 || ...
                    (t(k) - t(j)) <= 0.02 * (t(end) - t(j));
            end
            R = struct('b', b, 'on', on, 'r30', r30, 'r80', r80, ...
                'abrupt', abrupt);
        end
        
        function [tr, zj, wid, dtk] = autoRidge(~, D, t, sg, jA, jB, w)
            %autoRidge  Time of the steepest change of every wavelength
            %   between delays jA and jB: the "max |dA/dt|" ridge of the
            %   chirp window, with the slope taken over w delay steps
            %   (w = 1: neighbouring delays).
            %   A plain rise has one slope peak; its position is refined
            %   between the grid points (parabola through the three
            %   largest slopes, or the centre of gravity of the part
            %   above half height when the peak is wide).
            %   A coherent artifact (a Gaussian, or one of its
            %   derivatives, around time zero) has several slope peaks,
            %   symmetric about its centre, and the largest slope alone
            %   would sit on one of its flanks. When further slopes of
            %   30 % of the largest, and clear of the noise, are found
            %   within three peak widths, the position is the centre of
            %   gravity of the squared slope over all of them.
            %   zj   size of the step at the steepest point, in noise levels
            %   wid  width (FWHM) of the instrument response read from it
            %   dtk  delay step at the steepest point
            t = t(:).';
            M = size(D, 1);
            tr = NaN(M, 1);  zj = zeros(M, 1);  wid = NaN(M, 1);
            dtk = NaN(M, 1);
            ts = t(jA:jB);
            nR = numel(ts) - w;
            if nR < 3, return; end
            dt = ts(1+w:end) - ts(1:end-w);
            tm = (ts(1+w:end) + ts(1:end-w)) / 2;
            d1 = diff(ts);
            for i = 1:M
                if ~isfinite(sg(i)) || sg(i) <= 0, continue; end
                y = D(i, jA:jB);
                dA = y(1+w:end) - y(1:end-w);
                rate = abs(dA ./ dt);
                rate(~isfinite(rate)) = 0;
                [mx, k] = max(rate);
                if ~(mx > 0), continue; end
                zj(i) = abs(dA(k)) / (sg(i) * sqrt(2));
                a = k;
                while a > 1 && rate(a-1) >= mx / 2, a = a - 1; end
                b = k;
                while b < nR && rate(b+1) >= mx / 2, b = b + 1; end
                dtk(i) = d1(min(k, numel(d1)));
                tr(i) = tm(k);
                wLobe = tm(b) - tm(a) + dtk(i);     % steepest slope peak
                a3 = k;
                while a3 > 1 && rate(a3-1) >= 0.3 * mx, a3 = a3 - 1; end
                b3 = k;
                while b3 < nR && rate(b3+1) >= 0.3 * mx, b3 = b3 + 1; end
                % Further slope peaks count only when they stand clear
                % of the noise (step of 4 noise levels): next to a weak
                % rise, noise alone reaches 30 % of the largest slope.
                sel = abs(tm - tm(k)) <= 3 * wLobe & rate >= 0.3 * mx & ...
                    abs(dA) >= 4 * sqrt(2) * sg(i);
                sel(a3:b3) = true;
                several = any(sel(1:a3-1)) || any(sel(b3+1:end));
                if several
                    wgt = rate(sel).^2;
                    tr(i) = sum(wgt .* tm(sel)) / sum(wgt);
                elseif b - a + 1 > 4
                    tr(i) = sum(rate(a:b) .* tm(a:b)) / sum(rate(a:b));
                elseif k > 1 && k < nR && ...
                        abs((tm(k) - tm(k-1)) - (tm(k+1) - tm(k))) < ...
                        0.01 * (tm(k+1) - tm(k))
                    den = rate(k-1) - 2 * rate(k) + rate(k+1);
                    if den < 0
                        tr(i) = tm(k) + 0.5 * (rate(k-1) - rate(k+1)) / ...
                            den * (tm(k+1) - tm(k));
                    end
                end
                % Width of the response: that of the slope peak of a
                % plain rise (the part that comes from taking the slope
                % over w steps removed), or -- wider -- the extent of
                % all slope peaks of an artifact, which spans about 1.9
                % times the width of its Gaussian.
                wRise = sqrt(max(wLobe^2 - 0.46 * dt(k)^2, dtk(i)^2));
                iS = find(sel);
                ext = tm(iS(end)) - tm(iS(1)) + dtk(i);
                wid(i) = max(wRise, ext / 1.9);
            end
        end
        
        function P = autoPreprocess(app, wl, t, D, opt)
            %autoPreprocess  Proposal for the corrections of a data set.
            %   wl (M values), t (N values), D (M x N, NaN allowed).
            %   Nothing in the app is changed. The result P describes
            %   each step (background, mask, chirp, crop) with the
            %   numbers it rests on and holds the corrected data:
            %     P.Dfull        corrected data on the input axes
            %     P.wl, P.t, P.D the same after the crop (and time shift)
            if nargin < 5, opt = struct(); end
            opt = app.autoOptions(opt);
            wl = wl(:);  t = t(:).';
            [M, N] = size(D);
            if M < 5 || N < 10
                error(['The data set is too small for the automatic ' ...
                    'analysis (%d x %d).'], M, N);
            end
            tu = app.tUnitHdr();  xu = app.xUnitTxt();
            checks = {};

            % ---------- noise of every wavelength ----------
            sg = app.autoNoise(D, t);
            rowOK = isfinite(sg) & sg > 0;
            if nnz(rowOK) < 5
                error('Fewer than 5 wavelengths hold usable data.');
            end

            % ---------- where the signal starts ----------
            % The baseline is first taken from 3 delays and then
            % extended to all delays found to lie before the signal.
            nBase = 3;  jOn = NaN;  nBg = 0;
            for it = 1:6
                R = app.autoRise(D, t, sg, nBase, opt.zOnset);
                det = rowOK & isfinite(R.on);
                if nnz(det) < max(3, ceil(0.05 * nnz(rowOK)))
                    jOn = NaN;  nBg = 0;  break;
                end
                onS = sort(R.on(det));
                jOn = onS(max(1, ceil(0.05 * numel(onS))));
                nBg = jOn - 1 - max(2, ceil(0.15 * (jOn - 1)));
                if nBg <= nBase, break; end
                nBase = nBg;
            end
            R = app.autoRise(D, t, sg, nBase, opt.zOnset);
            det = rowOK & isfinite(R.on);
            hasT0 = isfinite(jOn);
            t0Why = '';
            if ~hasT0
                t0Why = 'no signal onset was found';
            elseif nnz(det & R.abrupt) < 0.2 * nnz(det)
                hasT0 = false;
                t0Why = ['the signal changes only gradually after the ' ...
                    'first delays'];
            end
            
            % ---------- channels with excess noise ----------
            % Size of the signal of every channel: its largest departure
            % from the baseline (median of three neighbouring delays, so
            % that one outlier does not count), in units of its noise.
            b0 = R.b;  b0(~isfinite(b0)) = 0;
            Xa = abs(D - b0 * ones(1, N));
            Xa(~isfinite(Xa)) = 0;
            if N >= 3
                X3 = sort(cat(3, Xa(:, 1:N-2), Xa(:, 2:N-1), Xa(:, 3:N)), 3);
                sigAmp = max(X3(:, :, 2), [], 2);
            else
                sigAmp = max(Xa, [], 2);
            end
            snr = sigAmp ./ sg;  snr(~rowOK) = NaN;
            maskRows = false(M, 1);  maskStatus = 'off';
            ratio = NaN(M, 1);  noisy = false(M, 1);
            if opt.doMask
                maskStatus = 'none';
                w = max(11, 2 * round(0.1 * M) + 1);
                sgf = sg;  sgf(~rowOK) = NaN;
                ratio = sg ./ app.autoRunMed(sgf, w);
                noisy = rowOK & ratio > opt.maskK;
                % A noisy channel is dropped only when the noise is
                % severe: the signal does not stand clear of it, or it
                % is far above the neighbours. A channel that is merely
                % noisier than its neighbours but shows a strong signal
                % (probe edge, pump region with a large bleach) is kept.
                % opt.maskSN = Inf drops every noisy channel.
                sev  = snr < opt.maskSN | ratio > opt.maskSevere;
                sevW = snr < opt.maskSN | ratio > opt.maskSevere / 2;
                flag = noisy & sev;
                if any(flag)
                    % Extend every flagged run over the neighbours that
                    % are still clearly noisier than normal (half the
                    % threshold, at least 1.5 times) and meet the same
                    % condition, plus one channel.
                    lowK = max(1.5, opt.maskK / 2);
                    fl = flag;
                    for q = 2:M
                        if fl(q-1) && rowOK(q) && ratio(q) > lowK && sevW(q)
                            fl(q) = true;
                        end
                    end
                    for q = M-1:-1:1
                        if fl(q+1) && rowOK(q) && ratio(q) > lowK && sevW(q)
                            fl(q) = true;
                        end
                    end
                    flag = fl;
                    fl(1:end-1) = fl(1:end-1) | flag(2:end);
                    fl(2:end)   = fl(2:end)   | flag(1:end-1);
                    idx = find(fl);
                    for q = 1:numel(idx) - 1
                        if idx(q+1) - idx(q) <= 3
                            fl(idx(q):idx(q+1)) = true;
                        end
                    end
                    if nnz(fl) > 0.25 * M
                        maskStatus = 'toomany';
                    else
                        maskStatus = 'applied';
                        maskRows = fl;
                    end
                end
            end
            % Runs of masked channels; a run that touches an end of the
            % spectrum is cut away by the crop instead of being masked.
            rs = find(diff([false; maskRows]) == 1);
            re = find(diff([maskRows; false]) == -1);
            edgeLo = 0;  edgeHi = 0;
            if opt.doCrop && ~isempty(rs)
                if rs(1) == 1, edgeLo = re(1); end
                if re(end) == M && (numel(rs) > 1 || edgeLo == 0)
                    edgeHi = rs(end);
                end
            end
            maskReg = zeros(0, 4);     % [from, to, largest noise ratio,
                                       %  smallest signal / noise]
            for q = 1:numel(rs)
                if (edgeLo > 0 && rs(q) == 1) || (edgeHi > 0 && re(q) == M)
                    continue;
                end
                seg = rs(q):re(q);
                maskReg(end+1, :) = [min(wl(seg)), max(wl(seg)), ...
                    max(ratio(seg)), min(snr(seg))]; %#ok<AGROW>
            end
            % Noisy channels that stay in the data
            keptReg = zeros(0, 4);     % same columns
            kn = noisy & ~maskRows;
            ks = find(diff([false; kn]) == 1);
            ke = find(diff([kn; false]) == -1);
            for q = 1:numel(ks)
                seg = ks(q):ke(q);
                keptReg(end+1, :) = [min(wl(seg)), max(wl(seg)), ...
                    max(ratio(seg)), min(snr(seg))]; %#ok<AGROW>
            end
            vr = rowOK & ~maskRows;
            
            % ---------- the rise of every wavelength ----------
            ch = struct('status', 'off', 'msg', '', 'pts', zeros(0, 2), ...
                'p', [], 'rms', NaN, 't0', zeros(M, 1), 'span', NaN, ...
                'fracIn', NaN, 'scatter', NaN, 'tol', NaN, ...
                'ridgeT', NaN(M, 1), 'good', false(M, 1), ...
                'win', [NaN NaN], 'dtZ', NaN, 'wStep', 1, 'cover', [NaN NaN], ...
                'model', '', ...
                'shift', 0, 'monotonic', true, 'nRejected', 0);
            fwGuess = NaN;  t0Guess = NaN;
            tr = NaN(M, 1);  good = false(M, 1);
            if hasT0
                use = vr & isfinite(R.on) & R.abrupt;
                if nnz(use) >= 5
                    % Search window: from just before the earliest onset
                    % to just after the latest rise. The slope is taken
                    % over as many delay steps as it needs to stand
                    % clearly above the noise (finely sampled or noisy
                    % data: more than one).
                    onU = sort(R.on(use));
                    j0  = onU(max(1, round(0.01 * numel(onU))));
                    r8  = sort(R.r80(use));
                    j1  = r8(max(1, ceil(0.98 * numel(r8))));
                    wList = [1 2 3 4 6 8 12 16 24 32];
                    nNeed = max(5, ceil(0.3 * nnz(use)));
                    bestN = -1;
                    for q = 1:numel(wList)
                        ws = wList(q);
                        jA = max(1, j0 - 5 - ws);
                        jB = min(N, j1 + 3 + ws);
                        if jB - jA < 4 * ws
                            if q > 1, break; end
                            jB = min(N, jA + 4);
                        end
                        [trQ, zjQ, widQ, dtkQ] = app.autoRidge(D, t, sg, jA, jB, ws);
                        nClear = nnz(use & zjQ > 6);
                        if nClear > bestN
                            bestN = nClear;
                            tr = trQ;  zj = zjQ;  wid = widQ;  dtk = dtkQ;
                            ch.win = [t(jA), t(jB)];  ch.wStep = ws;
                        end
                        % enough when 30 % of these wavelengths show a
                        % step of more than 6 noise levels
                        if nClear >= nNeed, break; end
                    end
                    good = use & isfinite(tr) & zj > 5;
                    ch.ridgeT = tr;  ch.good = good;
                    if nnz(good) >= 5
                        ch.dtZ  = median(dtk(good));
                        fwGuess = median(wid(good));
                        t0Guess = median(tr(good));
                    end
                end
                % A time zero inside the data shows as a rise after a flat
                % start, or -- with fewer than 3 delays before it -- as a
                % jump between neighbouring delays.
                if nnz(good) < max(5, ceil(0.05 * nnz(vr))) || ...
                        (nBg < 3 && ch.wStep > 1)
                    hasT0 = false;
                    t0Why = ['no sharp rise: the signal seems to be there ' ...
                        'from the first delays on'];
                    fwGuess = NaN;  t0Guess = NaN;
                    good = false(M, 1);  ch.good = good;
                end
            end
            
            % ---------- background ----------
            bg = zeros(M, 1);  bgStatus = 'off';  bSig = NaN;
            if opt.doBg
                if ~hasT0
                    bgStatus = 'notime0';
                elseif nBg < 3
                    bgStatus = 'few';
                else
                    b = app.autoNanMean(D(:, 1:nBg), 2);
                    e = app.autoNanMean(D(:, max(1, N-2):N), 2);
                    m = rowOK & isfinite(b) & isfinite(e);
                    bw = b(m) ./ sg(m);  ew = e(m) ./ sg(m);
                    dw = ew - bw;
                    nb = norm(bw);
                    bSig = nb / sqrt(max(1, nnz(m))) * sqrt(nBg);
                    cosv = (bw.' * dw) / max(nb * norm(dw), realmin);
                    if bSig > 3 && norm(ew) < 0.5 * nb && cosv < -0.8
                        % The first delays hold more signal than the last
                        % ones and the data move away from them towards
                        % zero: a decaying signal, not a background.
                        bgStatus = 'signal';
                    elseif bSig < 2
                        bgStatus = 'zero';
                    else
                        bgStatus = 'applied';
                        bg = b;  bg(~isfinite(bg)) = 0;
                    end
                end
            end
            Db = D - bg;
            
            % ---------- time zero along the spectrum (chirp) ----------
            if ~hasT0
                ch.status = 'notime0';
            else
                % Points of the curve: median ridge time in equal
                % wavelength intervals (robust against ridge jumps).
                nB = min(opt.nChirpPts, floor(nnz(good) / 3));
                pts = zeros(0, 2);
                if nB >= 5
                    edges = linspace(min(wl(vr)), max(wl(vr)), nB + 1);
                    for k = 1:nB
                        if k < nB
                            inb = wl >= edges(k) & wl <  edges(k+1);
                        else
                            inb = wl >= edges(k) & wl <= edges(k+1);
                        end
                        g = inb & good;
                        if nnz(g) >= max(3, ceil(0.3 * nnz(inb & vr)))
                            pts(end+1, :) = [median(wl(g)), median(tr(g))]; %#ok<AGROW>
                        end
                    end
                else
                    pts = [wl(good), tr(good)];
                end
                % Deviations that do not matter: two delay steps, or
                % 30 % of the width of the rise.
                ch.tol = max(2 * ch.dtZ, 0.3 * fwGuess);
                if size(pts, 1) < 5
                    ch.status = 'failed';
                    ch.msg = 'fewer than 5 usable points along the spectrum';
                else
                    % Curve 1: the chirp formula of the Chirp Corr window.
                    % Curve 2, only if the first gives nothing sensible
                    % (typically when the rise is seen in a small part
                    % of the spectrum): a parabola.
                    for attempt = 1:2
                        try
                            if attempt == 1
                                % First pass through the interval medians,
                                % then through every wavelength that agrees
                                % with that curve (the others: ridge jumps).
                                [p, rmsF] = app.fitChirpParams(pts);
                                r = tr - app.chirpModel(p, wl);
                                rg = r(good);
                                s = 1.4826 * median(abs(rg - median(rg)));
                                inl = good & abs(r - median(rg)) <= ...
                                    max(3 * s, 0.75 * ch.tol);
                                if nnz(inl) >= 8
                                    pts = [wl(inl), tr(inl)];
                                    [p, rmsF] = app.fitChirpParams(pts);
                                end
                                t0m = app.chirpModel(p, wl);
                                ch.model = 'formula';
                            else
                                xm = mean(wl(good));  xs = std(wl(good));
                                xg = (wl(good) - xm) / xs;
                                Ap = [ones(size(xg)), xg, xg.^2];
                                cp = Ap \ tr(good);
                                rg = tr(good) - Ap * cp;
                                s = 1.4826 * median(abs(rg - median(rg)));
                                keepP = abs(rg - median(rg)) <= ...
                                    max(3 * s, 0.75 * ch.tol);
                                if nnz(keepP) < 8, keepP = true(size(rg)); end
                                gi = find(good);
                                pts = [wl(gi(keepP)), tr(gi(keepP))];
                                cp = Ap(keepP, :) \ pts(:, 2);
                                xa = (wl - xm) / xs;
                                t0m = [ones(M, 1), xa, xa.^2] * cp;
                                p = [];
                                rmsF = sqrt(mean((pts(:, 2) - ...
                                    interp1(wl, t0m, pts(:, 1))).^2));
                                ch.model = 'parabola';
                            end
                            if any(~isfinite(t0m(vr)))
                                error('the fitted curve is not defined over the whole spectrum');
                            end
                            % Where no wavelength shows a rise the curve is
                            % an extrapolation. It may leave the range of
                            % the measured rises by half that range, not
                            % more.
                            lo = min(pts(:, 2));  hi = max(pts(:, 2));
                            slack = 0.5 * (hi - lo) + ch.tol;
                            if any(t0m(vr) < lo - slack | t0m(vr) > hi + slack)
                                error(['the fitted curve runs far outside the ' ...
                                    'measured rises (%.3g ... %.3g %s)'], lo, hi, tu);
                            end
                            res = tr(good) - t0m(good);
                            ch.nRejected = nnz(good) - size(pts, 1);
                            ch.cover = [min(pts(:, 1)), max(pts(:, 1))];
                            ch.pts = pts;  ch.p = p;  ch.rms = rmsF;
                            ch.fracIn  = nnz(abs(res) <= ch.tol) / numel(res);
                            ch.scatter = 1.4826 * median(abs(res - median(res)));
                            ch.span = max(t0m(vr)) - min(t0m(vr));
                            d1 = diff(t0m(vr));
                            ch.monotonic = all(d1 >= 0) || all(d1 <= 0);
                            % Is there a chirp at all? Judged where rises
                            % were measured, not where the curve is an
                            % extrapolation.
                            inC = vr & wl >= ch.cover(1) & wl <= ch.cover(2);
                            spanC = max(t0m(inC)) - min(t0m(inC));
                            if spanC < 0.75 * ch.tol
                                ch.status = 'none';
                                ch.span = spanC;
                            elseif ch.fracIn < 0.7
                                error(['only %.0f %% of the wavelengths lie ' ...
                                    'within %.2g %s of the fitted curve'], ...
                                    100 * ch.fracIn, ch.tol, tu);
                            elseif ~opt.doChirp
                                ch.status = 'found';
                                ch.t0 = t0m;
                            else
                                ch.status = 'applied';
                                ch.t0 = t0m;
                            end
                            ch.msg = '';
                            break;
                        catch errC
                            ch.status = 'failed';
                            if attempt == 1, ch.msg = errC.message; end
                        end
                    end
                end
            end
            
            Dc = Db;  tOut = t;
            switch ch.status
                case 'applied'
                    % Same operation as the chirp correction of the main
                    % window (recompute).
                    try
                        for i = 1:M
                            if ~isfinite(ch.t0(i)), Dc(i, :) = NaN;  continue; end
                            Dc(i, :) = interp1(t, Db(i, :), t + ch.t0(i), ...
                                'linear', NaN);
                        end
                        t0Guess = 0;
                    catch errI
                        % e.g. two identical delay values
                        Dc = Db;
                        ch.status = 'failed';
                        ch.msg = ['the correction could not be carried ' ...
                            'out: ' errI.message];
                    end
                case 'none'
                    % One time zero for all wavelengths: move it to 0
                    % when it is off by more than one delay step.
                    if opt.doChirp && isfinite(t0Guess) && ...
                            abs(t0Guess) > 0.5 * ch.tol
                        ch.shift = t0Guess;
                        tOut = t - ch.shift;
                        t0Guess = 0;
                    end
            end
            Dc(maskRows, :) = NaN;

            % ---------- crop ----------
            keepT = true(1, N);  keepW = true(M, 1);
            if opt.doCrop
                colOK = all(isfinite(Dc(vr, :)), 1);
                if any(colOK) && ~all(colOK)
                    j1 = find(colOK, 1);  j2 = find(colOK, 1, 'last');
                    if j2 - j1 + 1 >= 10
                        keepT(1:j1-1) = false;  keepT(j2+1:end) = false;
                    end
                end
                if edgeLo > 0, keepW(1:edgeLo) = false; end
                if edgeHi > 0, keepW(edgeHi:M) = false; end
            end

            % ---------- things to check by eye ----------
            if ~hasT0
                checks{end+1} = sprintf(['No time zero was found in the data ' ...
                    '(%s). Background and chirp were left alone, and the ' ...
                    'fits use no instrument response.'], t0Why);
            end
            switch bgStatus
                case 'few'
                    checks{end+1} = sprintf(['Only %d delay(s) lie safely before ' ...
                        'the signal: no background was subtracted. Set it by ' ...
                        'hand if there is one.'], max(nBg, 0));
                case 'signal'
                    checks{end+1} = ['The first delays look like a signal that ' ...
                        'decays later, not like a background: nothing was ' ...
                        'subtracted.'];
            end
            if strcmp(maskStatus, 'toomany')
                checks{end+1} = ['More than a quarter of the channels are much ' ...
                    'noisier than their neighbours: no mask was set. Look at ' ...
                    'the noise plot.'];
            end
            if ~isempty(keptReg)
                checks{end+1} = ['Some channels are noisier than their ' ...
                    'neighbours but were kept, because their signal is far ' ...
                    'above that noise (listed under Mask). The fits are not ' ...
                    'weighted, so such channels count more than quiet ones. ' ...
                    'If the residual shows stripes there, mask them: Mask ' ...
                    'rule "Every noisy channel", or Mask in the main window.'];
            end
            switch ch.status
                case 'failed'
                    checks{end+1} = sprintf(['The time zero along the spectrum ' ...
                        'could not be determined reliably (%s): no chirp ' ...
                        'correction. Use the Chirp Corr window.'], ch.msg);
                case {'applied', 'found'}
                    if strcmp(ch.model, 'parabola')
                        checks{end+1} = ['The chirp formula of the Chirp Corr ' ...
                            'window gave no sensible curve here; a parabola ' ...
                            'through the rises was used instead. Compare the ' ...
                            'curve with the map.'];
                    end
                    if ~ch.monotonic
                        checks{end+1} = ['The fitted time-zero curve is not ' ...
                            'monotonic along the spectrum. Chirp from material ' ...
                            'dispersion is: compare the curve with the map.'];
                    end
                    if ch.fracIn < 0.9
                        checks{end+1} = sprintf(['Only %.0f %% of the wavelengths ' ...
                            'lie within %.2g %s of the time-zero curve: ' ...
                            'compare the curve with the map.'], ...
                            100 * ch.fracIn, ch.tol, tu);
                    end
                    wr = max(wl(vr)) - min(wl(vr));
                    if ch.cover(1) - min(wl(vr)) > 0.1 * wr || ...
                            max(wl(vr)) - ch.cover(2) > 0.1 * wr
                        checks{end+1} = sprintf(['A sharp rise is seen only from ' ...
                            '%.4g to %.4g %s. Outside, the time-zero curve is an ' ...
                            'extrapolation.'], ch.cover(1), ch.cover(2), xu);
                    end
            end

            P = struct();
            P.opt = opt;
            P.wlIn = wl;  P.tIn = t;
            P.noise = sg;  P.rowOK = rowOK;  P.noiseRatio = ratio;
            P.hasT0 = hasT0;  P.t0Why = t0Why;
            P.jOn = jOn;
            if hasT0, P.tOn = t(jOn); else, P.tOn = NaN; end
            P.nBg = nBg;  P.bg = bg;  P.bgStatus = bgStatus;  P.bgSig = bSig;
            P.maskRows = maskRows;  P.maskReg = maskReg;
            P.maskStatus = maskStatus;
            P.keptReg = keptReg;  P.snr = snr;  P.noisy = noisy;
            P.chirp = ch;
            P.keepT = keepT;  P.keepW = keepW;
            P.cropT = [];  P.cropW = [];
            if any(~keepT), P.cropT = [t(find(keepT, 1)), t(find(keepT, 1, 'last'))]; end
            if any(~keepW), P.cropW = [min(wl(keepW)), max(wl(keepW))]; end
            P.Dfull = Dc;
            P.wl = wl(keepW);  P.t = tOut(keepT);  P.D = Dc(keepW, keepT);
            P.noiseOut = sg(keepW);
            P.t0Guess = t0Guess;  P.fwGuess = fwGuess;
            P.changed = strcmp(bgStatus, 'applied') || any(maskRows) || ...
                strcmp(ch.status, 'applied') || ch.shift ~= 0 || ...
                any(~keepT) || any(~keepW);
            P.checks = checks;
            P.tu = tu;  P.xu = xu;
        end

        function F = autoFitGA(app, D, t, tau0, t0s, fws, hasInf, irf, fwFixed, caK, irfFixed)
            %autoFitGA  One global fit for the automatic analysis; [] when
            %   it fails. Same fit as "Run Fit" of the Global Analysis
            %   window with tau, t0 and (unless fwFixed) FWHM free, and
            %   with caK coherent-artifact terms of the IRF width. With
            %   irfFixed, t0 and FWHM both stay at t0s and fws. The time
            %   constants come back in ascending order. The fit is the
            %   plain one, without the two-step repeat of
            %   fitGlobalAnalysis: a width below half a delay step is
            %   dealt with by the sequence itself (autoFitSequence).
            if nargin < 10, caK = 0; end
            if nargin < 11, irfFixed = false; end
            F = [];
            isNone = strcmp(irf, 'none');
            try
                [tau, ~, t0, fw, A, fitM, info] = app.fitGlobalAnalysis( ...
                    D, t, tau0, t0s, fws, false, isNone || irfFixed, ...
                    isNone || fwFixed || irfFixed, hasInf, [], [], irf, [], ...
                    struct('order', caK, 'mode', 'irf'), 'plain');
            catch
                return;
            end
            if ~isfinite(info.rss) || info.rss >= 1e29, return; end
            [tau, ix] = sort(tau(:).');
            nT = numel(tau);
            perm = ix;
            if hasInf, perm = [ix, nT + 1]; end
            F = struct('tau', tau, 'tauStd', info.tauStd(ix), ...
                't0', t0, 'fwhm', fw, 't0Std', info.t0Std, ...
                'fwhmStd', info.fwhmStd, 'A', A(:, perm), 'fit', fitM, ...
                'rss', info.rss, 'rms', info.rms, 'hasInf', hasInf, ...
                'caK', info.caOrder, 'caSpec', info.caSpec, ...
                'irfFixed', irfFixed && ~isNone);
        end
        
        function S = autoSummarize(app, P, opt, prog)
            %autoSummarize  Summary of the kinetics without a kinetic
            %   model, for the corrected data in P (see autoPreprocess):
            %   singular values, global fits with 1, 2, ... components
            %   (autoFitSequence), with coherent-artifact terms when the
            %   data call for them, choice of the number of components,
            %   and the model of that number (autoSelectModel).
            %   prog(fraction, text) is called as the work proceeds.
            if nargin < 3, opt = struct(); end
            if nargin < 4, prog = []; end
            opt = app.autoOptions(opt);
            wl = P.wl(:);  t = P.t(:).';  D = P.D;

            % Rectangular block without missing values
            rows = any(isfinite(D), 2);
            cols = all(isfinite(D(rows, :)), 1);
            if nnz(rows) < 5 || nnz(cols) < max(10, 0.5 * numel(t))
                error(['Too many missing values for the automatic ' ...
                    'analysis. Crop the data to the region that is ' ...
                    'filled completely.']);
            end
            Dv = D(rows, cols);  tv = t(cols);  wv = wl(rows);
            Nv = numel(tv);
            sgv = P.noiseOut(rows);
            S = struct();
            S.rows = rows;  S.cols = cols;  S.t = tv;  S.wl = wv;
            S.rmsNoise = sqrt(mean(sgv(isfinite(sgv)).^2));

            % ---------- singular values ----------
            [U, Sg, V] = svd(Dv, 'econ');
            sv = diag(Sg);
            nsv = numel(sv);
            fl = median(sv(ceil(nsv / 2):end));
            if ~(fl > 0), fl = eps; end
            kMax = min(10, nsv);
            acU = zeros(1, kMax);  acV = zeros(1, kMax);
            for k = 1:kMax
                acU(k) = app.autoAc1(U(:, k));
                acV(k) = app.autoAc1(V(:, k));
            end
            nStruct = 0;
            while nStruct < kMax && acU(nStruct+1) > 0.5 && acV(nStruct+1) > 0.5
                nStruct = nStruct + 1;
            end
            S.sv = sv;  S.svFloor = fl;
            S.nAbove = nnz(sv > 5 * fl);
            S.acU = acU;  S.acV = acV;  S.nStruct = nStruct;

            % ---------- instrument response ----------
            dAll = diff(tv);
            dtMin = min(dAll(dAll > 0));
            dtZ = P.chirp.dtZ;
            hasT0 = P.hasT0 && isfinite(P.t0Guess) && nnz(tv < P.t0Guess) >= 3;
            if hasT0
                irf = 'numerical';
                t0s = P.t0Guess;
                if ~isfinite(dtZ) || dtZ <= 0, dtZ = dtMin; end
                fws = P.fwGuess;
                if ~isfinite(fws) || fws <= 0, fws = 2 * dtZ; end
                fws = max(fws, dtZ);
                tauLo = max(0.5 * fws, dtZ);
            else
                irf = 'none';
                t0s = min(0, tv(1));
                dtZ = dtMin;
                fws = dtMin;
                tauLo = 2 * dtMin;
            end
            tEnd  = tv(end) - t0s;
            tauHi = max(0.5 * tEnd, 8 * tauLo);
            S.irf = irf;  S.tEnd = tEnd;  S.dtZ = dtZ;

            % Data reduced to the leading singular vectors: the same
            % kinetics at a fraction of the cost. Used to compare
            % starting values; every result quoted is then refitted on
            % the full data.
            kc = min(10, nsv);
            Dk = Sg(1:kc, 1:kc) * V(:, 1:kc).';

            C = struct('irf', irf, 't0s', t0s, 'fws', fws, 'dtZ', dtZ, ...
                'tauLo', tauLo, 'tauHi', tauHi, 'tEnd', tEnd, 'tu', P.tu);
            
            % ---------- global fits, without artifact terms ----------
            T = app.autoFitSequence(Dv, Dk, tv, C, opt, 0, prog, 0.15, 0.50);
            if T.nStar == 0
                error('The global fits of the automatic analysis failed.');
            end
            
            % ---------- coherent artifact ----------
            % Three artifact terms: a Gaussian of the IRF width at time
            % zero and its first two derivatives, each with a free
            % spectrum. 'auto' uses them in two cases.
            % (1) The fit chosen without them needs a time constant
            %     below opt.caTauMin x IRF FWHM. Convolved with the IRF,
            %     such a decay is a Gaussian at time zero: it cannot be
            %     told from an artifact, and is not counted as a
            %     component.
            % (2) Added to the chosen fit, the terms lower the RMS by
            %     opt.caGain or more: an artifact left out distorts the
            %     IRF of the fit.
            % The fit with the terms starts from the chosen fit (without
            % its time constants of case 1). The components are then
            % counted again with the terms in place (autoFitSequence).
            ca = struct('mode', opt.ca, 'use', 0, 'why', '', 'gain', NaN, ...
                'note', '', 'tauFast', NaN, 'fw0', NaN, 'n0', T.nStar, ...
                'rms0', T.Rms(T.nStar), 'nSeed', 0);
            nPre = 0;
            if ~strcmp(irf, 'none')
                ca.fw0 = T.F{T.nStar}.fwhm;
                nPre = nnz(tv < T.F{T.nStar}.t0 - 1.5 * ca.fw0);
            end
            if strcmp(irf, 'none')
                ca.note = 'not applicable: no time zero in the data';
            elseif strcmp(opt.ca, 'auto') && nPre < 5
                % Without a baseline before the rise, an offset left in
                % the data looks like a signal at time zero.
                ca.note = sprintf(['not tried: only %d delay(s) before the ' ...
                    'rise, so a baseline error could not be told from ' ...
                    'an artifact'], nPre);
            elseif ~strcmp(opt.ca, 'off')
                if ~isempty(prog), prog(0.52, 'Coherent artifact...'); end
                F0 = T.F{T.nStar};
                fast = F0.tau < opt.caTauMin * F0.fwhm;
                tauK = F0.tau(~fast);
                if isempty(tauK), tauK = F0.tau; end
                F3 = [];
                for pass = 1:3
                    F3 = app.autoFitGA(Dv, tv, tauK, F0.t0, F0.fwhm, ...
                        F0.hasInf, irf, T.fwFixed, 3);
                    if isempty(F3) || numel(F3.tau) < 2, break; end
                    % Next to the terms the fastest time constant may no
                    % longer be supported by the data: leave it out.
                    rel1 = F3.tauStd(1) / F3.tau(1);
                    if F3.tau(1) >= opt.caTauMin * F3.fwhm && ...
                            isfinite(rel1) && rel1 <= 0.5
                        break;
                    end
                    tauK = F3.tau(2:end);
                end
                if isempty(F3)
                    ca.note = ['the fit with artifact terms failed (too ' ...
                        'few delays around time zero?)'];
                elseif F3.t0 < tv(1) || F3.t0 > tv(end) || ...
                        F3.fwhm > 3 * F0.fwhm || F3.fwhm < F0.fwhm / 3
                    % The terms left time zero: with a width of their
                    % own choosing they stand in for kinetics.
                    ca.note = ['with artifact terms the fitted IRF does ' ...
                        'not stay at time zero'];
                else
                    ca.gain = 1 - F3.rms / F0.rms;
                    why = '';
                    if strcmp(opt.ca, 'on')
                        why = 'on';
                    elseif any(fast)
                        why = 'fast';
                        ca.tauFast = min(F0.tau(fast));
                    elseif ca.gain >= opt.caGain
                        why = 'gain';
                    end
                    if ~isempty(why)
                        C.seed = struct('n', numel(F3.tau) + double(F3.hasInf), ...
                            'F', F3, 'tau0', {T.Tau}, 'inf0', T.Inf, ...
                            'fwFixed', T.fwFixed, 'fwNote', T.fwNote);
                        T3 = app.autoFitSequence(Dv, Dk, tv, C, opt, 3, ...
                            prog, 0.55, 0.85);
                        if T3.nStar > 0
                            T = T3;  ca.use = 3;  ca.why = why;
                            ca.nSeed = C.seed.n;
                        else
                            ca.note = 'the fits with artifact terms failed';
                        end
                    end
                end
            end
            S.ca = ca;
            S.tabN = T.N;  S.tabRms = T.Rms;  S.tabGain = T.Gain;
            S.tabInf = T.Inf;  S.tabTau = T.Tau;  S.tabStd = T.Std;
            S.tabWhy = T.Why;  S.tabF = T.F;
            S.nStar = T.nStar;  S.fwNote = T.fwNote;  S.fwFixed = T.fwFixed;
            S.stopWhy = T.stopWhy;
            S.hitMax = isempty(T.stopWhy);
            S.minGain = opt.minGain;
            S.caGainMin = opt.caGain;
            
            if ~isempty(prog), prog(0.9, 'Chosen model...'); end
            S = app.autoSelectModel(P, S, S.nStar);
            if ~isempty(prog), prog(1, 'Done.'); end
        end
        
        function T = autoFitSequence(app, Dv, Dk, tv, C, opt, caK, prog, p0, p1)
            %autoFitSequence  Global fits with 1, 2, ... components.
            %   Dv data, Dk the same reduced to its leading singular
            %   vectors (used to compare starting values; every fit
            %   quoted is then refitted on Dv), C the settings prepared
            %   by autoSummarize, caK the number of coherent-artifact
            %   terms in every fit.
            %   A further component is rejected when it lowers the RMS
            %   by less than opt.minGain or brings a time constant the
            %   data do not determine; next to artifact terms also when
            %   a time constant falls below opt.caTauMin x IRF FWHM
            %   (the terms describe a decay that short equally well).
            %   T.nStar is the last number before the first rejection.
            %   The fits stop there unless opt.tryAll is set.
            %   With artifact terms (C.seed, see autoSummarize) the count
            %   starts at the fit the terms were found with, C.seed.F
            %   with C.seed.n components. Fits with fewer components are
            %   listed for comparison only, with t0 and FWHM kept at the
            %   values of that fit: in a model with too few components,
            %   free artifact terms leave time zero, widen and stand in
            %   for the missing kinetics.
            irf = C.irf;  dtZ = C.dtZ;  tEnd = C.tEnd;
            tauLo = C.tauLo;  tauHi = C.tauHi;
            nMax = opt.maxComp;
            tabN = [];  tabRms = [];  tabGain = [];  tabInf = [];
            tabTau = {};  tabStd = {};  tabWhy = {};  tabF = {};
            fwFixed = false;  fwNote = '';
            useInf = false;        % one component that does not decay
            prevTau = [];  t0c = C.t0s;  fwc = C.fws;
            nStar = 0;  stopWhy = '';  decided = false;
            seed = [];
            if caK > 0 && isfield(C, 'seed'), seed = C.seed; end
            n = 1;
            while n <= nMax
                if ~isempty(prog)
                    if caK > 0
                        msg = sprintf(['Global fit with %d component(s) ' ...
                            'and artifact terms...'], n);
                    else
                        msg = sprintf('Global fit with %d component(s)...', n);
                    end
                    prog(p0 + (p1 - p0) * (n - 1) / nMax, msg);
                end
                below = ~isempty(seed) && n < seed.n;
                if below
                    % comparison fit: structure of the fit without the
                    % terms, IRF of the seed fit
                    F = app.autoFitGA(Dv, tv, seed.tau0{n}, seed.F.t0, ...
                        seed.F.fwhm, logical(seed.inf0(n)), irf, true, ...
                        caK, true);
                    if isempty(F)
                        stopWhy = sprintf('the fit with %d components failed', n);
                        break;
                    end
                    gain = NaN;
                    if n > 1, gain = 1 - F.rms / tabRms(n-1); end
                    tabN(n) = n;  tabRms(n) = F.rms;  tabGain(n) = gain; %#ok<AGROW>
                    tabInf(n) = F.hasInf; %#ok<AGROW>
                    tabTau{n} = F.tau;  tabStd{n} = F.tauStd; %#ok<AGROW>
                    tabWhy{n} = '';  tabF{n} = F; %#ok<AGROW>
                    n = n + 1;
                    continue;
                end
                atSeed = ~isempty(seed) && n == seed.n;
                if atSeed
                    useInf = seed.F.hasInf;
                    fwFixed = seed.fwFixed;  fwNote = seed.fwNote;
                end
                m = n - double(useInf);     % time constants in this fit
                st = {};
                if atSeed
                    st{1} = seed.F.tau;
                    t0c = seed.F.t0;  fwc = seed.F.fwhm;
                elseif isempty(prevTau)
                    st{1} = max(tEnd / 5, 2 * tauLo);
                    st{2} = sqrt(2 * tauLo * tauHi);
                else
                    % Starting values: evenly spread on a log scale, and
                    % the previous result plus one new time constant (in
                    % its widest gap, below it, above it).
                    st{1} = logspace(log10(2 * tauLo), log10(tauHi), m);
                    pv = sort(prevTau);
                    gr = [tauLo, pv, max(2 * tEnd, 2 * pv(end))];
                    [~, g] = max(diff(log(gr)));
                    st{end+1} = sort([pv, sqrt(gr(g) * gr(g+1))]); %#ok<AGROW>
                    if pv(1) / 4 > tauLo
                        st{end+1} = [pv(1) / 4, pv]; %#ok<AGROW>
                    end
                    if pv(end) * 4 < 3 * tEnd
                        st{end+1} = [pv, pv(end) * 4]; %#ok<AGROW>
                    end
                end
                best = [];
                if atSeed, st = {}; end
                for s = 1:numel(st)
                    Fs = app.autoFitGA(Dk, tv, st{s}, t0c, fwc, useInf, ...
                        irf, fwFixed, caK);
                    if ~isempty(Fs) && (isempty(best) || Fs.rss < best.rss)
                        best = Fs;
                    end
                end
                F = [];
                if atSeed
                    F = seed.F;
                elseif ~isempty(best)
                    F = app.autoFitGA(Dv, tv, best.tau, best.t0, best.fwhm, ...
                        useInf, irf, fwFixed, caK);
                end
                if isempty(F)
                    if ~decided
                        stopWhy = sprintf('the fit with %d components failed', n);
                    end
                    break;
                end
                if strcmp(irf, 'numerical') && ~fwFixed && ~atSeed && ...
                        F.fwhm < 0.5 * dtZ
                    % The fitted response is narrower than half a delay
                    % step: not supported by the data. Fix the width at
                    % the last accepted value (first pass: the width of
                    % the rise) and fit this number of components again.
                    fwFixed = true;
                    fwNote = sprintf(['the fit drove the IRF width below ' ...
                        'half a delay step; it was fixed at %.3g %s'], ...
                        fwc, C.tu);
                    continue;
                end
                
                % Is every time constant supported by the data?
                relE = F.tauStd ./ F.tau;
                undet = ~isfinite(relE) | relE > 0.5;
                longT = F.tau > 2 * tEnd | (F.tau > 0.5 * tEnd & undet);
                if ~useInf && nnz(longT) == 1 && numel(F.tau) >= 2
                    % A time constant beyond the time range: describe it
                    % as a component that does not decay, if that fits as
                    % well. All later fits then keep that component.
                    Fi = app.autoFitGA(Dv, tv, F.tau(~longT), F.t0, ...
                        F.fwhm, true, irf, fwFixed, caK);
                    if ~isempty(Fi) && Fi.rms <= 1.02 * F.rms
                        F = Fi;  useInf = true;
                        relE = F.tauStd ./ F.tau;
                        undet = ~isfinite(relE) | relE > 0.5;
                        longT = F.tau > 2 * tEnd | (F.tau > 0.5 * tEnd & undet);
                    end
                end
                why = '';
                if any(undet & ~longT)
                    why = 'a time constant is not determined';
                elseif F.hasInf && any(longT)
                    why = ['a time constant beyond the time range besides ' ...
                        'the non-decaying component'];
                elseif nnz(longT) > 1
                    why = 'two time constants beyond the time range';
                elseif any(F.tau < 0.5 * dtZ)
                    why = 'a time constant below half a delay step';
                elseif caK > 0 && any(F.tau < opt.caTauMin * F.fwhm)
                    why = ['a time constant too far below the IRF width ' ...
                        'to be told from the artifact terms'];
                elseif ~isempty(seed) && (F.fwhm > 3 * seed.F.fwhm || ...
                        F.fwhm < seed.F.fwhm / 3 || ...
                        abs(F.t0 - seed.F.t0) > 2 * seed.F.fwhm)
                    why = 'the IRF does not stay at time zero';
                end
                gain = NaN;
                if n > 1, gain = 1 - F.rms / tabRms(n-1); end
                if n > 1 && isempty(why) && gain < opt.minGain
                    why = sprintf('RMS gain below %.0f %%', 100 * opt.minGain);
                end
                if n == 1 || atSeed, why = ''; end   % always kept: one
                                                     % component; the seed fit
                tabN(n) = n;  tabRms(n) = F.rms;  tabGain(n) = gain; %#ok<AGROW>
                tabInf(n) = F.hasInf; %#ok<AGROW>
                tabTau{n} = F.tau;  tabStd{n} = F.tauStd; %#ok<AGROW>
                tabWhy{n} = why;  tabF{n} = F; %#ok<AGROW>
                if ~decided
                    if isempty(why)
                        nStar = n;
                    else
                        stopWhy = why;  decided = true;
                    end
                end
                if decided && ~opt.tryAll, break; end
                prevTau = F.tau;  t0c = F.t0;  fwc = F.fwhm;
                n = n + 1;
            end
            T = struct();
            T.N = tabN;  T.Rms = tabRms;  T.Gain = tabGain;  T.Inf = tabInf;
            T.Tau = tabTau;  T.Std = tabStd;  T.Why = tabWhy;  T.F = tabF;
            T.nStar = nStar;  T.stopWhy = stopWhy;
            T.fwFixed = fwFixed;  T.fwNote = fwNote;  T.caK = caK;
        end
        
        function S = autoSelectModel(app, P, S, n)
            %autoSelectModel  Makes the fit with n components the model of
            %   the summary S (see autoSummarize): spectra, residual, the
            %   quantities quoted in the report, and the oscillations in
            %   its residual. n can be any number in the table of fits,
            %   also one that the automatic choice rejected.
            n = max(1, min(numel(S.tabF), round(n)));
            F = S.tabF{n};
            wl = P.wl(:);
            rows = S.rows;  tv = S.t;  wv = S.wl;
            Dv = P.D(rows, S.cols);
            tEnd = S.tEnd;
            relE = F.tauStd ./ F.tau;
            undet = ~isfinite(relE) | relE > 0.5;
            longT = F.tau > 2 * tEnd | (F.tau > 0.5 * tEnd & undet);
            longNote = '';
            if F.hasInf
                longNote = ['one component does not decay within the time ' ...
                    'range; its lifetime cannot be given'];
            elseif any(longT)
                kL = find(longT, 1);
                if undet(kL)
                    longNote = sprintf(['the time constant %.3g %s lies beyond ' ...
                        'the time range and is not determined'], F.tau(kL), P.tu);
                else
                    longNote = sprintf(['the time constant %.3g %s is longer ' ...
                        'than twice the time range: treat it with caution'], ...
                        F.tau(kL), P.tu);
                end
            end
            nc = numel(F.tau) + double(F.hasInf);
            DADS = NaN(numel(wl), nc);
            DADS(rows, :) = F.A;
            EADS = NaN(numel(wl), nc);
            try
                EADS(rows, :) = app.computeEADSFromDADS(F.A, F.tau, F.hasInf);
            catch
            end
            caSpec = NaN(numel(wl), F.caK);
            caSpec(rows, :) = F.caSpec;
            caPkW = NaN(1, F.caK);  caPkV = NaN(1, F.caK);
            for k = 1:F.caK
                [~, ii] = max(abs(F.caSpec(:, k)));
                caPkW(k) = wv(ii);  caPkV(k) = F.caSpec(ii, k);
            end
            Rv = Dv - F.fit;
            [Ur, ~, Vr] = svd(Rv, 'econ');
            simR = NaN(1, max(0, nc - 1));  simS = simR;
            for k = 1:nc - 1
                e1 = EADS(rows, k);  e2 = EADS(rows, k+1);
                if all(isfinite(e1)) && all(isfinite(e2)) && ...
                        norm(e1) > 0 && norm(e2) > 0
                    simR(k) = (e1.' * e2) / (norm(e1) * norm(e2));
                    simS(k) = (e1.' * e2) / (e1.' * e1);
                end
            end
            pkW = NaN(1, nc);  pkV = NaN(1, nc);
            for k = 1:nc
                [~, ii] = max(abs(F.A(:, k)));
                pkW(k) = wv(ii);  pkV(k) = F.A(ii, k);
            end
            S.model = struct('tau', F.tau, 'tauStd', F.tauStd, ...
                't0', F.t0, 't0Std', F.t0Std, 'fwhm', F.fwhm, ...
                'fwhmStd', F.fwhmStd, 'hasInf', F.hasInf, 'irf', S.irf, ...
                'fwFixed', S.fwFixed, 'irfSeed', F.irfFixed, ...
                'DADS', DADS, 'EADS', EADS, ...
                'rms', F.rms, 'fit', F.fit, 'resid', Rv, ...
                'acT', app.autoAc1(Vr(:, 1)), 'acW', app.autoAc1(Ur(:, 1)), ...
                'simR', simR, 'simS', simS, 'pkW', pkW, 'pkV', pkV, ...
                'longNote', longNote, 'caK', F.caK, 'caSpec', caSpec, ...
                'caPkW', caPkW, 'caPkV', caPkV);
            S.nUse = n;
            
            % ---------- oscillations in the residual ----------
            sgw = P.noiseOut(rows);
            if strcmp(S.irf, 'numerical')
                S.osc = app.autoOscillations(wv, tv, Rv, F.t0, F.fwhm, sgw);
            else
                S.osc = app.autoOscillations(wv, tv, Rv, NaN, NaN, sgw);
            end
        end
        
        function O = autoOscillations(app, wl, t, R, t0, fwhm, sgw)
            %autoOscillations  Periodic signals in a fit residual.
            %   The spectrum (|FFT|^2 summed over the wavelengths, Hann
            %   window) of the residual in the evenly sampled delays
            %   after the instrument response is searched for peaks.
            %   With sgw, the noise of every wavelength, each wavelength
            %   enters the sum divided by the square of its noise
            %   (relative to the median noise): a few noisy channels
            %   would otherwise set the floor of the whole spectrum.
            %   A peak must stand above the local floor (running median)
            %   by 10 times the scatter of the spectrum around that
            %   floor, and be at least twice the floor.
            O = struct('status', 'none', 'msg', '', 'win', [NaN NaN], ...
                'res', NaN, 'fMax', NaN, 'freq', [], 'power', [], ...
                'floor', [], 'pkF', zeros(1, 0), 'pkRatio', zeros(1, 0), ...
                'pkIdx', zeros(1, 0), 'pkWl', {{}}, 'pkNode', zeros(1, 0), ...
                'Pmap', []);
            if ~any(strcmp(app.timeUnit, {'fs', 'ps'}))
                O.status = 'skipped';
                O.msg = 'the time unit is not fs or ps';
                return;
            end
            if ~isfinite(t0)
                O.status = 'skipped';
                O.msg = 'no time zero in the data';
                return;
            end
            t = t(:).';
            jA = find(t >= t0 + 2.5 * fwhm, 1);
            if isempty(jA) || jA > numel(t) - 32
                O.status = 'short';  return;
            end
            d = diff(t(jA:end));
            d0 = median(d(1:min(5, numel(d))));
            k = 1;
            while k <= numel(d) && d(k) <= 1.5 * d0, k = k + 1; end
            jB = jA + k - 1;      % last delay of the evenly sampled part
            if jB - jA + 1 < 32
                O.status = 'short';  return;
            end
            O.win = [t(jA), t(jB)];
            try
                [Pm, fr, ci] = app.computeCoherence(R, t, struct( ...
                    'tMin', t(jA), 'tMax', t(jB), 'apod', 'hann', ...
                    'zeroPad', 4, 'detrend', 'linear', ...
                    'freqUnit', 'cm-1', 'normMode', 'none'));
            catch errO
                O.status = 'failed';  O.msg = errO.message;  return;
            end
            fr = fr(:).';
            wgt = ones(size(Pm, 1), 1);
            if nargin >= 7 && numel(sgw) == size(Pm, 1)
                sgw = sgw(:);
                okw = isfinite(sgw) & sgw > 0;
                if any(okw)
                    wgt(okw)  = (median(sgw(okw)) ./ sgw(okw)).^2;
                    wgt(~okw) = 0;
                end
            end
            Ps = sum(Pm .* (wgt * ones(1, size(Pm, 2))), 1);
            res = 2 * fr(end) / ci.Nwin;       % 1 / (length of the window)
            df  = fr(2) - fr(1);
            hw  = max(10, round(10 * res / df));
            flr = app.autoRunMed(Ps, 2 * hw + 1).';
            nF  = numel(fr);
            isMax = false(1, nF);
            isMax(2:nF-1) = Ps(2:nF-1) >= Ps(1:nF-2) & Ps(2:nF-1) > Ps(3:nF);
            inR = fr >= 2 * res;
            dev = Ps(inR) - flr(inR);
            spr = 1.4826 * median(abs(dev - median(dev)));
            cand = find(isMax & inR & Ps > flr + 10 * spr & Ps > 2 * flr);
            [~, ord] = sort(Ps(cand), 'descend');
            cand = cand(ord);
            keep = zeros(1, 0);
            for q = 1:numel(cand)
                if isempty(keep) || all(abs(fr(cand(q)) - fr(keep)) >= 2 * res)
                    keep(end+1) = cand(q); %#ok<AGROW>
                end
                if numel(keep) >= 6, break; end
            end
            keep = sort(keep);
            O.res = res;  O.fMax = fr(end);
            O.freq = fr;  O.power = Ps;  O.floor = flr;  O.Pmap = Pm;
            O.pkIdx = keep;
            O.pkF = fr(keep);
            O.pkRatio = Ps(keep) ./ flr(keep);
            O.pkWl = cell(1, numel(keep));
            O.pkNode = NaN(1, numel(keep));
            for q = 1:numel(keep)
                % Where along the spectrum the oscillation is strongest;
                % two maxima with a deep minimum between them point to a
                % node (phase jump) at that minimum.
                pr = Pm(:, keep(q));
                pr = (pr + [pr(1); pr(1:end-1)] + [pr(2:end); pr(end)]) / 3;
                [m1, i1] = max(pr);
                bestJ = 0;  bestV = 0;  bestNode = NaN;
                for j = 2:numel(pr) - 1
                    if j == i1 || pr(j) < pr(j-1) || pr(j) < pr(j+1), continue; end
                    if pr(j) < 0.25 * m1, continue; end
                    a = min(i1, j);  b = max(i1, j);
                    [dip, id] = min(pr(a:b));
                    if dip < 0.3 * pr(j) && pr(j) > bestV
                        bestJ = j;  bestV = pr(j);  bestNode = wl(a + id - 1);
                    end
                end
                if bestJ > 0
                    O.pkWl{q} = sort([wl(i1), wl(bestJ)]);
                    O.pkNode(q) = bestNode;
                else
                    O.pkWl{q} = wl(i1);
                end
            end
            if isempty(keep), O.status = 'nopeak'; else, O.status = 'found'; end
        end

        function L = autoReportLines(app, P, S, src)
            %autoReportLines  Text of the automatic-analysis report
            %   (cell array of lines, plain ASCII).
            if nargin < 4 || isempty(src), src = 'current data'; end
            tu = app.tUnitHdr();  xu = app.xUnitTxt();
            L = {};
            L{end+1} = 'AUTOMATIC ANALYSIS';
            L{end+1} = sprintf('Data: %s', src);
            L{end+1} = sprintf('      %d wavelengths (%.4g ... %.4g %s) x %d delays (%.4g ... %.4g %s)', ...
                numel(P.wlIn), min(P.wlIn), max(P.wlIn), xu, ...
                numel(P.tIn), min(P.tIn), max(P.tIn), tu);
            L{end+1} = '';
            L{end+1} = '1. CORRECTIONS (proposal)';
            ok = P.rowOK;
            L{end+1} = sprintf('  Noise        median %.3g, from %.3g to %.3g along the spectrum', ...
                median(P.noise(ok)), min(P.noise(ok)), max(P.noise(ok)));
            if P.hasT0
                L{end+1} = sprintf('  Signal       starts at t = %.4g %s (delay no. %d)', ...
                    P.tOn, tu, P.jOn);
            else
                L{end+1} = sprintf('  Signal       no time zero found: %s', P.t0Why);
            end
            switch P.bgStatus
                case 'applied'
                    [bm, ib] = max(abs(P.bg));
                    L{end+1} = sprintf(['  Background   mean of the first %d delays ' ...
                        '(t <= %.4g %s) subtracted'], P.nBg, P.tIn(P.nBg), tu);
                    L{end+1} = sprintf(['               largest value %.3g at %.4g %s; ' ...
                        'on average %.1f times its own uncertainty'], ...
                        bm, P.wlIn(ib), xu, P.bgSig);
                case 'zero'
                    L{end+1} = sprintf(['  Background   none: the mean of the first %d ' ...
                        'delays is zero within the noise'], P.nBg);
                case 'few'
                    L{end+1} = '  Background   not set: too few delays before the signal';
                case 'signal'
                    L{end+1} = '  Background   not set: the first delays look like signal';
                case 'notime0'
                    L{end+1} = '  Background   not set: no time zero';
                otherwise
                    L{end+1} = '  Background   switched off';
            end
            switch P.maskStatus
                case 'applied'
                    if isempty(P.maskReg)
                        L{end+1} = '  Mask         none inside the spectrum (noisy end: see Crop)';
                    end
                    for q = 1:size(P.maskReg, 1)
                        if isfinite(P.opt.maskSN)
                            L{end+1} = sprintf(['  Mask         %.4g ... %.4g %s ' ...
                                '(noise up to %.1f times the neighbours, signal ' ...
                                'down to %.0f times the noise)'], ...
                                P.maskReg(q, 1), P.maskReg(q, 2), xu, ...
                                P.maskReg(q, 3), P.maskReg(q, 4)); %#ok<AGROW>
                        else
                            L{end+1} = sprintf(['  Mask         %.4g ... %.4g %s ' ...
                                '(noise up to %.1f times the neighbours)'], ...
                                P.maskReg(q, 1), P.maskReg(q, 2), xu, ...
                                P.maskReg(q, 3)); %#ok<AGROW>
                        end
                    end
                case 'none'
                    if isempty(P.keptReg)
                        L{end+1} = sprintf(['  Mask         none: no channel is more ' ...
                            'than %.3g times noisier than its neighbours'], ...
                            P.opt.maskK);
                    else
                        L{end+1} = ['  Mask         none: the noisy channels show a ' ...
                            'clear signal'];
                    end
                case 'toomany'
                    L{end+1} = '  Mask         not set: too many noisy channels';
                otherwise
                    L{end+1} = '  Mask         switched off';
            end
            if ~isempty(P.keptReg)
                L{end+1} = '               noisier than the neighbours, but kept:';
                nK = size(P.keptReg, 1);
                for q = 1:min(nK, 6)
                    L{end+1} = sprintf(['                 %.4g ... %.4g %s  (noise up to ' ...
                        '%.1f times, signal at least %.0f times the noise)'], ...
                        P.keptReg(q, 1), P.keptReg(q, 2), xu, ...
                        P.keptReg(q, 3), P.keptReg(q, 4)); %#ok<AGROW>
                end
                if nK > 6
                    L{end+1} = sprintf('                 ... and %d more ranges', nK - 6);
                end
            end
            ch = P.chirp;
            vr = P.rowOK & ~P.maskRows;
            switch ch.status
                case {'applied', 'found'}
                    wv = P.wlIn(vr);  tv0 = ch.t0(vr);
                    L{end+1} = sprintf(['  Chirp        time zero runs from %.3g %s at ' ...
                        '%.4g %s to %.3g %s at %.4g %s (span %.3g %s)'], ...
                        tv0(1), tu, wv(1), xu, tv0(end), tu, wv(end), xu, ch.span, tu);
                    if strcmp(ch.model, 'parabola')
                        L{end+1} = '               (parabola: the chirp formula gave no sensible curve)';
                    end
                    L{end+1} = sprintf(['               curve through the rise of %d wavelengths (RMS %.2g %s); ' ...
                        '%.0f %% of all rises within %.2g %s of it'], ...
                        size(ch.pts, 1), ch.rms, tu, 100 * ch.fracIn, ch.tol, tu);
                    if strcmp(ch.status, 'applied')
                        L{end+1} = '               corrected: time zero is now 0 at every wavelength';
                    else
                        L{end+1} = '               NOT corrected (switched off)';
                    end
                case 'none'
                    L{end+1} = sprintf(['  Chirp        none: time zero varies by %.2g %s along ' ...
                        'the spectrum (negligible below %.2g %s)'], ...
                        ch.span, tu, 0.75 * ch.tol, tu);
                    if ch.shift ~= 0
                        L{end+1} = sprintf(['               time axis shifted by %.4g %s so that ' ...
                            'time zero is 0'], -ch.shift, tu);
                    end
                case 'failed'
                    L{end+1} = sprintf('  Chirp        not determined: %s', ch.msg);
                case 'notime0'
                    L{end+1} = '  Chirp        not determined: no time zero';
                otherwise
                    L{end+1} = '  Chirp        switched off';
            end
            if isempty(P.cropT) && isempty(P.cropW)
                if P.opt.doCrop
                    L{end+1} = '  Crop         nothing to remove';
                else
                    L{end+1} = '  Crop         switched off';
                end
            else
                if ~isempty(P.cropT)
                    L{end+1} = sprintf(['  Crop         t = %.4g ... %.4g %s (%d delays without ' ...
                        'data at all wavelengths removed)'], ...
                        P.cropT(1), P.cropT(2), tu, nnz(~P.keepT));
                end
                if ~isempty(P.cropW)
                    L{end+1} = sprintf(['  Crop         %.4g ... %.4g %s (%d noisy channels at ' ...
                        'the end of the spectrum removed)'], ...
                        P.cropW(1), P.cropW(2), xu, nnz(~P.keepW));
                end
            end
            L{end+1} = sprintf('  Result       %d wavelengths x %d delays', ...
                numel(P.wl), numel(P.t));
            if isempty(S), return; end

            L{end+1} = '';
            L{end+1} = '2. NUMBER OF COMPONENTS';
            kS = min(8, numel(S.sv));
            parts = cell(1, kS);
            for k = 1:kS
                parts{k} = sprintf('%.2g', S.sv(k) / S.sv(1));
            end
            L{end+1} = sprintf('  Singular values / first:  %s', strjoin(parts, '  '));
            L{end+1} = sprintf(['  above 5 x noise floor: %d;  with smooth vectors ' ...
                '(spectrum and time): %d'], S.nAbove, S.nStruct);
            switch S.ca.mode
                case 'off'
                    L{end+1} = '  Coherent artifact terms: switched off';
                otherwise
                    if ~isempty(S.ca.note)
                        L{end+1} = sprintf('  Coherent artifact terms: none (%s)', S.ca.note);
                    elseif S.ca.use > 0 && strcmp(S.ca.why, 'on')
                        L{end+1} = sprintf(['  Coherent artifact terms: %d, switched on ' ...
                            '(RMS without them: %.4g with %d component(s))'], ...
                            S.ca.use, S.ca.rms0, S.ca.n0);
                    elseif S.ca.use > 0 && strcmp(S.ca.why, 'fast')
                        L{end+1} = sprintf(['  Coherent artifact terms: %d (Gaussian of the IRF ' ...
                            'and its two derivatives).'], S.ca.use);
                        L{end+1} = sprintf(['    Without them the fits need a time constant of ' ...
                            '%.3g %s, below a third of the IRF width'], S.ca.tauFast, tu);
                        L{end+1} = sprintf(['    (%.3g %s). A signal that short cannot be told ' ...
                            'from a coherent artifact;'], S.ca.fw0, tu);
                        L{end+1} = ['    it is described by the terms and not counted as a ' ...
                            'component.'];
                        L{end+1} = sprintf(['    (RMS without the terms: %.4g with %d ' ...
                            'component(s).)'], S.ca.rms0, S.ca.n0);
                    elseif S.ca.use > 0
                        L{end+1} = sprintf(['  Coherent artifact terms: %d (Gaussian of the IRF ' ...
                            'and its two derivatives). They lower the RMS by %.1f %%;'], ...
                            S.ca.use, 100 * S.ca.gain);
                        L{end+1} = '    the components below were counted with them in place.';
                    else
                        L{end+1} = sprintf(['  Coherent artifact terms: not needed (three terms ' ...
                            'lower the RMS by %.1f %%, less than %.0f %%)'], ...
                            100 * S.ca.gain, 100 * S.caGainMin);
                    end
            end
            if S.ca.use > 0 && isfinite(S.ca.fw0)
                L{end+1} = sprintf(['    IRF FWHM: %.4g %s without the terms, %.4g %s ' ...
                    'with them.'], S.ca.fw0, tu, S.tabF{S.nStar}.fwhm, tu);
            end
            L{end+1} = '  Global fits with N components (time constants and IRF free):';
            L{end+1} = sprintf('    N   RMS        gain    time constants (%s)', tu);
            for n = 1:numel(S.tabN)
                tp = cell(1, numel(S.tabTau{n}));
                for k = 1:numel(tp)
                    tp{k} = sprintf('%.3g', S.tabTau{n}(k));
                end
                if S.tabInf(n), tp{end+1} = 'non-decaying'; end
                if n == 1
                    gs = '   -  ';
                else
                    gs = sprintf('%5.1f%%', 100 * S.tabGain(n));
                end
                mark = '';
                if ~isempty(S.tabWhy{n})
                    mark = ['   (' S.tabWhy{n} ')'];
                end
                if n == S.nStar && n == S.nUse
                    mark = '   <- chosen';
                elseif n == S.nStar
                    mark = '   <- automatic choice';
                elseif n == S.nUse
                    mark = [mark '   <- used (set by hand)'];
                end
                L{end+1} = sprintf('    %d   %-9.4g  %s  %s%s', n, S.tabRms(n), ...
                    gs, strjoin(tp, ', '), mark); %#ok<AGROW>
            end
            if S.ca.use > 0 && S.ca.nSeed > 1
                L{end+1} = sprintf(['  (N below %d: t0 and FWHM kept at the values of the ' ...
                    'fit with %d, for comparison.)'], S.ca.nSeed, S.ca.nSeed);
            end
            if S.hitMax
                L{end+1} = sprintf(['  Automatic choice: %d, the largest number tried. A ' ...
                    'further component may still improve the fit.'], S.nStar);
            else
                L{end+1} = sprintf('  Automatic choice: %d. One more was rejected: %s.', ...
                    S.nStar, S.stopWhy);
            end
            if S.nUse ~= S.nStar
                L{end+1} = sprintf(['  Used below: %d components, set by hand ' ...
                    '("Components used").'], S.nUse);
            end

            Mo = S.model;
            L{end+1} = '';
            if S.nUse == S.nStar
                L{end+1} = '3. CHOSEN FIT';
            else
                L{end+1} = sprintf('3. FIT WITH %d COMPONENTS (set by hand)', S.nUse);
            end
            nT = numel(Mo.tau);
            for k = 1:nT
                if isfinite(Mo.tauStd(k))
                    es = sprintf('+/- %.2g', Mo.tauStd(k));
                else
                    es = '[not determined]';
                end
                L{end+1} = sprintf(['  tau_%d = %-9.4g %-16s %s   DADS largest at ' ...
                    '%.4g %s (%s)'], k, Mo.tau(k), es, tu, Mo.pkW(k), xu, ...
                    app.pick(Mo.pkV(k) >= 0, 'positive', 'negative')); %#ok<AGROW>
            end
            if Mo.hasInf
                L{end+1} = sprintf(['  non-decaying component              DADS largest at ' ...
                    '%.4g %s (%s)'], Mo.pkW(end), xu, ...
                    app.pick(Mo.pkV(end) >= 0, 'positive', 'negative'));
            end
            if ~isempty(Mo.longNote)
                L{end+1} = sprintf('  Note: %s.', Mo.longNote);
            end
            if strcmp(Mo.irf, 'none')
                L{end+1} = '  IRF: none (pure exponentials; no time zero in the data)';
            elseif Mo.irfSeed
                L{end+1} = sprintf('  IRF t0   = %.3g %s (fixed)', Mo.t0, tu);
                L{end+1} = sprintf('  IRF FWHM = %.4g %s (fixed)', Mo.fwhm, tu);
                L{end+1} = sprintf(['  Note: t0 and FWHM are those of the fit with %d ' ...
                    'components. In a fit with too few'], S.ca.nSeed);
                L{end+1} = ['        components, free artifact terms leave time zero ' ...
                    'and stand in for the missing kinetics.'];
            else
                if isfinite(Mo.t0Std)
                    L{end+1} = sprintf('  IRF t0   = %.3g +/- %.2g %s', Mo.t0, Mo.t0Std, tu);
                else
                    L{end+1} = sprintf('  IRF t0   = %.3g %s', Mo.t0, tu);
                end
                if Mo.fwFixed
                    L{end+1} = sprintf('  IRF FWHM = %.4g %s (fixed)', Mo.fwhm, tu);
                    L{end+1} = sprintf('  Note: %s.', S.fwNote);
                elseif isfinite(Mo.fwhmStd)
                    L{end+1} = sprintf('  IRF FWHM = %.4g +/- %.2g %s', Mo.fwhm, Mo.fwhmStd, tu);
                else
                    L{end+1} = sprintf('  IRF FWHM = %.4g %s', Mo.fwhm, tu);
                end
            end
            if Mo.caK > 0
                caNm = {'Gaussian', '1st derivative', '2nd derivative'};
                L{end+1} = sprintf(['  Coherent artifact: %d terms of the IRF width, ' ...
                    'spectra free'], Mo.caK);
                for k = 1:Mo.caK
                    L{end+1} = sprintf('    %-14s largest %.3g at %.4g %s', ...
                        caNm{k}, Mo.caPkV(k), Mo.caPkW(k), xu); %#ok<AGROW>
                end
            end
            L{end+1} = sprintf('  RMS = %.4g   (noise estimate %.3g)', Mo.rms, S.rmsNoise);
            if Mo.acT > 0.5 && Mo.acW > 0.5
                L{end+1} = sprintf(['  Residual: structured (autocorrelation of its main ' ...
                    'part: %.2f in time, %.2f along the spectrum).'], Mo.acT, Mo.acW);
                L{end+1} = '            Not everything is described by the exponentials.';
            else
                L{end+1} = sprintf(['  Residual: no clear structure (autocorrelation %.2f ' ...
                    'in time, %.2f along the spectrum).'], Mo.acT, Mo.acW);
            end
            if numel(Mo.simR) >= 1
                L{end+1} = '  EADS (sequential picture 1 -> 2 -> ...), consecutive spectra:';
                for k = 1:numel(Mo.simR)
                    if ~isfinite(Mo.simR(k)), continue; end
                    if Mo.simR(k) > 0.95
                        hint = 'same shape';
                    elseif Mo.simR(k) > 0.8
                        hint = 'similar shape';
                    else
                        hint = 'different shape';
                    end
                    L{end+1} = sprintf(['    EADS %d -> %d: shape correlation %.2f (%s), ' ...
                        'size ratio %.2f'], k, k + 1, Mo.simR(k), hint, Mo.simS(k)); %#ok<AGROW>
                end
            end

            O = S.osc;
            L{end+1} = '';
            L{end+1} = '4. OSCILLATIONS IN THE RESIDUAL';
            switch O.status
                case {'found', 'nopeak'}
                    L{end+1} = sprintf(['  Window %.4g ... %.4g %s; detectable %.0f ... %.0f cm-1 ' ...
                        '(resolution about %.0f cm-1)'], O.win(1), O.win(2), tu, ...
                        2 * O.res, O.fMax, O.res);
                    if strcmp(O.status, 'nopeak')
                        L{end+1} = '  No significant peak.';
                    end
                    for q = 1:numel(O.pkF)
                        ws = O.pkWl{q};
                        if numel(ws) == 2
                            wtxt = sprintf(['strongest at %.4g and %.4g %s, ' ...
                                'node near %.4g %s'], ws(1), ws(2), xu, O.pkNode(q), xu);
                        else
                            wtxt = sprintf('strongest at %.4g %s', ws(1), xu);
                        end
                        L{end+1} = sprintf('  %6.0f cm-1   %5.1f x floor   %s', ...
                            O.pkF(q), O.pkRatio(q), wtxt); %#ok<AGROW>
                    end
                case 'short'
                    L{end+1} = '  Not analysed: fewer than 32 evenly spaced delays after time zero.';
                case 'skipped'
                    L{end+1} = sprintf('  Not analysed: %s.', O.msg);
                otherwise
                    L{end+1} = sprintf('  Not analysed: %s', O.msg);
            end

            % ---------- things to check ----------
            C = P.checks;
            if S.nUse ~= S.nStar && ~isempty(S.tabWhy{S.nUse})
                C{end+1} = sprintf(['The fit with %d components was set by hand. ' ...
                    'The automatic choice rejected it: %s.'], S.nUse, S.tabWhy{S.nUse});
            end
            if Mo.caK > 0
                C{end+1} = ['Artifact terms absorb every signal confined to time ' ...
                    'zero: a real artifact, a leftover chirp, and processes ' ...
                    'faster than the IRF. Time constants below the IRF width ' ...
                    'are not reliable next to them.'];
                fwA = S.tabF{S.nStar}.fwhm;
                if isfinite(S.ca.fw0) && (fwA > 1.3 * S.ca.fw0 || ...
                        fwA < S.ca.fw0 / 1.3)
                    C{end+1} = sprintf(['The artifact terms change the IRF FWHM from ' ...
                        '%.3g to %.3g %s. With an artifact in the data the ' ...
                        'second value is the better one. If the dashed ' ...
                        'spectra on the Components tab do not look like an ' ...
                        'artifact, run again with "Coherent artifact: ' ...
                        'Off" and compare.'], S.ca.fw0, fwA, tu);
                end
                if any(strcmp(P.chirp.status, {'applied', 'found'}))
                    C{end+1} = ['With an artifact in the data the steepest rise ' ...
                        'is not exactly time zero: where an artifact of one ' ...
                        'sign is as large as the signal it lies up to about ' ...
                        '0.3 IRF widths off, and the chirp curve of this run ' ...
                        'with it. If a scan of the solvent or substrate alone ' ...
                        'exists, take the chirp from it (Chirp Corr window).'];
                end
            end
            if ~isfinite(P.jOn)
                C{end+1} = ['The data do not change with time beyond the noise. ' ...
                    'The fit below then describes noise or a constant signal.'];
            end
            if S.nStruct > S.nUse
                % (counted against the fit in use; with artifact terms in
                % that fit, one of the smooth vectors may be the artifact)
                if Mo.caK > 0
                    extra = ['the coherent artifact (described by the ' ...
                        'artifact terms), oscillations'];
                else
                    extra = 'oscillations';
                end
                C{end+1} = sprintf(['%d singular vectors are smooth but %d ' ...
                    'component(s) were fitted. The extra ones can be ' ...
                    '%s, a time zero that still varies along the ' ...
                    'spectrum, or kinetics that are not exponential.'], ...
                    S.nStruct, S.nUse, extra);
            end
            if S.hitMax
                C{end+1} = ['The number of components stopped at the largest ' ...
                    'number tried, not because a further one was rejected.'];
            end
            if Mo.fwFixed
                C{end+1} = ['The IRF width could not be fitted and was fixed. ' ...
                    'Time constants close to it depend on that value.'];
            end
            if strcmp(Mo.irf, 'numerical') && Mo.caK > 0 && ...
                    any(Mo.tau < 2 * Mo.fwhm)
                C{end+1} = sprintf(['A time constant (%.3g %s) is close to the IRF ' ...
                    'width (%.3g %s) while artifact terms are in use. The ' ...
                    'two are hard to separate: it may describe what the ' ...
                    'terms leave of the artifact or of a leftover chirp. ' ...
                    'Compare with the fit with one component less ' ...
                    '("Components used").'], min(Mo.tau), tu, Mo.fwhm, tu);
            elseif strcmp(Mo.irf, 'numerical') && any(Mo.tau < Mo.fwhm)
                C{end+1} = sprintf(['A time constant (%.3g %s) is shorter than the ' ...
                    'IRF width (%.3g %s): it depends on the IRF shape assumed ' ...
                    '(Gaussian), and it can also be a coherent artifact at ' ...
                    'time zero instead of a population.'], ...
                    min(Mo.tau), tu, Mo.fwhm, tu);
            end
            L{end+1} = '';
            L{end+1} = '5. TO CHECK BY EYE';
            if isempty(C)
                L{end+1} = '  Nothing was flagged.';
            end
            for q = 1:numel(C)
                L{end+1} = sprintf('  - %s', C{q}); %#ok<AGROW>
            end
            L{end+1} = '';
            L{end+1} = 'The time constants describe the data; which states and steps';
            L{end+1} = 'they belong to is not decided here (see Target Analysis).';
        end
        
        % =================================================================
        % Drawing
        % =================================================================
        % ---- Spectral-axis unit helpers (mirror the time-unit ones) ----
        function s = xAxLabel(app)
            % Full axis label for the spectral axis (TeX interpreter).
            if strcmp(app.specUnit, 'cm-1')
                s = 'Wavenumber (cm^{-1})';
            else
                s = 'Wavelength (nm)';
            end
        end
        function s = xUnitAx(app)
            % Unit string for titles / legends (TeX interpreter).
            if strcmp(app.specUnit, 'cm-1'), s = 'cm^{-1}'; else, s = 'nm'; end
        end
        function s = xUnitTxt(app)
            % Unit string for plain-text UI labels (no interpreter).
            if strcmp(app.specUnit, 'cm-1'), s = 'cm-1'; else, s = 'nm'; end
        end
        function s = xSym(app)
            % Variable symbol for the spectral coordinate (TeX): lambda
            % for wavelength, nu-tilde for wavenumber.
            if strcmp(app.specUnit, 'cm-1'), s = '\nu'; else, s = '\lambda'; end
        end
        function s = xSymTxt(app)
            % The same symbol for plain-text UI labels (no interpreter).
            s = app.uiStr(app.xSym());
        end
        function s = uiStr(app, s)
            %UISTR  TeX-style markup -> plain Unicode text.
            %   For UI components without an interpreter (uilabel,
            %   uibutton, uicheckbox, dropdown items, tooltips, panel
            %   titles, alerts, text areas), where '\lambda' would be
            %   shown letter by letter. Axes titles / labels / legends
            %   are rendered by the TeX interpreter and must NOT go
            %   through this function.
            %   Accepts a char vector or a cell array of char vectors.
            %   The characters are built with char(), so this source
            %   file stays plain ASCII whatever its encoding.
            %   Do not pass file paths: a Windows path such as
            %   C:\data\nu1 contains "\nu".
            if iscell(s)
                for k = 1:numel(s), s{k} = app.uiStr(s{k}); end
                return;
            end
            if ~ischar(s) || isempty(s), return; end
            % Single-digit subscript after a letter: t_0, \tau_1, \omega_0
            for d = 0:9
                s = regexprep(s, sprintf('([A-Za-z])_%d(?![0-9])', d), ...
                    ['$1' char(8320 + d)]);
            end
            s = strrep(s, '^{-1}', '-1');                 % cm^{-1} -> cm-1
            s = regexprep(s, '\^2(?![0-9])', char(178));  % ^2 -> superscript 2
            s = regexprep(s, '([_^])\{([^{}]*)\}', '$1$2'); % t_{min} -> t_min
            map = { ...
                '\rightarrow', 8594; '\lambda', 955; '\sigma', 963; ...
                '\alpha', 945;  '\omega', 969;  '\Delta', 916; ...
                '\infty', 8734; '\beta', 946;   '\cdot', 183; ...
                '\tau', 964;    '\nu', 957;     '\mu', 956; ...
                '\pi', 960};
            for k = 1:size(map, 1)
                s = strrep(s, map{k, 1}, char(map{k, 2}));
            end
        end
        function s = undetNote(~, info, fld, k)
            % Result-text note for a parameter whose uncertainty came
            % back as Inf (see fitParamUncertainty).
            s = '';
            if isstruct(info) && isfield(info, fld)
                v = info.(fld);
                if numel(v) >= k && isinf(v(k))
                    s = '  [not determined by the data]';
                end
            end
        end
        function s = pick(~, cond, a, b)
            % Inline conditional for labels: pick(cond, a, b).
            if cond, s = a; else, s = b; end
        end
        
        function s = xHdr(app)
            % Column header for CSV exports.
            if strcmp(app.specUnit, 'cm-1'), s = 'wavenumber_cm-1'; else, s = 'wavelength_nm'; end
        end
        
        function s = tUnitAx(app)
            % Time-unit string for axes labels/titles (TeX interpreter).
            switch app.timeUnit
                case 'us', s = '\mus';
                case 'ns', s = 'ns';
                case 'fs', s = 'fs';
                otherwise, s = 'ps';
            end
        end
        function s = tUnitTxt(app)
            % Time-unit string for plain-text UI labels (no interpreter).
            switch app.timeUnit
                case 'us', s = [char(956) 's'];     % micro sign + s
                case 'ns', s = 'ns';
                case 'fs', s = 'fs';
                otherwise, s = 'ps';
            end
        end
        function s = tUnitHdr(app)
            % Time-unit string for filenames / CSV headers (ASCII-safe).
            switch app.timeUnit
                case 'us', s = 'us';
                case 'ns', s = 'ns';
                case 'fs', s = 'fs';
                otherwise, s = 'ps';
            end
        end
        
        function v = defaultChirpTMax(app)
            % Sensible default upper bound for the chirp-correction
            % display window, expressed in the same numeric units as
            % app.delay (which carries app.timeUnit).
            % Roughly 2 ps regardless of which unit the data uses:
            %   fs -> 500   (we use 500 because few-fs steps span ~ a ps;
            %                500 fs gives a clean view of the IRF region)
            %   ps -> 2
            %   ns -> 0.002
            %   us -> 0.000002
            switch app.timeUnit
                case 'fs', v = 500;
                case 'ns', v = 0.002;
                case 'us', v = 2e-6;
                otherwise, v = 2;       % ps
            end
        end
        
        function onTimeUnitChange(app, newUnit)
            if strcmp(app.timeUnit, newUnit), return; end
            app.timeUnit = newUnit;
            % Refresh the chirp-correction display window upper bound to
            % a sensible default for the new unit (e.g. 500 when fs,
            % 2 when ps). We only override when no data is loaded yet,
            % or when the previous value was the previous unit's
            % default and not a custom user value.  Here we just always
            % refresh -- if the user already customized, they can
            % re-enter their value in the chirp window.
            if ~isempty(app.delay)
                app.chirpViewTMax = min(app.defaultChirpTMax(), max(app.delay));
            else
                app.chirpViewTMax = app.defaultChirpTMax();
            end
            % Close any open child windows -- their static labels would
            % otherwise be inconsistent with the new unit. They will use
            % the new unit next time they are opened.
            app.closeChildWindows();
            % Results depend on the unit: coherence / time-frequency
            % wavenumbers are derived from it, and stored report text
            % quotes it. Drop them rather than restore them mislabelled.
            app.invalidateChildResults();
            app.refreshTimeUnitLabels();
            app.updateAll();
            app.savePrefs();
        end
        
        function onSpecUnitChange(app, newUnit)
            if strcmp(app.specUnit, newUnit), return; end
            app.specUnit = newUnit;
            % Child windows bake the unit into their labels, titles and
            % table headers when they open; close them so they come back
            % relabelled. Same policy as the time-unit switch.
            app.closeChildWindows();
            app.invalidateChildResults();   % stored labels quote the unit
            app.refreshSpecUnitLabels();
            app.updateAll();
            app.savePrefs();
        end
        
        function refreshSpecUnitLabels(app)
            % Toolbar labels carrying the spectral unit. Axes labels and
            % titles are re-applied inside the draw functions, which read
            % app.xAxLabel() / app.xUnitAx() fresh every call.
            if ~isempty(app.lblSelWL) && isvalid(app.lblSelWL)
                app.lblSelWL.Text = sprintf('%s (%s):', app.xSymTxt(), app.xUnitTxt());
            end
        end
        
        function refreshTimeUnitLabels(app)
            % Re-set Text on all stored toolbar labels that carry time
            % units.  Axes labels/titles are re-applied inside drawMap2D /
            % drawSpectrum / drawKinetics, which read app.tUnitAx() fresh
            % every call.
            u = app.tUnitTxt();
            if ~isempty(app.lblSplitThreshold) && isvalid(app.lblSplitThreshold)
                app.lblSplitThreshold.Text = sprintf('Split thr. (%s):', u);
            end
            if ~isempty(app.lblViewTMin) && isvalid(app.lblViewTMin)
                app.lblViewTMin.Text = sprintf('Display t_min (%s):', u);
            end
            if ~isempty(app.lblViewTMax) && isvalid(app.lblViewTMax)
                app.lblViewTMax.Text = sprintf('t_max (%s):', u);
            end
            if ~isempty(app.lblSelT) && isvalid(app.lblSelT)
                app.lblSelT.Text = sprintf('Delay (%s):', u);
            end
        end
        
        function pos = fitFigurePosition(~, desiredW, desiredH)
            % Return a [x y w h] Position for a uifigure that fits
            % inside the current screen. Caps the size to the screen
            % minus an OS-chrome margin (title bar, taskbar) and centers
            % it. Falls back gracefully if ScreenSize is unavailable.
            try
                ss = get(groot, 'ScreenSize');     % [x y w h]
                maxW = ss(3) - 40;
                maxH = ss(4) - 100;                 % room for taskbar
                figW = min(desiredW, maxW);
                figH = min(desiredH, maxH);
                figX = max(20, round((ss(3) - figW) / 2));
                figY = max(40, round((ss(4) - figH) / 2));
            catch
                figW = desiredW;  figH = desiredH;  figX = 80;  figY = 60;
            end
            pos = [figX, figY, figW, figH];
        end
        
        function applyLightTheme(~, fig)
            % Force the given uifigure to the built-in "light" theme so
            % the app looks consistent regardless of the OS theme the
            % user is running (this matters most for compiled apps that
            % would otherwise inherit Windows dark mode at runtime).
            % The Theme property and matlab.graphics.theme helpers were
            % introduced in R2025a; older releases don't have them, so
            % we fall back silently.
            if isempty(fig) || ~isvalid(fig), return; end
            try
                fig.Theme = 'light';
            catch
                % Theme property not supported on this release; ignore.
            end
        end
        
        function updateAll(app)
            app.recompute();
            app.drawMap2D();
            drawnow limitrate;
            app.drawSpectrum();
            app.drawKinetics();
        end
        
        function [lo, hi] = clipViewToMask(app, tSub, scaleStr)
            % Clip user's mainView[Tmin,Tmax] to data in tSub and scale limits.
            tmin = app.mainViewTMin;
            tmax = app.mainViewTMax;
            useView = ~isempty(tmin) && ~isempty(tmax) && tmin < tmax;
            dataLo = min(tSub);  dataHi = max(tSub);
            if useView
                lo = max(dataLo, tmin);
                hi = min(dataHi, tmax);
            else
                lo = dataLo;  hi = dataHi;
            end
            if strcmp(scaleStr, 'log')
                pos = tSub(tSub > 0);
                if ~isempty(pos), minPos = min(pos); else, minPos = eps; end
                lo = max(lo, minPos);
            end
            if hi <= lo
                lo = dataLo;  hi = dataHi;
                if strcmp(scaleStr, 'log')
                    pos = tSub(tSub > 0);
                    if ~isempty(pos), lo = min(pos); end
                end
            end
            if hi <= lo
                hi = lo + eps(lo)*100;
                if hi <= lo, hi = lo + 0.1; end
            end
        end
        
        function drawMap2D(app)
            if isempty(app.deltaA), return; end
            
            t   = app.delay;
            thr = app.splitThreshold;
            switch app.delayScaleMode
                case 'linear'
                    masks  = {true(size(t))};
                    scales = {'linear'};
                case 'log'
                    pos = t > 0;
                    masks  = {pos};
                    scales = {'log'};
                case 'split'
                    maskTop = t >= thr;
                    maskBot = t <= thr;
                    masks  = {maskTop, maskBot};
                    scales = {'log',   'linear'};
            end
            
            [cMin, cMax] = app.computeZLim();
            cmap = app.getColormap();
            
            app.hVLine = gobjects(0);
            app.hHLine = gobjects(0);
            for k = 1:numel(app.ax2D)
                ax = app.ax2D(k);
                cla(ax);  hold(ax, 'off');
                tSub = t(masks{k});
                if isempty(tSub), continue; end
                dSub = app.deltaA(:, masks{k});
                h = app.plotHeatmapOn(ax, tSub, dSub);
                h.ButtonDownFcn  = @(~,ev) app.on2DClick(ev);
                ax.ButtonDownFcn = @(~,ev) app.on2DClick(ev);
                ax.CLim  = [cMin, cMax];
                colormap(ax, cmap);
                set(ax, 'YScale', scales{k});
                
                % Apply user display range
                [ylo, yhi] = app.clipViewToMask(tSub, scales{k});
                ylim(ax, [ylo, yhi]);
                
                if k == 1 && isempty(app.cbarAx)
                    cb = colorbar(ax);  cb.Label.String = '\DeltaA or Intensity';
                end
                
                if k == numel(app.ax2D)
                    xlabel(ax, app.xAxLabel());
                else
                    xlabel(ax, '');
                    ax.XTickLabel = [];
                end
                ylabel(ax, sprintf('Delay time (%s)', app.tUnitAx()));
                title(ax, app.mapTitle(scales{k}));
                
                hold(ax, 'on');
                app.hVLine(k) = xline(ax, app.selWL, 'k-', ...
                    'LineWidth', 1, 'Alpha', 0.7, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                app.hHLine(k) = yline(ax, app.selT, 'k-', ...
                    'LineWidth', 1, 'Alpha', 0.7, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                hold(ax, 'off');
            end
            
            % In split mode, put the colorbar on a dedicated host axes so the
            % two data axes share the same x-axis width.
            if ~isempty(app.cbarAx) && isvalid(app.cbarAx)
                app.cbarAx.CLim = [cMin, cMax];
                colormap(app.cbarAx, cmap);
                cb = colorbar(app.cbarAx);
                cb.Label.String = '\DeltaA or Intensity';
            end
            % Let the renderer flush before any further graphics ops
            % (drawSpectrum / drawKinetics / overlays) are queued.
            drawnow limitrate;
        end
        
        function h = plotHeatmapOn(app, ax, t, data)
            wl = app.wavelength(:).';
            t  = t(:).';
            [W, T] = meshgrid(wl, t);
            C = data.';
            Z = zeros(size(C));
            h = surface(ax, W, T, Z, C, ...
                'EdgeColor', 'none', 'FaceColor', 'flat');
            if any(isnan(C(:)))
                h.AlphaData = double(~isnan(C));
                h.AlphaDataMapping = 'none';
                h.FaceAlpha = 'flat';
            end
            view(ax, 2);
            axis(ax, 'tight');
            set(ax, 'YDir', 'normal');
        end
        
        function drawSpectrum(app)
            if isempty(app.deltaA), return; end
            ax = app.axSpec;  cla(ax);
            hold(ax, 'on');
            
            cmap = lines(max(8, numel(app.specOverlays)));
            legH = gobjects(0);  legL = {};
            for k = 1:numel(app.specOverlays)
                [~, idx] = min(abs(app.delay - app.specOverlays(k)));
                hO = plot(ax, app.wavelength, app.deltaA(:, idx), '-', ...
                    'Color', cmap(k, :), 'LineWidth', 1.0);
                legH(end+1) = hO;  %#ok<AGROW>
                legL{end+1} = sprintf('t = %.3g %s', app.delay(idx), app.tUnitAx());
            end
            
            [~, idx] = min(abs(app.delay - app.selT));
            hCur = plot(ax, app.wavelength, app.deltaA(:, idx), 'b-', ...
                'LineWidth', 1.6);
            legH(end+1) = hCur;
            legL{end+1} = sprintf('t = %.3g %s (current)', app.delay(idx), app.tUnitAx());
            
            yline(ax, 0, 'k:');
            % Vertical line marking the currently monitored wavelength
            xline(ax, app.selWL, 'k-', 'LineWidth', 1.0, 'Alpha', 0.7, ...
                'HitTest', 'off', 'PickableParts', 'none');
            hold(ax, 'off');
            xlabel(ax, app.xAxLabel());
            ylabel(ax, '\DeltaA or Intensity');
            grid(ax, 'on');
            xlim(ax, [min(app.wavelength), max(app.wavelength)]);
            if numel(legH) > 1
                legend(ax, legH, legL, 'Location', 'best', 'FontSize', 8);
            else
                legend(ax, 'off');
                title(ax, legL{end});
            end
        end
        
        function drawKinetics(app)
            if isempty(app.deltaA), return; end
            
            cmap = lines(max(8, numel(app.kinOverlays)));
            
            switch app.delayScaleMode
                case 'linear'
                    masks = {true(size(app.delay))};
                    scales = {'linear'};
                case 'log'
                    pos = app.delay > 0;
                    if ~any(pos), pos = true(size(app.delay)); end
                    masks = {pos};
                    scales = {'log'};
                case 'split'
                    thr = app.splitThreshold;
                    maskL = app.delay <= thr;
                    maskR = app.delay >= thr;
                    masks  = {maskL, maskR};
                    scales = {'linear', 'log'};
            end
            
            for k = 1:numel(app.axKin)
                ax = app.axKin(k);  cla(ax);  hold(ax, 'on');
                m  = masks{k};
                tSub = app.delay(m);
                if isempty(tSub), continue; end
                
                legH = gobjects(0);  legL = {};
                for j = 1:numel(app.kinOverlays)
                    [~, idx] = min(abs(app.wavelength - app.kinOverlays(j)));
                    hO = plot(ax, tSub, app.deltaA(idx, m), '-', ...
                        'Color', cmap(j, :), 'LineWidth', 1.0);
                    if k == 1
                        legH(end+1) = hO;  %#ok<AGROW>
                        legL{end+1} = sprintf('%s = %.1f %s', app.xSym(), ...
                            app.wavelength(idx), app.xUnitAx());
                    end
                end
                
                [~, idx] = min(abs(app.wavelength - app.selWL));
                hCur = plot(ax, tSub, app.deltaA(idx, m), 'r-', ...
                    'LineWidth', 1.6);
                if k == 1
                    legH(end+1) = hCur;
                    legL{end+1} = sprintf('%s = %.1f %s (current)', ...
                        app.xSym(), app.wavelength(idx), app.xUnitAx());
                end
                
                yline(ax, 0, 'k:');
                % Vertical line marking the currently monitored delay.
                % In split-view the log half ignores t <= 0, so only
                % draw the line where it falls inside the visible mask.
                if app.selT >= min(tSub) && app.selT <= max(tSub)
                    xline(ax, app.selT, 'k-', 'LineWidth', 1.0, ...
                        'Alpha', 0.7, 'HitTest', 'off', ...
                        'PickableParts', 'none');
                end
                hold(ax, 'off');
                xlabel(ax, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(ax, '\DeltaA or Intensity');
                grid(ax, 'on');
                set(ax, 'XScale', scales{k});
                [xlo, xhi] = app.clipViewToMask(tSub, scales{k});
                xlim(ax, [xlo, xhi]);
                % On log scale the kinetics panel is narrow, so MATLAB's
                % auto tick algorithm may drop alternate decades. Force
                % every decade in [xlo, xhi] to be tick-labeled.
                if strcmp(scales{k}, 'log') && xlo > 0 && xhi > xlo
                    decLo = ceil(log10(xlo));
                    decHi = floor(log10(xhi));
                    if decHi >= decLo
                        ax.XTick = 10.^(decLo:decHi);
                    end
                end
                if k == 1 && numel(legH) > 1
                    legend(ax, legH, legL, 'Location','best','FontSize',8);
                else
                    legend(ax, 'off');
                end
                if k == 1 && numel(app.axKin) == 1
                    title(ax, sprintf('%s = %.1f %s', app.xSym(), ...
                        app.wavelength(idx), app.xUnitAx()));
                end
            end
            
            % In split-view (two stacked panels), force the same y-range
            % on both panels so that the linear and log delay halves can
            % be compared consistently.  Compute the range from the data
            % directly (not from queried ylim, which can be stale or
            % "frozen" from a previous manual ylim() call when the
            % wavelength selection changes).
            if numel(app.axKin) == 2
                yLo =  Inf;  yHi = -Inf;
                idxList = unique([app.kinOverlays(:); app.selWL]);
                anyMask = false(size(app.delay));
                for k = 1:numel(masks)
                    anyMask = anyMask | masks{k};
                end
                for jj = 1:numel(idxList)
                    [~, idx] = min(abs(app.wavelength - idxList(jj)));
                    yj = app.deltaA(idx, anyMask);
                    yj = yj(isfinite(yj));
                    if isempty(yj), continue; end
                    yLo = min(yLo, min(yj));
                    yHi = max(yHi, max(yj));
                end
                if isfinite(yLo) && isfinite(yHi) && yLo < yHi
                    pad = 0.05 * (yHi - yLo);
                    for k = 1:numel(app.axKin)
                        ylim(app.axKin(k), [yLo - pad, yHi + pad]);
                    end
                end
            end
        end
        
        function updateCrosshairs(app)
            for k = 1:numel(app.hVLine)
                if isvalid(app.hVLine(k)), app.hVLine(k).Value = app.selWL; end
            end
            for k = 1:numel(app.hHLine)
                if isvalid(app.hHLine(k)), app.hHLine(k).Value = app.selT; end
            end
        end
        
        function str = mapTitle(app, scaleStr)
            parts = {};
            if app.bgApplied,    parts{end+1} = sprintf('BG(N=%d)', app.bgN); end
            if app.chirpApplied, parts{end+1} = 'chirp(Sellmeier)'; end
            if isempty(parts), base = 'raw'; else, base = strjoin(parts, ', '); end
            if app.isCropped(), base = [base ', cropped']; end
            if nargin > 1 && ~isempty(scaleStr)
                str = sprintf('2D \\DeltaA or Intensity  (%s, %s scale)', base, scaleStr);
            else
                str = sprintf('2D \\DeltaA or Intensity  (%s)', base);
            end
        end
        
        % =================================================================
        % Interaction
        % =================================================================
        function on2DClick(app, ev)
            if isempty(app.deltaA), return; end
            % Remember which axes the drag was started on.
            app.drag2DAx = [];
            if isprop(ev, 'Source') && isa(ev.Source, ...
                    'matlab.ui.control.UIAxes')
                app.drag2DAx = ev.Source;
            end
            % Set drag flag BEFORE apply, so a quick mouse-up that
            % fires during apply's drawnow can correctly turn it off
            % without being overwritten.
            app.dragging2D = true;
            app.applyMain2DPoint(ev.IntersectionPoint);
        end
        
        function on2DMotion(app)
            if ~app.dragging2D, return; end
            % Safety: if the cursor has left the figure window, the
            % WindowButtonUp event may never reach us, so force-end the
            % drag here instead of tracking a stale CurrentPoint.
            if ~app.mouseInsideFigure(app.fig)
                app.on2DUp();
                return;
            end
            ax = app.drag2DAx;
            if isempty(ax) || ~isvalid(ax)
                for k = 1:numel(app.ax2D)
                    if isvalid(app.ax2D(k))
                        ax = app.ax2D(k);  break;
                    end
                end
            end
            if isempty(ax), return; end
            cp = ax.CurrentPoint;
            if isempty(cp) || size(cp, 1) < 1, return; end
            pt = cp(1, 1:2);
            app.applyMain2DPoint(pt);
        end
        
        function on2DUp(app)
            % Handlers stay attached permanently; we just flip the flag.
            app.dragging2D = false;
        end
        
        function applyMain2DPoint(app, pt)
            if numel(pt) < 2, return; end
            wl = pt(1);  t = pt(2);
            wl = min(max(wl, min(app.wavelength)), max(app.wavelength));
            t  = min(max(t,  min(app.delay)),       max(app.delay));
            app.setSelWL(wl, false);
            app.setSelT (t,  false);
            app.updateCrosshairs();
            app.drawSpectrum();
            app.drawKinetics();
            drawnow limitrate;
        end
        
        function inside = mouseInsideFigure(~, figHandle)
            % Return true when the OS mouse cursor is currently inside
            % the given uifigure's drawable area. We use the system-wide
            % pointer location instead of figure.CurrentPoint because
            % the latter can freeze at the last inside position when the
            % cursor leaves the figure, which would defeat the drag
            % safety check.
            inside = false;
            if isempty(figHandle) || ~isvalid(figHandle), return; end
            try
                mp = get(groot, 'PointerLocation');  % screen pixels
                po = figHandle.Position;             % screen pixels
                if numel(mp) < 2 || numel(po) < 4, return; end
                inside = mp(1) >= po(1) && mp(1) <= po(1) + po(3) && ...
                         mp(2) >= po(2) && mp(2) <= po(2) + po(4);
            catch
                inside = false;
            end
        end
        
        function setSelWL(app, val, doUpdate)
            if nargin < 3, doUpdate = true; end
            if isempty(app.wavelength), return; end
            val = min(max(val, min(app.wavelength)), max(app.wavelength));
            app.selWL = val;
            app.edSelWL.Value = val;
            if doUpdate
                app.updateCrosshairs();
                app.drawKinetics();
            end
        end
        
        function setSelT(app, val, doUpdate)
            if nargin < 3, doUpdate = true; end
            if isempty(app.delay), return; end
            val = min(max(val, min(app.delay)), max(app.delay));
            app.selT = val;
            app.edSelT.Value = val;
            if doUpdate
                app.updateCrosshairs();
                app.drawSpectrum();
            end
        end
        
        function onScaleModeChange(app, val)
            app.delayScaleMode = lower(val);
            app.build2DAxes();
            app.buildKinAxes();
            app.updateAll();
            app.savePrefs();
        end
        
        function onThresholdChange(app, val)
            if isnan(val) || val <= 0, return; end
            app.splitThreshold = val;
            if strcmp(app.delayScaleMode, 'split')
                app.updateAll();
            end
            app.savePrefs();
        end
        
        function onMainViewChange(app, which, val)
            switch which
                case 'min', app.mainViewTMin = val;
                case 'max', app.mainViewTMax = val;
            end
            app.drawMap2D();
            app.drawKinetics();
        end
        
        function resetMainView(app)
            if isempty(app.delay)
                app.mainViewTMin = [];  app.mainViewTMax = [];
                return;
            end
            app.mainViewTMin = app.delay(1);
            app.mainViewTMax = app.delay(end);
            app.edViewTMin.Value = app.mainViewTMin;
            app.edViewTMax.Value = app.mainViewTMax;
            app.drawMap2D();
            app.drawKinetics();
        end
        
        % ---------- Colormap & Z-range handling ----------
        function cmap = getColormap(app)
            N = 256;
            switch lower(app.mapColormap)
                case 'turbo',   cmap = turbo(N);
                case 'parula',  cmap = parula(N);
                case 'jet',     cmap = jet(N);
                case 'hot',     cmap = hot(N);
                case 'cool',    cmap = cool(N);
                case 'gray',    cmap = gray(N);
                case 'rdbu'
                    cmap = app.divergingColormap(N, ...
                        [0.02 0.19 0.38], [1 1 1], [0.40 0 0.12]);
                case 'bwr'
                    cmap = app.divergingColormap(N, ...
                        [0 0 1], [1 1 1], [1 0 0]);
                case 'seismic'
                    cmap = app.divergingColormap(N, ...
                        [0 0 0.30], [1 1 1], [0.55 0 0]);
                otherwise
                    cmap = turbo(N);
            end
        end
        
        function cmap = divergingColormap(~, N, c1, cMid, c2)
            half = floor(N/2);
            t1 = linspace(0, 1, half).';
            t2 = linspace(0, 1, N - half).';
            part1 = (1 - t1).*c1 + t1.*cMid;
            part2 = (1 - t2).*cMid + t2.*c2;
            cmap = [part1; part2];
        end
        
        function [cMin, cMax] = computeZLim(app)
            if ~isempty(app.mapZMin) && ~isempty(app.mapZMax) && ...
                    app.mapZMin < app.mapZMax
                cMin = app.mapZMin;  cMax = app.mapZMax;
            else
                if isempty(app.deltaA)
                    cMin = -1;  cMax = 1;
                else
                    % Robust: use 99.5th percentile of |deltaA| so a
                    % handful of outliers (pump scattering, chirp-edge
                    % spikes) don't flatten the color contrast for the
                    % rest of the map.
                    vals = abs(app.deltaA(:));
                    vals = vals(isfinite(vals));
                    if isempty(vals)
                        cl = 1;
                    else
                        % (from the sorted values: prctile would need
                        % the Statistics Toolbox)
                        vs = sort(vals);
                        cl = vs(max(1, ceil(0.995 * numel(vs))));
                        if ~isfinite(cl) || cl == 0
                            cl = max(vals, [], 'omitnan');
                        end
                    end
                    if isempty(cl) || cl == 0 || isnan(cl), cl = 1; end
                    cMin = -cl;  cMax = cl;
                end
            end
        end
        
        function onColormapChange(app, val)
            app.mapColormap = val;
            cmap = app.getColormap();
            for k = 1:numel(app.ax2D)
                if isvalid(app.ax2D(k)), colormap(app.ax2D(k), cmap); end
            end
            if ~isempty(app.cbarAx) && isvalid(app.cbarAx)
                colormap(app.cbarAx, cmap);
            end
            % Sweep every axes inside open child windows so all 2D maps
            % follow the main-window colormap consistently. Axes that
            % don't show a colormap image are unaffected visibly.
            childFigs = app.childFigList();
            for kk = 1:numel(childFigs)
                cf = childFigs{kk};
                if isempty(cf) || ~isvalid(cf), continue; end
                axList = findall(cf, 'Type', 'axes');
                for ai = 1:numel(axList)
                    try, colormap(axList(ai), cmap); catch, end
                end
            end
            app.savePrefs();
        end
        
        function onZChange(app, which, val)
            if strcmp(which, 'min')
                newMin = val;                 newMax = app.edZMax.Value;
            else
                newMin = app.edZMin.Value;    newMax = val;
            end
            if newMin >= newMax
                % Invalid entry: silently revert field to current applied limit
                [cMin, cMax] = app.computeZLim();
                app.edZMin.Value = cMin;  app.edZMax.Value = cMax;
                return;
            end
            app.mapZMin = newMin;  app.mapZMax = newMax;
            app.applyZLim();
        end
        
        function applyZLim(app)
            [cMin, cMax] = app.computeZLim();
            for k = 1:numel(app.ax2D)
                if isvalid(app.ax2D(k)), app.ax2D(k).CLim = [cMin, cMax]; end
            end
            if ~isempty(app.cbarAx) && isvalid(app.cbarAx)
                app.cbarAx.CLim = [cMin, cMax];
            end
        end
        
        function resetZRange(app)
            app.mapZMin = [];  app.mapZMax = [];
            [cMin, cMax] = app.computeZLim();
            app.edZMin.Value = cMin;
            app.edZMax.Value = cMax;
            app.applyZLim();
        end
        
        % =================================================================
        % Session save / load
        % =================================================================
        function names = sessionFields(~)
            % Single source of truth for which properties round-trip
            % through a session file. Keep save and load symmetric.
            names = { ...
                ... % data
                'wavelength', 'delay', 'deltaA', 'deltaA_raw', ...
                'originalWavelength', 'originalDelay', ...
                'originalDeltaA', 'dataSourceDesc', ...
                ... % corrections
                'bgApplied', 'bgSpectrum', 'bgN', ...
                'chirpApplied', 'chirpPts', 'chirpParams', ...
                'chirpFitRMS', 'chirpViewTMin', 'chirpViewTMax', ...
                'ridgeMethod', 'ridgeTmin', 'ridgeTmax', ...
                'ridgeSmoothN', 'ridgeWL', 'ridgeT', ...
                'tZeroShift', 'maskedRegions', ...
                ... % selection / overlays / view
                'selWL', 'selT', 'specOverlays', 'kinOverlays', ...
                'delayScaleMode', 'splitThreshold', ...
                'mainViewTMin', 'mainViewTMax', ...
                'mapColormap', 'mapZMin', 'mapZMax', 'timeUnit', 'specUnit', ...
                ... % single-wavelength fit settings
                'kfNComp', 'kfTauInit', 'kfTauFixed', ...
                'kfBetaInit', 'kfBetaFixed', 'kfStretchOn', ...
                'kfHasInf', 'kfT0', 'kfT0Fixed', ...
                'kfFWHM', 'kfFWHMFixed', 'kfIRFTreat', 'kfFitTMin', 'kfFitTMax', ...
                ... % global analysis settings + results
                'gaNComp', 'gaTauInit', 'gaTauFixed', ...
                'gaBetaInit', 'gaStretchOn', 'gaIRFTreat', 'gaSeqModel', ...
                'gaHasInf', 'gaT0', 'gaT0Fixed', ...
                'gaFWHM', 'gaFWHMFixed', 'gaCA', 'gaResultCA', ...
                'gaResultTau', 'gaResultBeta', 'gaResultT0', ...
                'gaResultFWHM', 'gaResultHasInf', 'gaResultDADS', ...
                'gaResultEADS', 'gaResultFit', 'gaResultRMS', ...
                'gaResultTauStd', 'gaFitTMin', 'gaFitTMax', ...
                'gaResultStretchOn', 'gaResultIRF', 'gaResultTRange', ...
                ... % persisted child-window analysis results
                'cohSaved', 'ldaSaved', 'mcrSaved', 'svdSaved', ...
                'kfSaved', 'tfSaved', 'tgSet', 'tgSaved'};
        end
        
        function saveSession(app)
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');  return;
            end
            defName = 'TRSpec_session.mat';
            base = app.dataStem();
            if ~isempty(base)
                defName = ['TRSpec_' base '_session.mat'];
            end
            [fn, fp] = app.putFile('*.mat', 'Save session', defName);
            if isequal(fn, 0), return; end
            
            s = struct();
            flds = app.sessionFields();
            for k = 1:numel(flds)
                s.(flds{k}) = app.(flds{k});
            end
            s.meta = struct( ...
                'app',     'TRSpecAnalyzer', ...
                'version', 1, ...
                'saved',   datestr(now, 'yyyy-mm-dd HH:MM:SS'));  %#ok<TNOW1,DATST>
            
            try
                save(fullfile(fp, fn), '-struct', 's', '-v7.3');
            catch ME
                uialert(app.fig, sprintf('Save failed:\n%s', ...
                    ME.message), 'Session');
                return;
            end
            app.lblStatus.Text = sprintf('Session saved: %s', fn);
        end
        
        function L = childFigList(app)
            % Single source of truth for every child window handle.
            % Add new windows HERE so closing, colormap sweeps and
            % unit switches all see them.
            L = {app.bgFig, app.chirpFig, app.cropFig, app.loadFig, ...
                 app.gaFig, app.maskFig, app.kinFitFig, app.svdFig, ...
                 app.ldaFig, app.mcrFig, app.coherenceFig, app.tfFig, ...
                 app.targetFig, app.autoFig};
        end
        
        function closeAnalysisWindows(app)
            % Close the ANALYSIS windows only (not the correction tools
            % such as Crop / Chirp / Background / Mask / Load & Average,
            % one of which is usually the window that triggered the data
            % change and is still executing its callback).
            % The Auto Analysis window is closed as well: its report
            % describes the old data. (While it applies its own
            % corrections it takes itself off this list.)
            L = {app.tfFig, app.coherenceFig, app.gaFig, app.kinFitFig, ...
                 app.svdFig, app.ldaFig, app.mcrFig, app.targetFig, ...
                 app.autoFig};
            for kk = 1:numel(L)
                cf = L{kk};
                if ~isempty(cf) && isvalid(cf)
                    try, close(cf); catch, end
                end
            end
        end
        
        function dataChanged(app)
            % The working data (values or axes) changed. Every analysis
            % window holds local state computed from the old data (cached
            % results, edit-field limits, maps of the old size) and every
            % stored result describes the old data, so close the former
            % and drop the latter. Settings stored outside the results
            % (GA / single-wavelength fit setup) are kept.
            app.closeAnalysisWindows();
            app.invalidateChildResults();
        end
        
        function m = irfModeNorm(~, m)
            % IRF-model name as used by the fit windows: 'none' or
            % 'numerical'. Anything else -- e.g. 'skip' in a session
            % saved by an older version, an option that no longer
            % exists -- becomes 'numerical'.
            isNone = (ischar(m) || (isstring(m) && isscalar(m))) && ...
                strcmpi(m, 'none');
            if isNone
                m = 'none';
            else
                m = 'numerical';
            end
        end
        
        function [lo, hi] = clampFitRange(app, lo, hi)
            % Stored fit-time range -> a valid range on the CURRENT delay
            % axis. Empty, non-finite or non-overlapping -> full range.
            dLo = min(app.delay);  dHi = max(app.delay);
            if isempty(lo) || isempty(hi) || ~isfinite(lo) || ~isfinite(hi)
                lo = dLo;  hi = dHi;  return;
            end
            if lo > hi, tmp = lo; lo = hi; hi = tmp; end
            lo = max(lo, dLo);  hi = min(hi, dHi);
            if ~(hi > lo), lo = dLo;  hi = dHi; end
        end
        
        function [lo, hi] = fitRangeToStore(app, lo, hi)
            % A range that covers the whole delay axis is stored as
            % "no restriction" ([]), so it stays "full" when a longer
            % dataset is loaded later instead of silently becoming a
            % restriction to the old dataset's extent.
            dLo = min(app.delay);  dHi = max(app.delay);
            span = dHi - dLo;
            a = min(lo, hi);  b = max(lo, hi);
            if ~isfinite(a) || ~isfinite(b) || span <= 0 || ...
                    ((a - dLo) <= 1e-9 * span && (dHi - b) <= 1e-9 * span)
                lo = [];  hi = [];
            else
                lo = a;  hi = b;
            end
        end
        
        function closeChildWindows(app)
            % Close every open analysis / correction child window.
            % Used when the underlying dataset is replaced wholesale
            % (e.g. loading a session) so no window keeps operating on
            % stale data or axes.
            childFigs = app.childFigList();
            for kk = 1:numel(childFigs)
                cf = childFigs{kk};
                if ~isempty(cf) && isvalid(cf)
                    try, close(cf); catch, end
                end
            end
        end
        
        function claAll(~, ax)
            % cla() only removes children whose HandleVisibility is 'on';
            % reference lines created with 'HandleVisibility','off'
            % (xline / yline kept out of legends) survive it and pile up
            % on every redraw -- visibly when their position changes,
            % silently (as a leak) when it does not. Delete every child.
            % Then run the regular cla so everything it also does
            % (resetting the colour-order index, etc.) still happens.
            if isempty(ax) || ~isvalid(ax), return; end
            try
                delete(allchild(ax));
            catch
            end
            cla(ax);
        end
        
        % =================================================================
        % Preferences (persist across MATLAB restarts via prefdir)
        % =================================================================
        function pth = prefsFile(~)
            pth = fullfile(prefdir, 'TRSpecAnalyzer_prefs.mat');
        end
        
        function loadPrefs(app)
            % Read stored defaults into the app PROPERTIES only. Called
            % at the very top of createMainGUI, before any control or
            % axes exists, so the toolbar controls are created with
            % these values and the axes are built exactly once in the
            % stored scale mode. (Applying prefs through the UI change
            % handlers after construction would delete and re-create
            % the axes while the uifigure renderer is still starting
            % up, which leaves them blank.) Every step is guarded: a
            % missing or corrupt prefs file leaves the built-in defaults.
            try
                if ~exist(app.prefsFile(), 'file'), return; end
                p = load(app.prefsFile());
            catch
                return;
            end
            try
                if isfield(p, 'lastDir') && ischar(p.lastDir) && ...
                        ~isempty(p.lastDir) && isfolder(p.lastDir)
                    app.prefLastDir = p.lastDir;
                end
            catch, end
            try
                if isfield(p, 'timeUnit') && ischar(p.timeUnit) && ...
                        any(strcmp(p.timeUnit, {'fs', 'ps', 'ns', 'us'}))
                    app.timeUnit = p.timeUnit;
                    app.chirpViewTMax = app.defaultChirpTMax();
                end
            catch, end
            try
                if isfield(p, 'mapColormap') && ischar(p.mapColormap) && ...
                        any(strcmp(p.mapColormap, {'turbo', 'parula', 'jet', ...
                        'hot', 'cool', 'gray', 'RdBu', 'BWR', 'seismic'}))
                    app.mapColormap = p.mapColormap;
                end
            catch, end
            try
                if isfield(p, 'splitThreshold') && isnumeric(p.splitThreshold) ...
                        && isscalar(p.splitThreshold) && ...
                        isfinite(p.splitThreshold) && p.splitThreshold > 0
                    app.splitThreshold = p.splitThreshold;
                end
            catch, end
            try
                if isfield(p, 'delayScaleMode') && ischar(p.delayScaleMode) && ...
                        any(strcmp(p.delayScaleMode, {'linear', 'log', 'split'}))
                    app.delayScaleMode = p.delayScaleMode;
                end
            catch, end
            try
                if isfield(p, 'specUnit') && ischar(p.specUnit) && ...
                        any(strcmp(p.specUnit, {'nm', 'cm-1'}))
                    app.specUnit = p.specUnit;
                end
            catch, end
        end
        
        function savePrefs(app)
            % Best-effort; never interrupts the user if prefdir is
            % read-only.
            try
                lastDir        = app.prefLastDir;        %#ok<NASGU>
                timeUnit       = app.timeUnit;           %#ok<NASGU>
                mapColormap    = app.mapColormap;        %#ok<NASGU>
                delayScaleMode = app.delayScaleMode;     %#ok<NASGU>
                splitThreshold = app.splitThreshold;     %#ok<NASGU>
                specUnit       = app.specUnit;           %#ok<NASGU>
                save(app.prefsFile(), 'lastDir', 'timeUnit', ...
                    'mapColormap', 'delayScaleMode', 'splitThreshold', ...
                    'specUnit');
            catch
            end
        end
        
        function lbl = scaleModeLabel(app)
            % 'linear' | 'log' | 'split'  ->  dropdown item text
            switch lower(app.delayScaleMode)
                case 'linear', lbl = 'Linear';
                case 'log',    lbl = 'Log';
                otherwise,     lbl = 'Split';
            end
        end
        
        function d = startDir(app)
            % Folder a file dialog should open in: the last one used,
            % falling back to the current folder.
            d = app.prefLastDir;
            if isempty(d) || ~ischar(d) || ~isfolder(d), d = pwd; end
        end
        
        function rememberDir(app, pth)
            if ischar(pth) && ~isempty(pth) && isfolder(pth)
                app.prefLastDir = pth;
                app.savePrefs();
            end
        end
        
        function nm = dataStem(app)
            %dataStem  Short name of the loaded data, for the default file
            %   names of the save dialogs: the file name without folder
            %   and extension; for an average of several files the first
            %   file name followed by _avgN. Characters that a file name
            %   cannot hold, and blanks, become '_'. Empty when no name
            %   is known.
            nm = '';
            d = app.dataSourceDesc;
            if ~ischar(d) || isempty(d), return; end
            sfx = '';
            tok = regexp(d, '^Avg of (\d+) files \((.*)\.\.\.\)$', ...
                'tokens', 'once');
            if ~isempty(tok)
                d = tok{2};  sfx = ['_avg' tok{1}];
            elseif strncmp(d, 'Loaded: ', 8)
                d = d(9:end);
            end
            [~, nm] = fileparts(d);
            nm = regexprep([nm sfx], '[\\/:*?"<>|\s]+', '_');
        end

        function [file, path] = putFile(app, filterSpec, dlgTitle, defName)
            % uiputfile that starts in the last-used folder and
            % remembers the folder the user picks.
            [file, path] = uiputfile(filterSpec, dlgTitle, ...
                fullfile(app.startDir(), defName));
            if ~isequal(file, 0), app.rememberDir(path); end
        end
        
        function [nSug, sigma, txt] = suggestNComponents(app)
            % Rank estimate from the singular values of the current
            % data: count components that stand clearly above the
            % noise floor (median of the lower half of the spectrum).
            % Returns nSug (1..8), the singular values, and a short
            % text summary of the leading ratios sigma_k / sigma_1.
            nSug = 1;  sigma = [];  txt = '';
            if isempty(app.deltaA), return; end
            D = app.deltaA;  D(~isfinite(D)) = 0;
            sigma = svd(D, 'econ');
            n = numel(sigma);
            if n < 4, nSug = max(1, min(n, 2)); return; end
            noise = median(sigma(ceil(n/2):end));
            if ~(noise > 0), noise = eps; end
            nSug = nnz(sigma > 5 * noise);
            nSug = max(1, min(nSug, 8));
            kShow = min(8, n);
            parts = cell(1, kShow);
            for k = 1:kShow
                parts{k} = sprintf('%d:%.3g', k, sigma(k) / sigma(1));
            end
            txt = sprintf('sigma_k/sigma_1 = %s  (noise ~ %.2g)', ...
                strjoin(parts, '  '), noise / sigma(1));
        end
        
        function clearGaResults(app)
            app.gaResultTau    = [];
            app.gaResultBeta   = [];
            app.gaResultT0     = [];
            app.gaResultFWHM   = [];
            app.gaResultHasInf = false;
            app.gaResultDADS   = [];
            app.gaResultEADS   = [];
            app.gaResultFit    = [];
            app.gaResultRMS    = NaN;
            app.gaResultTauStd = [];
            app.gaResultStretchOn = [];
            app.gaResultIRF    = '';
            app.gaResultTRange = [];
            app.gaResultCA     = [];
        end
        
        function invalidateChildResults(app)
            % Drop all persisted analysis results. Called whenever the
            % underlying dataset changes shape or content in a way that
            % makes previous results meaningless (new file load, crop,
            % reset-all). Includes the GA results: gaResultDADS /
            % gaResultFit are sized to the wavelength / delay axes and
            % would otherwise crash plotDADS on the next window open.
            app.clearGaResults();
            app.cohSaved = [];
            app.ldaSaved = [];
            app.mcrSaved = [];
            app.kfSaved  = [];
            app.tfSaved  = [];
            app.tgSaved  = [];
            % tgSet (scheme and starting values) is a setting: kept.
            % svdSaved holds only view settings (N, nShow, scale),
            % which remain sensible across datasets -- keep it.
        end
        
        function loadSessionFile(app)
            [fn, fp] = uigetfile('*.mat', 'Load session', app.startDir());
            if isequal(fn, 0), return; end
            app.rememberDir(fp);
            try
                s = load(fullfile(fp, fn));
            catch ME
                uialert(app.fig, sprintf('Load failed:\n%s', ...
                    ME.message), 'Session');
                return;
            end
            if ~isfield(s, 'meta') || ~isfield(s.meta, 'app') || ...
                    ~strcmp(s.meta.app, 'TRSpecAnalyzer')
                uialert(app.fig, ['This .mat file does not look like ' ...
                    'a TRSpecAnalyzer session (missing meta tag).'], ...
                    'Session');
                return;
            end
            
            % Restore every recognized field; silently skip fields a
            % future/older session may not have so versions stay
            % mutually loadable. Close child windows first -- they
            % would otherwise keep referencing the previous dataset --
            % and drop any in-memory child results so a session file
            % without those fields does not inherit stale ones.
            app.closeChildWindows();
            app.invalidateChildResults();
            flds = app.sessionFields();
            for k = 1:numel(flds)
                if isfield(s, flds{k})
                    app.(flds{k}) = s.(flds{k});
                end
            end
            % Sessions of older versions: the "Skip pre-IRF region" IRF
            % option no longer exists (-> the Gaussian-IRF model), and
            % a session without a "Stretched step" entry was computed
            % with the default reading of stretched steps.
            app.kfIRFTreat = app.irfModeNorm(app.kfIRFTreat);
            app.gaIRFTreat = app.irfModeNorm(app.gaIRFTreat);
            if ~isempty(app.gaResultIRF)
                app.gaResultIRF = app.irfModeNorm(app.gaResultIRF);
            end
            if ~isfield(s, 'gaCA'),       app.gaCA = [];        end
            if ~isfield(s, 'gaResultCA'), app.gaResultCA = [];  end
            if ~isfield(s, 'gaSeqModel') || ...
                    ~any(strcmp(app.gaSeqModel, {'age', 't0'}))
                app.gaSeqModel = 'age';
            end
            
            % --- Sync toolbar controls to restored state ---
            scMap = struct('linear', 'Linear', 'log', 'Log', ...
                           'split', 'Split');
            if isfield(scMap, app.delayScaleMode)
                app.ddScale.Value = scMap.(app.delayScaleMode);
            end
            app.edThreshold.Value = app.splitThreshold;
            if ~isempty(app.mainViewTMin), app.edViewTMin.Value = app.mainViewTMin;
            else,                          app.edViewTMin.Value = 0; end
            if ~isempty(app.mainViewTMax), app.edViewTMax.Value = app.mainViewTMax;
            else,                          app.edViewTMax.Value = 0; end
            try, app.ddColormap.Value = app.mapColormap; catch, end
            try, app.ddTimeUnit.Value = app.timeUnit; catch, end
            app.refreshTimeUnitLabels();
            try, app.ddSpecUnit.Value = app.specUnit; catch, end
            app.refreshSpecUnitLabels();
            % Selection fields: set the LIMITS for the restored axes
            % first -- they still describe the previously loaded dataset,
            % and assigning a Value outside them is an error.
            if ~isempty(app.wavelength) && ~isempty(app.delay)
                wlLo = min(app.wavelength);  wlHi = max(app.wavelength);
                tLo  = min(app.delay);       tHi  = max(app.delay);
                if isempty(app.selWL) || ~isfinite(app.selWL)
                    app.selWL = app.wavelength(round(end/2));
                end
                if isempty(app.selT) || ~isfinite(app.selT)
                    app.selT = app.delay(1);
                end
                app.selWL = min(max(app.selWL, wlLo), wlHi);
                app.selT  = min(max(app.selT,  tLo),  tHi);
                app.edSelWL.Limits = [wlLo, wlHi];
                app.edSelT.Limits  = [tLo,  tHi ];
                app.edSelWL.Value  = app.selWL;
                app.edSelT.Value   = app.selT;
            end
            % Z-range fields show the applied limits (custom or auto).
            try
                [cMin, cMax] = app.computeZLim();
                app.edZMin.Value = cMin;  app.edZMax.Value = cMax;
            catch
            end
            
            % Rebuild axes layout for the restored scale mode, then
            % redraw everything. recompute() rebuilds the working data
            % from the restored raw data + corrections; hold the change
            % detection so a last-bit numerical difference (other MATLAB
            % version / machine) cannot wipe the restored results.
            app.holdDataChange = true;
            try
                app.build2DAxes();
                app.buildKinAxes();
                app.updateStatusLabel();
                app.updateAll();
            catch ME
                app.holdDataChange = false;
                rethrow(ME);
            end
            app.holdDataChange = false;
            drawnow;
            app.lblStatus.Text = sprintf('%s   [session: %s]', ...
                app.lblStatus.Text, fn);
        end
        
        function exportCorrected(app)
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');  return;
            end
            
            % Build suggestion filename reflecting current pipeline state
            parts = {};
            if app.bgApplied,    parts{end+1} = sprintf('BG%d', app.bgN);     end
            if app.chirpApplied, parts{end+1} = 'chirp';                     end
            if app.isCropped(),  parts{end+1} = 'cropped';                   end
            if isempty(parts),   suffix = 'raw'; else, suffix = strjoin(parts, '_'); end
            
            nameOnly = app.dataStem();
            if isempty(nameOnly), nameOnly = 'data'; end
            suggested = sprintf('%s_%s.csv', nameOnly, suffix);
            
            [file, path] = app.putFile( ...
                {'*.csv','CSV (comma-delimited, *.csv)'; ...
                 '*.dat','DAT (tab-delimited, *.dat)'; ...
                 '*.txt','TXT (tab-delimited, *.txt)'}, ...
                'Export corrected data', suggested);
            if isequal(file, 0), return; end
            fullpath = fullfile(path, file);
            
            [~, ~, ext] = fileparts(fullpath);
            if strcmpi(ext, '.csv'), delim = ','; else, delim = sprintf('\t'); end
            
            try
                app.writeDataFile(fullpath, ...
                    app.wavelength, app.delay, app.deltaA, delim);
                uialert(app.fig, ...
                    sprintf(['Saved %d wavelengths x %d delays.\n\n' ...
                             'Pipeline: %s\n\n' ...
                             'File: %s'], ...
                             numel(app.wavelength), numel(app.delay), ...
                             suffix, fullpath), ...
                    'Export complete', 'Icon', 'success');
            catch ME
                uialert(app.fig, sprintf('Export failed:\n%s', ME.message), ...
                    'Export error');
            end
        end
        
        function writeDataFile(~, fullpath, wl, t, A, delim)
            wl = wl(:);  t = t(:).';
            M = numel(wl);  N = numel(t);
            if ~isequal(size(A), [M N])
                error('Data size mismatch: A is %dx%d, expected %dx%d', ...
                    size(A, 1), size(A, 2), M, N);
            end
            out = zeros(M+1, N+1);
            out(1, 1)         = 0;          % corner cell (ignored on reading)
            out(1, 2:end)     = t;
            out(2:end, 1)     = wl;
            out(2:end, 2:end) = A;
            writematrix(out, fullpath, 'Delimiter', delim, ...
                'FileType', 'text');
        end
        
        % =================================================================
        % Overlays
        % =================================================================
        function pinSpec(app)
            if isempty(app.deltaA), return; end
            [~, idx] = min(abs(app.delay - app.selT));
            app.specOverlays(end+1, 1) = app.delay(idx);
            app.specOverlays = unique(app.specOverlays, 'stable');
            app.drawSpectrum();
        end
        function pinKin(app)
            if isempty(app.deltaA), return; end
            [~, idx] = min(abs(app.wavelength - app.selWL));
            app.kinOverlays(end+1, 1) = app.wavelength(idx);
            app.kinOverlays = unique(app.kinOverlays, 'stable');
            app.drawKinetics();
        end
        function clearSpecOverlays(app)
            app.specOverlays = zeros(0,1);  app.drawSpectrum();
        end
        function clearKinOverlays(app)
            app.kinOverlays = zeros(0,1);  app.drawKinetics();
        end
        
        % =================================================================
        % Reset
        % =================================================================
        function resetCorrections(app)
            % Restore the working data all the way back to the raw
            % snapshot taken at load time. Undoes crop, BG, chirp, mask,
            % and t=0 shift -- equivalent to re-loading the file.
            if isempty(app.originalDeltaA), return; end
            
            % Previous analysis results were computed on the corrected
            % data and no longer apply.
            app.dataChanged();
            
            % --- Working data back to the load-time snapshot ---
            app.wavelength = app.originalWavelength;
            app.delay      = app.originalDelay;
            app.deltaA_raw = app.originalDeltaA;
            
            % --- Clear all correction state ---
            app.bgApplied    = false;  app.bgSpectrum = [];  app.bgN = 0;
            app.chirpApplied = false;
            app.chirpPts     = zeros(0,2);
            app.chirpParams  = [];     app.chirpFitRMS = NaN;
            app.ridgeWL      = [];     app.ridgeT      = [];
            app.maskedRegions = cell(0, 3);
            app.tZeroShift   = 0;
            
            % --- Clamp the user selection and view to the new range ---
            wlLo = min(app.wavelength);  wlHi = max(app.wavelength);
            tLo  = min(app.delay);       tHi  = max(app.delay);
            app.selWL = min(max(app.selWL, wlLo), wlHi);
            app.selT  = min(max(app.selT,  tLo),  tHi);
            if ~isempty(app.edSelWL) && isvalid(app.edSelWL)
                app.edSelWL.Limits = [wlLo, wlHi];
                app.edSelWL.Value  = app.selWL;
            end
            if ~isempty(app.edSelT) && isvalid(app.edSelT)
                app.edSelT.Limits = [tLo, tHi];
                app.edSelT.Value  = app.selT;
            end
            app.mainViewTMin = tLo;  app.mainViewTMax = tHi;
            if ~isempty(app.edViewTMin) && isvalid(app.edViewTMin)
                app.edViewTMin.Value = app.mainViewTMin;
                app.edViewTMax.Value = app.mainViewTMax;
            end
            
            % --- Auto Z range based on raw data ---
            app.mapZMin = [];  app.mapZMax = [];
            if ~isempty(app.edZMin) && isvalid(app.edZMin)
                cl = max(abs(app.originalDeltaA(:)), [], 'omitnan');
                if isempty(cl) || cl == 0 || isnan(cl), cl = 1; end
                app.edZMin.Value = -cl;
                app.edZMax.Value =  cl;
            end
            
            % --- Clip pinned overlays to the new full range ---
            app.specOverlays = app.specOverlays( ...
                app.specOverlays >= tLo & app.specOverlays <= tHi);
            app.kinOverlays  = app.kinOverlays( ...
                app.kinOverlays  >= wlLo & app.kinOverlays  <= wlHi);
            
            app.updateAll();
        end
        
        % =================================================================
        % Zero-time shift (simple constant shift for ns-TA data)
        % =================================================================
        function applyZeroTimeShift(app, delta)
            % Shift the delay axis so that the time value `delta` in the
            % current axis becomes t = 0. delta is cumulative: applying
            % 2 then 0.5 results in a total 2.5 ps/us shift from the
            % original loaded data.
            if ~isfinite(delta) || delta == 0 || isempty(app.delay), return; end
            % The delay values are about to change under every stored
            % result (fit t0, coherence / time-frequency time axes,
            % restored traces): close the analysis windows and drop the
            % results first, while the state is still self-consistent.
            app.dataChanged();
            app.delay = app.delay - delta;
            app.tZeroShift = app.tZeroShift + delta;
            % Things that live on the delay axis move with it.
            app.specOverlays = app.specOverlays - delta;
            if ~isempty(app.gaFitTMin), app.gaFitTMin = app.gaFitTMin - delta; end
            if ~isempty(app.gaFitTMax), app.gaFitTMax = app.gaFitTMax - delta; end
            if ~isempty(app.kfFitTMin), app.kfFitTMin = app.kfFitTMin - delta; end
            if ~isempty(app.kfFitTMax), app.kfFitTMax = app.kfFitTMax - delta; end
            if isstruct(app.tgSet) && isfield(app.tgSet, 'fitTMin') && ...
                    ~isempty(app.tgSet.fitTMin)
                app.tgSet.fitTMin = app.tgSet.fitTMin - delta;
                app.tgSet.fitTMax = app.tgSet.fitTMax - delta;
            end
            % Keep related edit fields / limits consistent with the new axis.
            app.selT = app.selT - delta;
            if ~isempty(app.mainViewTMin)
                app.mainViewTMin = app.mainViewTMin - delta;
                app.mainViewTMax = app.mainViewTMax - delta;
            end
            if ~isempty(app.edSelT) && isvalid(app.edSelT)
                app.edSelT.Limits = [min(app.delay), max(app.delay)];
                newVal = app.selT;
                newVal = max(min(newVal, max(app.delay)), min(app.delay));
                app.edSelT.Value = newVal;
                app.selT = newVal;
            end
            if ~isempty(app.edViewTMin) && isvalid(app.edViewTMin)
                app.edViewTMin.Value = app.mainViewTMin;
                app.edViewTMax.Value = app.mainViewTMax;
            end
            app.updateAll();
        end
        
        function resetZeroTimeShift(app)
            if app.tZeroShift == 0, return; end
            app.applyZeroTimeShift(-app.tZeroShift);
            % applyZeroTimeShift leaves tZeroShift at exactly 0 already,
            % but guard against floating-point drift.
            app.tZeroShift = 0;
        end
        
        
        function applyZeroTimeAtCrosshair(app)
            % Shift the delay axis so that the currently selected
            % crosshair position on the main 2D map becomes t = 0.
            % The crosshair is moved by clicking on the 2D map (which
            % sets app.selT) or by editing the Delay field directly.
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');
                return;
            end
            delta = app.selT;
            if ~isfinite(delta)
                return;
            end
            if delta == 0
                uialert(app.fig, ...
                    ['Crosshair is already at t = 0. Click somewhere ' ...
                     'else on the 2D map first to pick a new t_0 ' ...
                     'location.'], 'No shift needed');
                return;
            end
            app.applyZeroTimeShift(delta);
            if ~isempty(app.lblStatus) && isvalid(app.lblStatus)
                app.lblStatus.Text = sprintf( ...
                    't_0 shift applied: %.4g %s  (total: %.4g %s)', ...
                    delta, app.tUnitTxt(), ...
                    app.tZeroShift, app.tUnitTxt());
            end
        end
        
        % =================================================================
        % Wavelength region masking (pump scattering etc.)
        % =================================================================
        function openMaskWindow(app)
            if isempty(app.deltaA_raw)
                uialert(app.fig, 'Load data first.', 'No data');  return;
            end
            if ~isempty(app.maskFig) && isvalid(app.maskFig)
                figure(app.maskFig);  return;
            end
            
            f = uifigure('Name', 'Mask Wavelength Regions', ...
                'Position', app.fitFigurePosition(960, 680));
            app.applyLightTheme(f);
            app.maskFig = f;
            
            gl = uigridlayout(f, [4, 1]);
            gl.RowHeight   = {'fit', '1x', 'fit', 180};
            gl.Padding     = [12 12 12 12];
            gl.RowSpacing  = 8;
            
            uilabel(gl, 'Text', ['Mask wavelength regions to exclude ' ...
                'from analysis. Click once on the 2D map to pick the ' ...
                'start wavelength (From), click again for the end ' ...
                '(To). Then choose a mode and press Add region. ' ...
                '"Set to NaN" removes the region from fits and shows ' ...
                'it as a transparent gap; "Set to 0" forces values to ' ...
                '0 while keeping them in the wavelength axis (not ' ...
                'recommended for fitting).'], ...
                'FontAngle', 'italic', 'FontColor', [0.3 0.3 0.3], ...
                'WordWrap', 'on');
            
            % 2D map preview
            mapPanel = uipanel(gl, 'BorderType', 'none');
            mgl = uigridlayout(mapPanel, [1, 1]);
            mgl.Padding = [0 0 0 0];
            ax = uiaxes(mgl);
            
            % Add-region controls row
            addRow = uigridlayout(gl, [1, 9]);
            addRow.ColumnWidth = {'fit', 90, 'fit', 90, 'fit', 150, ...
                                   120, '1x', 'fit'};
            addRow.ColumnSpacing = 6;  addRow.Padding = [0 0 0 0];
            uilabel(addRow, 'Text', sprintf('From (%s):', app.xUnitTxt()), ...
                'HorizontalAlignment', 'right');
            edWL1 = uieditfield(addRow, 'numeric', 'Value', 0);
            uilabel(addRow, 'Text', sprintf('To (%s):', app.xUnitTxt()), ...
                'HorizontalAlignment', 'right');
            edWL2 = uieditfield(addRow, 'numeric', 'Value', 0);
            uilabel(addRow, 'Text', 'Mode:', ...
                'HorizontalAlignment', 'right');
            ddMode = uidropdown(addRow, ...
                'Items',     {'Set to NaN (mask out)', 'Set to 0'}, ...
                'ItemsData', {'nan', 'zero'}, ...
                'Value', 'nan');
            btnAdd = uibutton(addRow, 'Text', 'Add region', ...
                'BackgroundColor', [0.30 0.55 0.80], 'FontColor', 'w');
            uilabel(addRow, 'Text', '');   % spacer
            lblPick = uilabel(addRow, 'Text', ...
                'Click map to pick From / To', ...
                'FontAngle', 'italic', 'FontColor', [0.45 0.45 0.45]);
            
            % Regions table + row buttons
            tblRow = uigridlayout(gl, [1, 2]);
            tblRow.ColumnWidth = {'1x', 150};
            tblRow.ColumnSpacing = 8;  tblRow.Padding = [0 0 0 0];
            tbl = uitable(tblRow, ...
                'ColumnName',    {sprintf('From (%s)', app.xUnitTxt()), ...
                                  sprintf('To (%s)', app.xUnitTxt()), 'Mode'}, ...
                'ColumnFormat',  {'numeric', 'numeric', 'char'}, ...
                'ColumnWidth',   {100, 100, 200}, ...
                'RowName', 'numbered', ...
                'ColumnEditable', [false false false]);
            btnCol = uigridlayout(tblRow, [3, 1]);
            btnCol.RowHeight = {36, 36, '1x'};
            btnCol.Padding   = [0 0 0 0];  btnCol.RowSpacing = 6;
            btnRemove = uibutton(btnCol, 'Text', 'Remove selected');
            btnClear  = uibutton(btnCol, 'Text', 'Clear all');
            uilabel(btnCol, 'Text', '');   % spacer
            
            % State within window
            pickStage = 1;   % next click sets From (1) or To (2)
            selRow    = [];  % currently highlighted table row
            
            drawMap();
            refreshTable();
            
            edWL1.ValueChangedFcn = @(~,~) updatePickStage();
            edWL2.ValueChangedFcn = @(~,~) updatePickStage();
            btnAdd.ButtonPushedFcn = @(~,~) doAdd();
            btnRemove.ButtonPushedFcn = @(~,~) doRemove();
            btnClear.ButtonPushedFcn = @(~,~) doClear();
            tbl.CellSelectionCallback = @(~,ev) onTableSelect(ev);
            f.CloseRequestFcn = @(~,~) closeAndClear();
            
            % ------------ nested helpers ------------
            function drawMap()
                app.claAll(ax);   % mask patches are hidden-handle objects
                wl = app.wavelength(:).';
                t  = app.delay(:).';
                [W, T] = meshgrid(wl, t);
                cl = max(abs(app.deltaA(:)), [], 'omitnan');
                if isempty(cl) || cl == 0 || isnan(cl), cl = 1; end
                sM = surf(ax, W, T, zeros(size(W)), app.deltaA.', ...
                    'EdgeColor', 'none');
                shading(ax, 'flat');
                view(ax, 2);
                colormap(ax, app.getColormap());
                try, clim(ax, [-cl cl]); catch, caxis(ax, [-cl cl]); end %#ok<CAXIS>
                colorbar(ax);
                xlabel(ax, app.xAxLabel());
                ylabel(ax, sprintf('Delay time (%s)', app.tUnitAx()));
                xlim(ax, [min(wl), max(wl)]);
                ylim(ax, [min(t), max(t)]);
                grid(ax, 'on');
                title(ax, 'Masked regions preview  (click to pick wavelengths)');
                ax.NextPlot = 'add';
                
                % Draw currently-applied mask regions as vertical bands
                yl = [min(t), max(t)];
                for k = 1:size(app.maskedRegions, 1)
                    w1 = app.maskedRegions{k, 1};
                    w2 = app.maskedRegions{k, 2};
                    col = [0.8 0.2 0.2];
                    if strcmp(app.maskedRegions{k, 3}, 'zero')
                        col = [0.2 0.5 0.8];
                    end
                    patch(ax, [w1 w2 w2 w1], [yl(1) yl(1) yl(2) yl(2)], ...
                        col, 'FaceAlpha', 0.18, 'EdgeColor', col, ...
                        'LineWidth', 1.2, 'HandleVisibility', 'off');
                end
                
                % Click handlers
                ax.ButtonDownFcn = @(~, ev) onMapClick(ev);
                kids = ax.Children;
                for kk = 1:numel(kids)
                    kids(kk).HitTest = 'on';
                    kids(kk).PickableParts = 'all';
                    try, kids(kk).ButtonDownFcn = @(~, ev) onMapClick(ev); catch, end
                end
            end
            
            function refreshTable()
                if isempty(app.maskedRegions)
                    tbl.Data = cell(0, 3);
                    selRow = [];
                    return;
                end
                d = cell(size(app.maskedRegions, 1), 3);
                for k = 1:size(app.maskedRegions, 1)
                    d{k, 1} = app.maskedRegions{k, 1};
                    d{k, 2} = app.maskedRegions{k, 2};
                    if strcmp(app.maskedRegions{k, 3}, 'nan')
                        d{k, 3} = 'Set to NaN';
                    else
                        d{k, 3} = 'Set to 0';
                    end
                end
                tbl.Data = d;
            end
            
            function onMapClick(ev)
                if isempty(ev) || ~isprop(ev, 'IntersectionPoint'), return; end
                p = ev.IntersectionPoint;
                if numel(p) < 1, return; end
                wlClick = p(1);
                if pickStage == 1
                    edWL1.Value = wlClick;
                    pickStage = 2;
                    lblPick.Text = 'Next click sets the To (end) wavelength';
                else
                    edWL2.Value = wlClick;
                    pickStage = 1;
                    lblPick.Text = ...
                        'Next click restarts picking from the From edge';
                end
            end
            
            function updatePickStage()
                % Manual edit -> reset pick stage
                pickStage = 1;
                lblPick.Text = 'Click map to pick From / To';
            end
            
            function doAdd()
                wl1 = edWL1.Value;  wl2 = edWL2.Value;
                if wl1 > wl2, tmp = wl1; wl1 = wl2; wl2 = tmp; end
                if ~isfinite(wl1) || ~isfinite(wl2) || wl1 >= wl2
                    uialert(f, ...
                        'Invalid wavelength range. From must be < To.', ...
                        'Input error');
                    return;
                end
                wlLo = min(app.wavelength);  wlHi = max(app.wavelength);
                if wl2 < wlLo || wl1 > wlHi
                    uialert(f, sprintf( ...
                        'Region [%.1f, %.1f] %s is outside data range [%.1f, %.1f] %s.', ...
                        wl1, wl2, app.xUnitTxt(), wlLo, wlHi, app.xUnitTxt()), 'Input error');
                    return;
                end
                app.maskedRegions(end+1, :) = {wl1, wl2, ddMode.Value};
                refreshTable();
                app.updateAll();
                drawMap();
                pickStage = 1;
                lblPick.Text = ...
                    sprintf('Added region %d.  Click map to pick the next From / To.', ...
                        size(app.maskedRegions, 1));
            end
            
            function onTableSelect(ev)
                if isempty(ev.Indices)
                    selRow = [];
                else
                    selRow = ev.Indices(1, 1);
                end
            end
            
            function doRemove()
                if isempty(selRow) || ...
                        selRow > size(app.maskedRegions, 1)
                    uialert(f, 'Select a row in the table first.', ...
                        'Nothing selected');
                    return;
                end
                app.maskedRegions(selRow, :) = [];
                selRow = [];
                refreshTable();
                app.updateAll();
                drawMap();
            end
            
            function doClear()
                if isempty(app.maskedRegions), return; end
                app.maskedRegions = cell(0, 3);
                refreshTable();
                app.updateAll();
                drawMap();
            end
            
            function closeAndClear()
                delete(f);
                app.maskFig = [];
            end
        end
        
        % =================================================================
        % Single-trace kinetic fit (one wavelength at a time)
        % =================================================================
        function openKineticFitWindow(app)
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');
                return;
            end
            if ~isempty(app.kinFitFig) && isvalid(app.kinFitFig)
                figure(app.kinFitFig);
                return;
            end
            
            f = uifigure('Name', 'Single-Wavelength Kinetic Fit', ...
                'Position', app.fitFigurePosition(1180, 720));
            app.applyLightTheme(f);
            app.kinFitFig = f;
            
            gl = uigridlayout(f, [1, 2]);
            gl.ColumnWidth = {340, '1x'};
            gl.Padding     = [10 10 10 10];
            gl.ColumnSpacing = 8;
            
            % ============ Left: setup ============
            left = uipanel(gl, 'Title', 'Fit setup');
            left.Layout.Column = 1;
            lgl = uigridlayout(left, [12, 1]);
            lgl.RowHeight = {30, 30, 30, 22, 150, 30, 30, 30, 30, 32, 22, '1x'};
            lgl.Padding = [8 8 8 8];  lgl.RowSpacing = 6;
            
            % Wavelength row
            rowWL = uigridlayout(lgl, [1, 4]);
            rowWL.ColumnWidth = {'fit', 90, 'fit', '1x'};
            rowWL.Padding = [0 0 0 0];  rowWL.ColumnSpacing = 6;
            uilabel(rowWL, 'Text', sprintf('%s (%s):', app.xSymTxt(), app.xUnitTxt()));
            edWL = uieditfield(rowWL, 'numeric', 'Value', app.selWL, ...
                'Limits', [min(app.wavelength), max(app.wavelength)]);
            btnUseSel = uibutton(rowWL, 'Text', 'Use crosshair');
            uilabel(rowWL, 'Text', '');
            
            % Components + tau=inf
            rowC = uigridlayout(lgl, [1, 3]);
            rowC.ColumnWidth = {'fit', 60, 'fit'};
            rowC.Padding = [0 0 0 0];  rowC.ColumnSpacing = 6;
            uilabel(rowC, 'Text', 'Components:');
            ddN = uidropdown(rowC, 'Items', {'1','2','3','4','5'}, ...
                'Value', num2str(app.kfNComp));
            cbInf = uicheckbox(rowC, 'Text', app.uiStr('Include \tau = \infty'), ...
                'Value', app.kfHasInf);
            
            % Trace scale
            rowSc = uigridlayout(lgl, [1, 4]);
            rowSc.ColumnWidth = {'fit', 90, 'fit', '1x'};
            rowSc.Padding = [0 0 0 0];  rowSc.ColumnSpacing = 6;
            uilabel(rowSc, 'Text', 'Trace scale:');
            ddScale = uidropdown(rowSc, 'Items', {'Linear', 'Log'}, ...
                'Value', 'Log');
            uilabel(rowSc, 'Text', '');  uilabel(rowSc, 'Text', '');
            
            uilabel(lgl, 'Text', sprintf( ...
                'Initial time constants (%s, positive):', app.tUnitTxt()), ...
                'FontAngle', 'italic', 'FontColor', [0.3 0.3 0.3]);
            
            tauTable = uitable(lgl, ...
                'ColumnName',   app.uiStr({ ...
                                 sprintf('\\tau_{init} (%s)', app.tUnitTxt()), ...
                                 '\tau fixed', '\beta_{init}', 'Stretched'}), ...
                'ColumnFormat', {'numeric', 'logical', 'numeric', 'logical'}, ...
                'ColumnEditable', [true true true true], ...
                'ColumnWidth',  {85, 70, 70, 70}, ...
                'RowName', 'numbered');
            
            % IRF t0
            rowT0 = uigridlayout(lgl, [1, 3]);
            rowT0.ColumnWidth = {'fit', 80, 'fit'};
            rowT0.Padding = [0 0 0 0];  rowT0.ColumnSpacing = 6;
            uilabel(rowT0, 'Text', sprintf('IRF t_0 (%s):', app.tUnitTxt()));
            edT0 = uieditfield(rowT0, 'numeric', 'Value', app.kfT0);
            cbT0fix = uicheckbox(rowT0, 'Text', 'Fixed', ...
                'Value', app.kfT0Fixed);
            
            % IRF FWHM
            rowFW = uigridlayout(lgl, [1, 3]);
            rowFW.ColumnWidth = {'fit', 80, 'fit'};
            rowFW.Padding = [0 0 0 0];  rowFW.ColumnSpacing = 6;
            uilabel(rowFW, 'Text', sprintf('IRF FWHM (%s):', app.tUnitTxt()));
            edFW = uieditfield(rowFW, 'numeric', 'Value', app.kfFWHM, ...
                'Limits', [1e-6 Inf]);
            cbFWfix = uicheckbox(rowFW, 'Text', 'Fixed', ...
                'Value', app.kfFWHMFixed);
            
            % IRF model: Gaussian IRF included ("Numerical conv.") or
            % none (pure exponentials, for a fit restricted to the
            % post-IRF region). A value left by an older version
            % ('skip') is mapped to the IRF model first: the dropdown
            % would reject it.
            app.kfIRFTreat = app.irfModeNorm(app.kfIRFTreat);
            rowIRF = uigridlayout(lgl, [1, 3]);
            rowIRF.ColumnWidth = {'fit', 165, '1x'};
            rowIRF.Padding = [0 0 0 0];  rowIRF.ColumnSpacing = 6;
            uilabel(rowIRF, 'Text', 'IRF model:');
            ddIRF = uidropdown(rowIRF, ...
                'Items',     {'None (pure exponentials)', ...
                              'Numerical conv.'}, ...
                'ItemsData', {'none', 'numerical'}, ...
                'Value',     app.kfIRFTreat, ...
                'Tooltip', ['Numerical conv.: the model is convolved ' ...
                    'with the Gaussian IRF (t0, FWHM). Standard exp ' ...
                    'components use the closed-form erfcx expression, ' ...
                    'stretched components a numerical quadrature at ' ...
                    'each delay.' newline ...
                    'None: no IRF convolution -- fit a plain ' ...
                    'sum of exponentials exp(-(t-t0)/tau)*H(t-t0). ' ...
                    'Combine with "Fit t range" to fit only the ' ...
                    'post-IRF part of the trace; t0 and FWHM are ' ...
                    'excluded from the fit.'], ...
                'ValueChangedFcn', @(~,~) onIRFModeChange());
            uilabel(rowIRF, 'Text', '');
            
            % Fit time range: points outside are excluded from the fit.
            rowFitT = uigridlayout(lgl, [1, 5]);
            rowFitT.ColumnWidth = {'fit', 62, 'fit', 62, 'fit'};
            rowFitT.Padding = [0 0 0 0];  rowFitT.ColumnSpacing = 4;
            uilabel(rowFitT, 'Text', sprintf('Fit t range (%s):', ...
                app.tUnitTxt()));
            edFitTmin = uieditfield(rowFitT, 'numeric', ...
                'Value', min(app.delay), ...
                'Tooltip', 'Lower bound of the fitted delay range.');
            uilabel(rowFitT, 'Text', '~', ...
                'HorizontalAlignment', 'center');
            edFitTmax = uieditfield(rowFitT, 'numeric', ...
                'Value', max(app.delay), ...
                'Tooltip', 'Upper bound of the fitted delay range.');
            uibutton(rowFitT, 'Text', 'Full', ...
                'Tooltip', 'Reset to the full delay range.', ...
                'ButtonPushedFcn', @(~,~) setFitTFull());
            
            % Run / Reset / Export
            rowBtn = uigridlayout(lgl, [1, 3]);
            rowBtn.ColumnWidth = {'1x', '1x', '1x'};
            rowBtn.Padding = [0 0 0 0];  rowBtn.ColumnSpacing = 5;
            btnRun = uibutton(rowBtn, 'Text', 'Run Fit', ...
                'BackgroundColor', [0.30 0.65 0.35], 'FontColor', 'w');
            btnReset = uibutton(rowBtn, 'Text', 'Reset');
            btnExport = uibutton(rowBtn, 'Text', 'Export...');
            
            lblStatus = uilabel(lgl, 'Text', 'Ready.', ...
                'FontColor', [0.2 0.2 0.2]);
            
            txtResult = uitextarea(lgl, 'Editable', 'off', ...
                'Value', {'Run a fit to see results.'}, ...
                'FontName', 'Monospaced');
            
            % ============ Right: data + fit + residual ============
            right = uipanel(gl, 'Title', 'Trace and fit');
            right.Layout.Column = 2;
            rgl = uigridlayout(right, [2, 1]);
            rgl.RowHeight = {'2x', '1x'};
            rgl.Padding = [6 6 6 6];  rgl.RowSpacing = 6;
            
            axTrace = uiaxes(rgl);
            axTrace.Layout.Row = 1;
            xlabel(axTrace, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axTrace, '\DeltaA or Intensity');
            grid(axTrace, 'on');
            title(axTrace, 'Click "Run Fit" to fit the kinetic trace');
            
            axRes = uiaxes(rgl);
            axRes.Layout.Row = 2;
            xlabel(axRes, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axRes, 'Residual (\DeltaA or Intensity)');
            grid(axRes, 'on');
            title(axRes, 'Residual');
            
            % Keep the two axes' x-limits in sync.
            try, linkaxes([axTrace, axRes], 'x'); catch, end
            
            % --- state for last fit (in closure) ---
            lastFit = struct('wlActual', NaN, 'iw', NaN, ...
                'data', [], 't', [], ...
                'tau', [], 't0', NaN, 'fwhm', NaN, 'hasInf', false, ...
                'A', [], 'fitVec', [], 'res', [], 'rms', NaN);
            
            fillTable();
            onIRFModeChange();
            % Restore a previously used fit time range.
            [edFitTmin.Value, edFitTmax.Value] = app.clampFitRange( ...
                app.kfFitTMin, app.kfFitTMax);
            drawTraceOnly();
            
            % ---- Restore previous fit (if any) ----
            % Fit settings persist in app.kf* by design; the result
            % snapshot in app.kfSaved lets a reopened window show the
            % last fit again. Only restore if the delay axis still
            % matches (crop / re-load invalidates the trace).
            if ~isempty(app.kfSaved) && isstruct(app.kfSaved) && ...
                    ~isempty(app.kfSaved.lastFit.fitVec) && ...
                    numel(app.kfSaved.lastFit.t) == numel(app.delay)
                try
                    lastFit = app.kfSaved.lastFit;
                    edWL.Value = lastFit.wlActual;
                    drawWithFit();
                    txtResult.Value = app.kfSaved.lines;
                    lblStatus.Text  = sprintf('%s  [restored]', ...
                        app.kfSaved.status);
                catch
                    lblStatus.Text = 'Stored result could not be restored -- run again.';
                end
            end
            
            ddN.ValueChangedFcn      = @(~,~) onNCompChange();
            ddScale.ValueChangedFcn  = @(~,~) refreshAxesScale();
            edWL.ValueChangedFcn     = @(~,~) drawTraceOnly();
            btnUseSel.ButtonPushedFcn = @(~,~) useCrosshair();
            btnRun.ButtonPushedFcn    = @(~,~) doRun();
            btnReset.ButtonPushedFcn  = @(~,~) doReset();
            btnExport.ButtonPushedFcn = @(~,~) doExport();
            f.CloseRequestFcn         = @(~,~) closeAndClear();
            
            % --- nested helpers ---
            function fillTable()
                N = app.kfNComp;
                if numel(app.kfTauInit)  ~= N
                    app.kfTauInit  = adjustVec(app.kfTauInit, N, 1);
                end
                if numel(app.kfTauFixed) ~= N
                    app.kfTauFixed = adjustVec(app.kfTauFixed, N, false);
                end
                if numel(app.kfBetaInit) ~= N
                    app.kfBetaInit  = adjustVec(app.kfBetaInit, N, 1);
                end
                if numel(app.kfBetaFixed) ~= N
                    app.kfBetaFixed = adjustVec(app.kfBetaFixed, N, true);
                end
                if numel(app.kfStretchOn) ~= N
                    app.kfStretchOn = adjustVec(app.kfStretchOn, N, false);
                end
                D = cell(N, 4);
                for ii = 1:N
                    D{ii, 1} = app.kfTauInit(ii);
                    D{ii, 2} = logical(app.kfTauFixed(ii));
                    D{ii, 3} = app.kfBetaInit(ii);
                    D{ii, 4} = logical(app.kfStretchOn(ii));
                end
                tauTable.Data = D;
            end
            
            function v = adjustVec(v0, N, pad)
                v = v0(:).';
                if numel(v) < N
                    v = [v, repmat(pad, 1, N - numel(v))];
                else
                    v = v(1:N);
                end
            end
            
            function readTable()
                N = size(tauTable.Data, 1);
                tau = zeros(1, N);  fix = false(1, N);
                bet = zeros(1, N);  son = false(1, N);
                for ii = 1:N
                    tau(ii) = tauTable.Data{ii, 1};
                    fix(ii) = logical(tauTable.Data{ii, 2});
                    bet(ii) = tauTable.Data{ii, 3};
                    son(ii) = logical(tauTable.Data{ii, 4});
                end
                app.kfTauInit   = tau;
                app.kfTauFixed  = fix;
                app.kfBetaInit  = bet;
                app.kfStretchOn = son;
                % beta-fixed semantics for the fitter: if the user typed
                % beta_init=1 AND has Stretched off, beta is irrelevant
                % (we use plain exp). If Stretched is on, beta is treated
                % as a free parameter unless the user enters exactly an
                % integer-1 placeholder in beta_init (heuristic). We
                % expose a separate "beta fixed" toggle in a future
                % version; for now keep it free whenever stretched is on.
                app.kfBetaFixed = ~son;   % only "fixed" when stretched off
            end
            
            function onNCompChange()
                app.kfNComp = str2double(ddN.Value);
                fillTable();
            end
            
            function onIRFModeChange()
                % 'none' removes the IRF from the model: FWHM is
                % meaningless and t0 is degenerate with amplitudes.
                isNone = strcmp(ddIRF.Value, 'none');
                if isNone
                    edFW.Enable    = 'off';
                    cbFWfix.Enable = 'off';
                    cbT0fix.Value  = true;
                    cbT0fix.Enable = 'off';
                else
                    edFW.Enable    = 'on';
                    cbFWfix.Enable = 'on';
                    cbT0fix.Enable = 'on';
                end
            end
            
            function setFitTFull()
                edFitTmin.Value = min(app.delay);
                edFitTmax.Value = max(app.delay);
            end
            
            function fr = getFitTRangeIfRestricted()
                % [tmin tmax] when the fit range is meaningfully
                % narrower than the full delay axis; [] otherwise.
                fr = [];
                lo = min(edFitTmin.Value, edFitTmax.Value);
                hi = max(edFitTmin.Value, edFitTmax.Value);
                if ~isfinite(lo) || ~isfinite(hi), return; end
                dLo = min(app.delay);  dHi = max(app.delay);
                span = dHi - dLo;
                if span <= 0, return; end
                if (lo - dLo) > 1e-9 * span || (dHi - hi) > 1e-9 * span
                    fr = [lo, hi];
                end
            end
            
            function useCrosshair()
                edWL.Value = app.selWL;
                drawTraceOnly();
            end
            
            function refreshAxesScale()
                isLog = strcmpi(ddScale.Value, 'Log') && any(app.delay > 0);
                if isLog
                    set(axTrace, 'XScale', 'log');
                    set(axRes,   'XScale', 'log');
                    pos = app.delay(app.delay > 0);
                    if ~isempty(pos)
                        xl = [min(pos), max(app.delay)];
                    else
                        xl = [1, 10];
                    end
                else
                    set(axTrace, 'XScale', 'linear');
                    set(axRes,   'XScale', 'linear');
                    xl = [min(app.delay), max(app.delay)];
                end
                % When the fit is restricted to a delay sub-range, zoom
                % the display onto that range.
                fr = getFitTRangeIfRestricted();
                if ~isempty(fr)
                    if isLog, fr(1) = max(fr(1), xl(1)); end
                    if fr(1) < fr(2), xl = fr; end
                end
                xlim(axTrace, xl);
                xlim(axRes,   xl);
            end
            
            function drawTraceOnly()
                % Show data only (no fit yet, or before re-fitting).
                wlReq = edWL.Value;
                [~, iw] = min(abs(app.wavelength - wlReq));
                wlActual = app.wavelength(iw);
                d = app.deltaA(iw, :);
                app.claAll(axTrace);
                plot(axTrace, app.delay, d, 'k.', 'MarkerSize', 8, ...
                    'DisplayName', 'Data');
                yline(axTrace, 0, 'k:', 'HandleVisibility', 'off');
                grid(axTrace, 'on');
                xlabel(axTrace, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axTrace, '\DeltaA or Intensity');
                title(axTrace, sprintf( ...
                    'Trace at %s = %.2f %s  (requested %.2f %s)', ...
                    app.xSym(), wlActual, app.xUnitAx(), wlReq, app.xUnitAx()));
                legend(axTrace, 'Location', 'best');
                cla(axRes);
                title(axRes, 'Residual (no fit yet)');
                grid(axRes, 'on');
                xlabel(axRes, sprintf('Delay time (%s)', app.tUnitAx()));
                refreshAxesScale();
            end
            
            function doReset()
                app.kfNComp     = 3;
                app.kfTauInit   = [1 10 100];
                app.kfTauFixed  = [false false false];
                app.kfBetaInit  = [1 1 1];
                app.kfBetaFixed = [true true true];
                app.kfStretchOn = [false false false];
                app.kfHasInf    = false;
                app.kfT0        = 0;       app.kfT0Fixed   = true;
                app.kfFWHM      = 0.15;    app.kfFWHMFixed = true;
                app.kfIRFTreat  = 'numerical';
                ddN.Value       = num2str(app.kfNComp);
                cbInf.Value     = app.kfHasInf;
                edT0.Value      = app.kfT0;     cbT0fix.Value = app.kfT0Fixed;
                edFW.Value      = app.kfFWHM;   cbFWfix.Value = app.kfFWHMFixed;
                ddIRF.Value     = app.kfIRFTreat;
                app.kfFitTMin   = [];
                app.kfFitTMax   = [];
                setFitTFull();
                fillTable();
                onIRFModeChange();
                lblStatus.Text  = 'Reset to defaults.';
                app.kfSaved = [];
            end
            
            function doRun()
                readTable();
                app.kfHasInf    = cbInf.Value;
                app.kfT0        = edT0.Value;
                app.kfT0Fixed   = cbT0fix.Value;
                app.kfFWHM      = edFW.Value;
                app.kfFWHMFixed = cbFWfix.Value;
                app.kfIRFTreat  = ddIRF.Value;
                if strcmp(app.kfIRFTreat, 'none')
                    app.kfT0Fixed   = true;
                    app.kfFWHMFixed = true;
                end
                [app.kfFitTMin, app.kfFitTMax] = app.fitRangeToStore( ...
                    edFitTmin.Value, edFitTmax.Value);
                tSelChk = app.delay >= min([edFitTmin.Value, edFitTmax.Value]) & ...
                          app.delay <= max([edFitTmin.Value, edFitTmax.Value]);
                if nnz(tSelChk) < (app.kfNComp + 2)
                    uialert(f, ['Fit t range contains too few delay ' ...
                        'points. Widen the range.'], 'Invalid range');
                    return;
                end
                
                if any(app.kfTauInit <= 0)
                    uialert(f, 'All tau values must be positive.', ...
                        'Invalid input');  return;
                end
                if any(app.kfBetaInit <= 0) || any(app.kfBetaInit > 2)
                    uialert(f, ['beta_init must be in (0, 2]. ' ...
                        'Use 1 for ordinary exponential.'], ...
                        'Invalid input');  return;
                end
                if app.kfFWHM <= 0 && ~strcmp(app.kfIRFTreat, 'none')
                    uialert(f, 'IRF FWHM must be positive.', ...
                        'Invalid input');  return;
                end
                
                wlReq = edWL.Value;
                [~, iw] = min(abs(app.wavelength - wlReq));
                wlActual = app.wavelength(iw);
                d = app.deltaA(iw, :);
                if all(~isfinite(d))
                    uialert(f, 'Trace at this wavelength is all NaN.', ...
                        'No data');  return;
                end
                
                lblStatus.Text = 'Fitting...';
                drawnow;
                try
                    [tauOut, betaOut, t0Out, fwhmOut, A, fitVec, ...
                        resVec, info] = app.fitSingleTrace( ...
                            app.delay, d, ...
                            app.kfTauInit, app.kfTauFixed, ...
                            app.kfBetaInit, app.kfBetaFixed, ...
                            app.kfStretchOn, ...
                            app.kfT0, app.kfT0Fixed, ...
                            app.kfFWHM, app.kfFWHMFixed, ...
                            app.kfHasInf, app.kfIRFTreat, ...
                            [app.kfFitTMin, app.kfFitTMax]);
                catch ME
                    uialert(f, sprintf('Fit failed:\n%s', ME.message), ...
                        'Fit Error');
                    lblStatus.Text = 'Fit failed.';
                    lastFit = struct('wlActual', NaN, 'iw', NaN, ...
                        'data', [], 't', [], 'mask', [], ...
                        'tau', [], 'beta', [], 'stretchOn', [], ...
                        't0', NaN, 'fwhm', NaN, ...
                        'hasInf', false, ...
                        'A', [], 'fitVec', [], 'res', [], 'rms', NaN);
                    return;
                end
                
                lastFit.wlActual  = wlActual;
                lastFit.iw        = iw;
                lastFit.data      = d(:);
                lastFit.t         = app.delay(:);
                lastFit.mask      = info.mask;
                lastFit.tau       = tauOut;
                lastFit.beta      = betaOut;
                lastFit.stretchOn = app.kfStretchOn;
                lastFit.t0        = t0Out;
                lastFit.fwhm      = fwhmOut;
                lastFit.hasInf    = app.kfHasInf;
                lastFit.A         = A(:).';
                lastFit.fitVec    = fitVec(:);
                lastFit.res       = resVec(:);
                lastFit.rms       = info.rms;
                
                drawWithFit();
                
                % Result text
                lines = {};
                lines{end+1} = sprintf('%-11s: %.2f %s', ...
                    ternary(strcmp(app.specUnit, 'cm-1'), 'Wavenumber', 'Wavelength'), ...
                    wlActual, app.xUnitTxt());
                lines{end+1} = sprintf('Iterations : %d', info.iter);
                lines{end+1} = sprintf('Final   RMS: %.4g', info.rms);
                if info.exitFlag <= 0
                    lines{end+1} = '';
                    lines{end+1} = '** NOT CONVERGED **';
                end
                lines{end+1} = '';
                lines{end+1} = 'Time constants & beta:';
                for ii = 1:numel(tauOut)
                    if app.kfStretchOn(ii)
                        tag = sprintf(' [stretched, beta=%.3f]', betaOut(ii));
                        if isfield(info, 'betaStd') && ...
                                isfinite(info.betaStd(ii))
                            tag = sprintf( ...
                                ' [stretched, beta=%.3f+/-%.3f]', ...
                                betaOut(ii), info.betaStd(ii));
                        end
                    else
                        tag = '';
                    end
                    if app.kfTauFixed(ii)
                        tag = [tag, ' (tau fixed)']; %#ok<AGROW>
                    end
                    tag = [tag, app.undetNote(info, 'tauStd', ii)]; %#ok<AGROW>
                    if isfield(info, 'tauStd') && isfinite(info.tauStd(ii))
                        lines{end+1} = sprintf( ...
                            '  tau_%d = %10.4g +/- %.2g %s%s', ...
                            ii, tauOut(ii), info.tauStd(ii), ...
                            app.tUnitTxt(), tag); %#ok<AGROW>
                    else
                        lines{end+1} = sprintf('  tau_%d = %10.4g %s%s', ...
                            ii, tauOut(ii), app.tUnitTxt(), tag); %#ok<AGROW>
                    end
                end
                if app.kfHasInf
                    lines{end+1} = '  tau_inf  (offset)';
                end
                lines{end+1} = '';
                lines{end+1} = 'Amplitudes:';
                nMain = numel(tauOut);
                nTot  = nMain + double(app.kfHasInf);
                for ii = 1:nMain
                    lines{end+1} = sprintf('  A_%d   = %12.4g', ...
                        ii, lastFit.A(ii)); %#ok<AGROW>
                end
                if app.kfHasInf
                    lines{end+1} = sprintf('  A_inf = %12.4g', ...
                        lastFit.A(nTot));
                end
                lines{end+1} = '';
                if isfield(info, 't0Std') && isfinite(info.t0Std)
                    lines{end+1} = sprintf('IRF t0   = %.4g +/- %.2g %s', ...
                        t0Out, info.t0Std, app.tUnitTxt());
                else
                    lines{end+1} = sprintf('IRF t0   = %.4g %s%s', ...
                        t0Out, app.tUnitTxt(), ...
                        [ternary(app.kfT0Fixed, ' (fixed)', ''), ...
                         app.undetNote(info, 't0Std', 1)]);
                end
                if isfield(info, 'fwhmAtLimit') && info.fwhmAtLimit
                    lines{end+1} = sprintf('IRF FWHM = %.4g %s  [not resolved]', ...
                        fwhmOut, app.tUnitTxt());
                elseif isfield(info, 'fwhmStd') && isfinite(info.fwhmStd)
                    lines{end+1} = sprintf('IRF FWHM = %.4g +/- %.2g %s', ...
                        fwhmOut, info.fwhmStd, app.tUnitTxt());
                else
                    lines{end+1} = sprintf('IRF FWHM = %.4g %s%s', ...
                        fwhmOut, app.tUnitTxt(), ...
                        [ternary(app.kfFWHMFixed, ' (fixed)', ''), ...
                         app.undetNote(info, 'fwhmStd', 1)]);
                end
                lines = [lines, app.irfLimitNote(info)];
                lines{end+1} = ['IRF model: ' ternary( ...
                    strcmp(app.kfIRFTreat, 'none'), ...
                    'none (pure exponentials)', 'numerical conv.')];
                txtResult.Value = lines;
                lblStatus.Text  = sprintf('Done. RMS = %.3g', info.rms);
                
                % Persist so the fit survives closing this window.
                % Settings live in app.kf* already; only the result
                % and its report text need snapshotting.
                app.kfSaved = struct( ...
                    'lastFit', lastFit, ...
                    'status',  lblStatus.Text, ...
                    'lines',   {txtResult.Value});
            end
            
            function drawWithFit()
                app.claAll(axTrace);  hold(axTrace, 'on');
                plot(axTrace, lastFit.t, lastFit.data, 'k.', ...
                    'MarkerSize', 7, 'DisplayName', 'Data');
                plot(axTrace, lastFit.t, lastFit.fitVec, 'r-', ...
                    'LineWidth', 1.6, 'DisplayName', 'Fit');
                yline(axTrace, 0, 'k:', 'HandleVisibility', 'off');
                hold(axTrace, 'off');  grid(axTrace, 'on');
                xlabel(axTrace, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axTrace, '\DeltaA or Intensity');
                title(axTrace, sprintf( ...
                    'Fit at %s = %.2f %s  (RMS = %.3g)', ...
                    app.xSym(), lastFit.wlActual, app.xUnitAx(), lastFit.rms));
                legend(axTrace, 'Location', 'best');
                
                cla(axRes);  hold(axRes, 'on');
                plot(axRes, lastFit.t, lastFit.res, '-', ...
                    'Color', [0.45 0.45 0.45], 'LineWidth', 0.9);
                yline(axRes, 0, 'k:');
                hold(axRes, 'off');  grid(axRes, 'on');
                xlabel(axRes, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axRes, 'Residual (\DeltaA or Intensity)');
                title(axRes, sprintf('Residual (RMS = %.3g)', lastFit.rms));
                refreshAxesScale();
            end
            
            function doExport()
                if isempty(lastFit.fitVec)
                    uialert(f, 'Run a fit first.', 'No results');
                    return;
                end
                mat = [lastFit.t, lastFit.data, lastFit.fitVec, lastFit.res];
                hdr = {sprintf('delay_%s', app.tUnitHdr()), ...
                       'data_dA', 'fit_dA', 'residual_dA'};
                
                % Append a parameter summary as comment-like header rows
                % via a separate "params" file. To keep things simple,
                % save just the trace as CSV/DAT.
                suggested = sprintf('kinfit_%.1f%s.csv', lastFit.wlActual, app.xUnitTxt());
                [file, path] = app.putFile( ...
                    {'*.csv','CSV (*.csv)'; '*.dat','DAT (tab, *.dat)'}, ...
                    sprintf('Export kinetic fit @ %.2f %s', ...
                        lastFit.wlActual, app.xUnitTxt()), suggested);
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                fid = -1;
                try
                    T = array2table(mat, 'VariableNames', hdr);
                    writetable(T, full, 'Delimiter', dlm);
                    % Also write parameters next to it
                    [pathOut, baseOut] = fileparts(full);
                    paramPath = fullfile(pathOut, [baseOut '_params.txt']);
                    fid = fopen(paramPath, 'w');
                    if fid < 0
                        error('Cannot write %s', paramPath);
                    end
                    fprintf(fid, 'Single-trace kinetic fit\n');
                    fprintf(fid, '%s: %.4f %s\n', ...
                        ternary(strcmp(app.specUnit, 'cm-1'), ...
                            'Wavenumber', 'Wavelength'), ...
                        lastFit.wlActual, app.xUnitTxt());
                    fprintf(fid, 'RMS       : %.6g\n\n', lastFit.rms);
                    fprintf(fid, 'Time constants:\n');
                    for ii = 1:numel(lastFit.tau)
                        if ii <= numel(lastFit.stretchOn) && lastFit.stretchOn(ii)
                            fprintf(fid, '  tau_%d = %.6g %s   beta_%d = %.4f  [stretched]\n', ...
                                ii, lastFit.tau(ii), app.tUnitHdr(), ...
                                ii, lastFit.beta(ii));
                        else
                            fprintf(fid, '  tau_%d = %.6g %s\n', ...
                                ii, lastFit.tau(ii), app.tUnitHdr());
                        end
                    end
                    fprintf(fid, '\nAmplitudes:\n');
                    nMain = numel(lastFit.tau);
                    for ii = 1:nMain
                        fprintf(fid, '  A_%d   = %.6g\n', ...
                            ii, lastFit.A(ii));
                    end
                    if lastFit.hasInf
                        fprintf(fid, '  A_inf = %.6g\n', ...
                            lastFit.A(nMain + 1));
                    end
                    fprintf(fid, '\nIRF t0   = %.6g %s\n', ...
                        lastFit.t0,   app.tUnitHdr());
                    fprintf(fid, 'IRF FWHM = %.6g %s\n', ...
                        lastFit.fwhm, app.tUnitHdr());
                    fclose(fid);
                    fid = -1;
                    uialert(f, sprintf( ...
                        ['Saved trace to:\n  %s\n\nParameters in:\n  %s'], ...
                        full, paramPath), 'Exported', 'Icon', 'success');
                catch ME
                    % (do not leave the parameter file open)
                    if fid > 2
                        try, fclose(fid); catch, end
                    end
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Export error');
                end
            end
            
            function closeAndClear()
                delete(f);
                app.kinFitFig = [];
            end
            
            function s = ternary(cond, a, b)
                if cond, s = a; else, s = b; end
            end
        end
        
        % =================================================================
        % Crop
        % =================================================================
        function tf = isCropped(app)
            tf = ~isempty(app.originalWavelength) && ...
                 (numel(app.wavelength) ~= numel(app.originalWavelength) || ...
                  numel(app.delay)      ~= numel(app.originalDelay));
        end
        
        function applyCropByRange(app, wlMin, wlMax, tMin, tMax)
            % Crop the CURRENT working data and bake all previously-
            % applied corrections (BG, chirp, t=0 shift, masks) into
            % the new raw baseline. This way the chirp's NaN edge
            % truly disappears when the user crops it away -- if we
            % kept chirpApplied = true, the subsequent recompute()
            % would re-run interp1 on the cropped axes and produce a
            % fresh NaN edge inside the crop window, defeating the
            % whole point.
            %
            % The pristine raw snapshot (originalDeltaA + original*)
            % is left untouched so "Reset all" can still restore it.
            if isempty(app.deltaA), return; end
            if wlMin > wlMax, tmp = wlMin; wlMin = wlMax; wlMax = tmp; end
            if tMin  > tMax,  tmp = tMin;  tMin  = tMax;  tMax  = tmp;  end
            
            wlMask = app.wavelength >= wlMin & ...
                     app.wavelength <= wlMax;
            tMask  = app.delay      >= tMin  & ...
                     app.delay      <= tMax;
            
            if ~any(wlMask) || ~any(tMask)
                uialert(app.fig, 'Empty crop range.', 'Invalid range');
                return;
            end
            
            if nnz(wlMask) < 2 || nnz(tMask) < 2
                uialert(app.fig, ['The crop range must keep at least 2 ' ...
                    'points on each axis.'], 'Invalid range');
                return;
            end
            
            % Previous analysis results were computed on the pre-crop
            % axes and no longer apply.
            app.dataChanged();
            
            % Slice axes and the displayed (corrected) data, then
            % freeze that result as the new baseline.
            app.wavelength = app.wavelength(wlMask);
            app.delay      = app.delay(tMask);
            cropped        = app.deltaA(wlMask, tMask);
            app.deltaA     = cropped;
            app.deltaA_raw = cropped;
            
            % All prior correction state is now baked into deltaA_raw;
            % clear the flags so recompute() will not re-apply them.
            % Any future correction the user adds (a new chirp pass,
            % a new background, etc.) will stack on top of the cropped
            % + corrected baseline.
            app.bgApplied    = false;  app.bgSpectrum = [];  app.bgN = 0;
            app.chirpApplied = false;
            app.chirpPts     = zeros(0,2);
            app.chirpParams  = [];     app.chirpFitRMS = NaN;
            app.ridgeWL      = [];     app.ridgeT      = [];
            app.maskedRegions = cell(0, 3);
            app.tZeroShift   = 0;
            
            % Clip overlays and selection to new range
            wlLo = min(app.wavelength);  wlHi = max(app.wavelength);
            tLo  = min(app.delay);       tHi  = max(app.delay);
            app.specOverlays = app.specOverlays( ...
                app.specOverlays >= tLo & app.specOverlays <= tHi);
            app.kinOverlays  = app.kinOverlays( ...
                app.kinOverlays  >= wlLo & app.kinOverlays  <= wlHi);
            app.selWL = min(max(app.selWL, wlLo), wlHi);
            app.selT  = min(max(app.selT,  tLo),  tHi);
            
            app.edSelWL.Limits = [wlLo, wlHi];
            app.edSelT.Limits  = [tLo,  tHi ];
            app.edSelWL.Value  = app.selWL;
            app.edSelT.Value   = app.selT;
            
            % Reset main display range to new full range
            app.mainViewTMin = tLo;
            app.mainViewTMax = tHi;
            if ~isempty(app.edViewTMin) && isvalid(app.edViewTMin)
                app.edViewTMin.Value = app.mainViewTMin;
                app.edViewTMax.Value = app.mainViewTMax;
            end
            
            % Reset Z range to auto based on cropped data
            app.mapZMin = [];  app.mapZMax = [];
            if ~isempty(app.edZMin) && isvalid(app.edZMin)
                cl = max(abs(app.deltaA_raw(:)), [], 'omitnan');
                if isempty(cl) || cl == 0 || isnan(cl), cl = 1; end
                app.edZMin.Value = -cl;
                app.edZMax.Value =  cl;
            end
            
            app.updateStatusLabel();
            app.updateAll();
        end
        
        % =================================================================
        function openSvdWindow(app)
            % SVD analysis: singular-value spectrum, spectral and
            % temporal basis vectors, and N-component reconstruction.
            % NaN entries (from wavelength masks) are replaced with 0
            % before the decomposition (option (b) policy).
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');  return;
            end
            if ~isempty(app.svdFig) && isvalid(app.svdFig)
                figure(app.svdFig);  return;
            end
            
            f = uifigure('Name', 'SVD Analysis', ...
                'Position', app.fitFigurePosition(1280, 820));
            app.applyLightTheme(f);
            app.svdFig = f;
            f.CloseRequestFcn = @(~,~) onClose();
            
            gl = uigridlayout(f, [1, 2]);
            gl.ColumnWidth = {320, '1x'};
            gl.Padding = [10 10 10 10];  gl.ColumnSpacing = 10;
            
            % ---------- Left panel: controls + singular values plot ----------
            left = uipanel(gl, 'Title', 'SVD setup');
            left.Layout.Column = 1;
            lgl = uigridlayout(left, [6, 1]);
            lgl.RowHeight = {'fit', 'fit', 'fit', 'fit', '1x', 'fit'};
            lgl.Padding = [10 10 10 10];  lgl.RowSpacing = 8;
            
            % nKeep selector (number of components for reconstruction)
            rowN = uigridlayout(lgl, [1, 3]);
            rowN.ColumnWidth = {'fit', 90, '1x'};
            rowN.Padding = [0 0 0 0];  rowN.ColumnSpacing = 6;
            uilabel(rowN, 'Text', 'Components to keep (N):', ...
                'HorizontalAlignment', 'right');
            edN = uieditfield(rowN, 'numeric', 'Value', 3, ...
                'Limits', [1 50], 'RoundFractionalValues', 'on', ...
                'ValueChangedFcn', @(~,~) refreshAll());
            uilabel(rowN, 'Text', '');
            
            % nShow selector (how many basis vectors to plot)
            rowS = uigridlayout(lgl, [1, 3]);
            rowS.ColumnWidth = {'fit', 90, '1x'};
            rowS.Padding = [0 0 0 0];  rowS.ColumnSpacing = 6;
            uilabel(rowS, 'Text', 'Basis vectors to show:', ...
                'HorizontalAlignment', 'right');
            edShow = uieditfield(rowS, 'numeric', 'Value', 3, ...
                'Limits', [1 12], 'RoundFractionalValues', 'on', ...
                'ValueChangedFcn', @(~,~) refreshAll());
            uilabel(rowS, 'Text', '');
            
            % Temporal axis scale (Linear or Log)
            rowT = uigridlayout(lgl, [1, 3]);
            rowT.ColumnWidth = {'fit', 90, '1x'};
            rowT.Padding = [0 0 0 0];  rowT.ColumnSpacing = 6;
            uilabel(rowT, 'Text', 'Temporal axis scale:', ...
                'HorizontalAlignment', 'right');
            ddTScale = uidropdown(rowT, ...
                'Items', {'Linear', 'Log'}, ...
                'Value', 'Log', ...
                'ValueChangedFcn', @(~,~) refreshAll());
            uilabel(rowT, 'Text', '');
            
            % Status / info text
            infoLbl = uilabel(lgl, 'Text', '', ...
                'WordWrap', 'on', 'FontColor', [0.3 0.3 0.3]);
            
            % Singular value plot (top-N)
            sigPanel = uipanel(lgl, 'Title', 'Singular values (log scale)');
            sigPanel.Layout.Row = 5;
            ggSig = uigridlayout(sigPanel, [1, 1]);
            ggSig.Padding = [4 4 4 4];
            axSig = uiaxes(ggSig);
            axSig.YScale = 'log';  grid(axSig, 'on');
            xlabel(axSig, 'Component index');
            ylabel(axSig, '\sigma_i');
            
            % Export buttons
            rowE = uigridlayout(lgl, [1, 3]);
            rowE.ColumnWidth = {'fit', 'fit', 'fit'};
            rowE.Padding = [0 0 0 0];  rowE.ColumnSpacing = 6;
            uibutton(rowE, 'Text', 'Export U / S / V', ...
                'ButtonPushedFcn', @(~,~) exportAll());
            uibutton(rowE, 'Text', 'Export reconstruction', ...
                'ButtonPushedFcn', @(~,~) exportReconstruction());
            uilabel(rowE, 'Text', '');
            
            % ---------- Right panel: plots ----------
            right = uipanel(gl, 'Title', 'Bases and reconstruction');
            right.Layout.Column = 2;
            rgl = uigridlayout(right, [3, 2]);
            rgl.RowHeight   = {'1x', '1x', '1.2x'};
            rgl.ColumnWidth = {'1x', '1x'};
            rgl.Padding = [10 10 10 10];
            rgl.RowSpacing = 8;  rgl.ColumnSpacing = 10;
            
            pnlU = uipanel(rgl, 'Title', ...
                'Spectral basis  U(:, k)  vs wavelength');
            pnlU.Layout.Row = 1;  pnlU.Layout.Column = 1;
            ggU = uigridlayout(pnlU, [1, 1]);
            ggU.Padding = [4 4 4 4];
            axU = uiaxes(ggU);
            grid(axU, 'on');
            xlabel(axU, app.xAxLabel());  ylabel(axU, 'U_k');
            
            pnlV = uipanel(rgl, 'Title', ...
                'Temporal basis  V(:, k)  vs delay');
            pnlV.Layout.Row = 1;  pnlV.Layout.Column = 2;
            ggV = uigridlayout(pnlV, [1, 1]);
            ggV.Padding = [4 4 4 4];
            axV = uiaxes(ggV);
            grid(axV, 'on');
            xlabel(axV, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axV, 'V_k');
            
            pnlExp = uipanel(rgl, 'Title', 'Original data');
            pnlExp.Layout.Row = 2;  pnlExp.Layout.Column = 1;
            ggE = uigridlayout(pnlExp, [1, 1]);
            ggE.Padding = [4 4 4 4];
            axExp = uiaxes(ggE);
            
            pnlRec = uipanel(rgl, 'Title', 'Reconstruction (N components)');
            pnlRec.Layout.Row = 2;  pnlRec.Layout.Column = 2;
            ggR = uigridlayout(pnlRec, [1, 1]);
            ggR.Padding = [4 4 4 4];
            axRec = uiaxes(ggR);
            
            pnlRes = uipanel(rgl, 'Title', 'Residual = original - reconstruction');
            pnlRes.Layout.Row = 3;  pnlRes.Layout.Column = [1 2];
            ggRs = uigridlayout(pnlRes, [1, 1]);
            ggRs.Padding = [4 4 4 4];
            axRes = uiaxes(ggRs);
            
            % ---------- Compute SVD once on open; refresh views ----------
            % Replace NaN with 0 (option-b policy), keep a mask for info.
            D = app.deltaA;            % rows = wavelength, cols = delay
            wl = app.wavelength(:);
            t  = app.delay(:).';
            nanMask = ~isfinite(D);
            nNaN = nnz(nanMask);
            Dz = D;  Dz(nanMask) = 0;
            
            try
                [U, S, V] = svd(Dz, 'econ');
            catch ME
                uialert(f, sprintf( ...
                    ['SVD failed:\n%s\n\nTry reducing data size ' ...
                     '(crop) or check for non-finite values.'], ...
                    ME.message), 'SVD error');
                return;
            end
            sigma = diag(S);
            nMax  = numel(sigma);
            
            % ---- Restore previous view settings (if any) ----
            % SVD itself recomputes instantly on open, so only the
            % three view settings need to survive a close / reopen.
            if ~isempty(app.svdSaved) && isstruct(app.svdSaved)
                try
                    edN.Value      = min(app.svdSaved.N,    nMax);
                    edShow.Value   = min(app.svdSaved.show, nMax);
                    ddTScale.Value = app.svdSaved.tscale;
                catch
                end
            end
            
            refreshAll();
            
            % =====================  nested helpers  =====================
            
            function refreshAll()
                if ~isvalid(f), return; end
                nKeep = min(round(edN.Value), nMax);
                nShow = min(round(edShow.Value), nMax);
                
                % Snapshot the view settings so a reopened window
                % starts where the user left off.
                app.svdSaved = struct('N', nKeep, 'show', nShow, ...
                    'tscale', ddTScale.Value);
                
                % Singular value plot
                cla(axSig);
                kAll = 1:min(nMax, max(20, 2*nShow));
                semilogy(axSig, kAll, sigma(kAll), 'o-', ...
                    'MarkerSize', 5, 'LineWidth', 1.2);
                hold(axSig, 'on');
                semilogy(axSig, 1:nShow, sigma(1:nShow), 'o', ...
                    'MarkerFaceColor', [0.20 0.45 0.80], ...
                    'MarkerEdgeColor', [0.20 0.45 0.80], 'MarkerSize', 6);
                xline(axSig, nKeep + 0.5, 'r--', ...
                    sprintf('N=%d', nKeep), ...
                    'LabelOrientation', 'horizontal');
                hold(axSig, 'off');
                grid(axSig, 'on');
                xlabel(axSig, 'Component index');
                ylabel(axSig, '\sigma_i');
                
                % Spectral basis vectors (U columns)
                app.claAll(axU);  hold(axU, 'on');
                cmapLines = lines(max(nShow, 6));
                for k = 1:nShow
                    plot(axU, wl, U(:, k), '-', 'LineWidth', 1.2, ...
                        'Color', cmapLines(k, :), ...
                        'DisplayName', sprintf('k=%d  (\\sigma=%.3g)', ...
                            k, sigma(k)));
                end
                yline(axU, 0, 'k:', 'HandleVisibility', 'off');
                hold(axU, 'off');
                grid(axU, 'on');
                xlabel(axU, app.xAxLabel());  ylabel(axU, 'U_k');
                xlim(axU, [min(wl), max(wl)]);
                legend(axU, 'Location', 'best', 'FontSize', 8);
                
                % Temporal basis vectors (V columns) with chosen scale
                app.claAll(axV);  hold(axV, 'on');
                tScaleStr = lower(ddTScale.Value);
                for k = 1:nShow
                    plot(axV, t, V(:, k), '-', 'LineWidth', 1.2, ...
                        'Color', cmapLines(k, :), ...
                        'DisplayName', sprintf('k=%d', k));
                end
                yline(axV, 0, 'k:', 'HandleVisibility', 'off');
                hold(axV, 'off');
                grid(axV, 'on');
                xlabel(axV, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axV, 'V_k');
                set(axV, 'YScale', 'linear');   % V_k is signed; never log
                applyXScale(axV, t, tScaleStr);
                legend(axV, 'Location', 'best', 'FontSize', 8);
                
                % Reconstruction with N components and residual
                Dn   = U(:, 1:nKeep) * diag(sigma(1:nKeep)) * V(:, 1:nKeep)';
                Rn   = Dz - Dn;
                cl   = max(abs(Dz(:)), [], 'omitnan');
                if isempty(cl) || cl == 0 || ~isfinite(cl), cl = 1; end
                clr  = max(abs(Rn(:)), [], 'omitnan');
                if isempty(clr) || clr == 0 || ~isfinite(clr), clr = cl; end
                
                drawHeatmap(axExp, t, wl, Dz,  [-cl  cl ], tScaleStr);
                drawHeatmap(axRec, t, wl, Dn,  [-cl  cl ], tScaleStr);
                drawHeatmap(axRes, t, wl, Rn,  [-clr clr], tScaleStr);
                
                % Info line
                varExpl = sum(sigma(1:nKeep).^2) / sum(sigma.^2) * 100;
                infoTxt = sprintf( ...
                    ['Data: %d wavelengths x %d delays\n' ...
                     'NaN entries replaced with 0: %d\n' ...
                     'N=%d components -> %.4f %% variance explained\n' ...
                     '\\sigma_1 / \\sigma_N = %.3g / %.3g\n' ...
                     'Tip: use the knee of the \\sigma plot as the ' ...
                     'starting Components value in Global Analysis.'], ...
                    size(D, 1), size(D, 2), nNaN, ...
                    nKeep, varExpl, sigma(1), sigma(nKeep));
                infoLbl.Text = app.uiStr(infoTxt);
            end
            
            function drawHeatmap(ax, tAxis, wlAxis, C, clim, scaleStr)
                % C is wavelength x delay (rows x cols).
                cla(ax);
                [W, T] = meshgrid(wlAxis, tAxis);
                surface(ax, W, T, zeros(size(C.')), C.', ...
                    'EdgeColor', 'none', 'FaceColor', 'flat');
                view(ax, 2);  axis(ax, 'tight');
                set(ax, 'YDir', 'normal');
                ax.CLim = clim;
                colormap(ax, app.getColormap());
                xlabel(ax, app.xAxLabel());
                ylabel(ax, sprintf('Delay time (%s)', app.tUnitAx()));
                xlim(ax, [min(wlAxis), max(wlAxis)]);
                applyYScale(ax, tAxis, scaleStr);
                colorbar(ax);
            end
            
            function applyXScale(ax, tAxis, scaleStr)
                if strcmp(scaleStr, 'log')
                    pos = tAxis(tAxis > 0);
                    if isempty(pos)
                        set(ax, 'XScale', 'linear');
                        xlim(ax, [min(tAxis), max(tAxis)]);
                    else
                        set(ax, 'XScale', 'log');
                        xlim(ax, [min(pos), max(tAxis)]);
                    end
                else
                    set(ax, 'XScale', 'linear');
                    xlim(ax, [min(tAxis), max(tAxis)]);
                end
            end
            
            function applyYScale(ax, tAxis, scaleStr)
                if strcmp(scaleStr, 'log')
                    pos = tAxis(tAxis > 0);
                    if isempty(pos)
                        set(ax, 'YScale', 'linear');
                        ylim(ax, [min(tAxis), max(tAxis)]);
                    else
                        set(ax, 'YScale', 'log');
                        ylim(ax, [min(pos), max(tAxis)]);
                    end
                else
                    set(ax, 'YScale', 'linear');
                    ylim(ax, [min(tAxis), max(tAxis)]);
                end
            end
            
            function exportAll()
                [fn, fp] = app.putFile({'*.csv', 'CSV files'}, ...
                    'Choose base filename for U / S / V exports', ...
                    'svd_result.csv');
                if isequal(fn, 0), return; end
                [~, base, ~] = fileparts(fn);
                outDir = fp;
                nKeep = min(round(edN.Value), nMax);
                
                try
                    % U: wavelength x nKeep
                    writematrix([NaN, 1:nKeep; wl, U(:, 1:nKeep)], ...
                        fullfile(outDir, [base '_U.csv']));
                    % V: delay x nKeep
                    writematrix([NaN, 1:nKeep; t.', V(:, 1:nKeep)], ...
                        fullfile(outDir, [base '_V.csv']));
                    % sigma: full vector
                    writematrix([(1:nMax).', sigma], ...
                        fullfile(outDir, [base '_sigma.csv']));
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Export error');
                    return;
                end
                
                uialert(f, sprintf(['Saved:\n  %s_U.csv\n  %s_V.csv\n  ' ...
                    '%s_sigma.csv\n\nFolder:\n%s'], ...
                    base, base, base, outDir), 'Export complete', ...
                    'Icon', 'success');
            end
            
            function exportReconstruction()
                [fn, fp] = app.putFile({'*.csv', 'CSV files'}, ...
                    'Save reconstruction', 'svd_reconstruction.csv');
                if isequal(fn, 0), return; end
                nKeep = min(round(edN.Value), nMax);
                Dn = U(:, 1:nKeep) * diag(sigma(1:nKeep)) * V(:, 1:nKeep)';
                app.writeDataFile(fullfile(fp, fn), wl, t, Dn, ',');
                uialert(f, ['Saved: ' fullfile(fp, fn)], ...
                    'Export complete', 'Icon', 'success');
            end
            
            function onClose()
                app.svdFig = [];
                delete(f);
            end
        end
        
        function openLdaWindow(app)
            %OPENLDAWINDOW  Lifetime Density Analysis window.
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');
                return;
            end
            if ~isempty(app.ldaFig) && isvalid(app.ldaFig)
                figure(app.ldaFig);  return;
            end
            
            f = uifigure('Name', 'Lifetime Density Analysis', ...
                'Position', app.fitFigurePosition(1380, 900));
            app.applyLightTheme(f);
            app.ldaFig = f;
            f.CloseRequestFcn = @(~,~) onClose();
            % Drag motion/up handlers are permanently attached and
            % no-op until "dragging" becomes true on a map click.
            f.WindowButtonMotionFcn = @(~,~) onLdaMotion();
            f.WindowButtonUpFcn     = @(~,~) onLdaUp();
            
            gl = uigridlayout(f, [1, 2]);
            gl.ColumnWidth = {360, '1x'};
            gl.Padding = [10 10 10 10];  gl.ColumnSpacing = 10;
            
            % ============== Left: setup ==============
            left = uipanel(gl, 'Title', 'LDA setup');
            left.Layout.Column = 1;
            lgl = uigridlayout(left, [17, 1]);
            % rows: tauLabel, tauMin, tauMax, K, irfLabel, t0, FWHM,
            %       noIRF, fitTRange, regLabel, regType, alpha, buttons,
            %       expLabel, expBtns, statusLbl, resultArea.  17 rows.
            lgl.RowHeight = {22, 30, 30, 30, 22, 30, 30, 26, 30, ...
                             22, 30, 30, 32, 22, 30, 22, '1x'};
            lgl.Padding = [8 8 8 8];  lgl.RowSpacing = 6;
            
            % --- Tau grid ---
            uilabel(lgl, 'Text', 'Lifetime grid (log-spaced):', ...
                'FontWeight', 'bold');
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', app.uiStr(sprintf('\\tau_{min} (%s):', app.tUnitTxt())));
            edTauMin = uieditfield(row, 'numeric', ...
                'Limits', [eps Inf], ...
                'Value', max(0.5 * app.gaFWHM, eps));
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', app.uiStr(sprintf('\\tau_{max} (%s):', app.tUnitTxt())));
            tspan = max(app.delay) - min(app.delay);
            if tspan <= 0, tspan = 1000; end
            tauMaxDef = 3 * max(app.delay);
            if ~isfinite(tauMaxDef) || tauMaxDef <= 0, tauMaxDef = 5 * tspan; end
            edTauMax = uieditfield(row, 'numeric', ...
                'Limits', [eps Inf], 'Value', tauMaxDef);
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Grid points (K):');
            edK = uieditfield(row, 'numeric', 'Value', 100, ...
                'Limits', [10 500], 'RoundFractionalValues', 'on');
            
            % --- IRF ---
            uilabel(lgl, 'Text', 'IRF parameters:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('t_0 (%s):', app.tUnitTxt()));
            if ~isempty(app.gaResultT0), t0Init = app.gaResultT0; else, t0Init = app.gaT0; end
            edT0 = uieditfield(row, 'numeric', 'Value', t0Init);
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('FWHM (%s):', app.tUnitTxt()));
            if ~isempty(app.gaResultFWHM), fwInit = app.gaResultFWHM; else, fwInit = app.gaFWHM; end
            edFW = uieditfield(row, 'numeric', ...
                'Limits', [eps Inf], 'Value', fwInit);
            
            cbNoIRF = uicheckbox(lgl, 'Text', 'No IRF (pure exponentials)', ...
                'Value', false, ...
                'Tooltip', ['Basis becomes exp(-(t-t0)/tau)*H(t-t0) ' ...
                    'without IRF convolution; FWHM is ignored. ' ...
                    'Intended for data (or a Fit t range) starting ' ...
                    'after the IRF.'], ...
                'ValueChangedFcn', @(~,~) onNoIRFChange());
            
            row = uigridlayout(lgl, [1, 5]);
            row.ColumnWidth = {'fit', 58, 'fit', 58, 'fit'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 4;
            uilabel(row, 'Text', sprintf('Fit t range (%s):', ...
                app.tUnitTxt()));
            edFitTmin = uieditfield(row, 'numeric', ...
                'Value', min(app.delay), ...
                'Tooltip', ['Delay points outside this range are ' ...
                    'excluded from the LDA solve (reconstruction is ' ...
                    'still shown on the full axis).' newline ...
                    'NOTE: when restricting to t >= T, tau values ' ...
                    'much shorter than T are no longer constrained ' ...
                    'by the data -- raise tau_min accordingly and do ' ...
                    'not interpret the fast-tau part of the map.']);
            uilabel(row, 'Text', '~', 'HorizontalAlignment', 'center');
            edFitTmax = uieditfield(row, 'numeric', ...
                'Value', max(app.delay), ...
                'Tooltip', 'Upper bound of the fitted delay range.');
            uibutton(row, 'Text', 'Full', ...
                'Tooltip', 'Reset to the full delay range.', ...
                'ButtonPushedFcn', @(~,~) setFitTFull());
            
            % --- Regularization ---
            uilabel(lgl, 'Text', 'Regularization:', 'FontWeight', 'bold');
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Type:');
            ddReg = uidropdown(row, ...
                'Items', app.uiStr({'L2 (||A||^2)', ...
                    '1st-deriv (smooth in log-\tau)'}), ...
                'ItemsData', {'l2', 'l2deriv'}, ...
                'Value', 'l2deriv');
            
            row = uigridlayout(lgl, [1, 3]);
            row.ColumnWidth = {'fit', 100, '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', app.uiStr('\alpha:'));
            edAlpha = uieditfield(row, 'numeric', ...
                'Value', 1e-2, 'Limits', [0 Inf]);
            uilabel(row, 'Text', '');
            
            % --- Buttons ---
            btnRow = uigridlayout(lgl, [1, 3]);
            btnRow.ColumnWidth = {'1x', '1x', '1x'};
            btnRow.Padding = [0 0 0 0];  btnRow.ColumnSpacing = 6;
            btnRun = uibutton(btnRow, 'Text', 'Run LDA', ...
                'BackgroundColor', [0.20 0.55 0.40], ...
                'FontColor', 'w', 'FontWeight', 'bold');
            btnLcurve = uibutton(btnRow, 'Text', 'L-curve...', ...
                'Tooltip', ['Sweep alpha from 1e-6 to 1e2 and find ' ...
                    'the corner of the L-curve. The corner alpha ' ...
                    'is then loaded into the alpha field.']);
            btnReset = uibutton(btnRow, 'Text', 'Reset');
            
            % --- Export ---
            uilabel(lgl, 'Text', 'Export:', 'FontWeight', 'bold');
            expRow = uigridlayout(lgl, [1, 3]);
            expRow.ColumnWidth = {'1x', '1x', '1x'};
            expRow.Padding = [0 0 0 0];  expRow.ColumnSpacing = 6;
            btnExpAmap = uibutton(expRow, 'Text', app.uiStr('A(\lambda,\tau)'), ...
                'Enable', 'off');
            btnExpRec  = uibutton(expRow, 'Text', 'Recon.', ...
                'Enable', 'off');
            btnExp1D   = uibutton(expRow, 'Text', app.uiStr('A(\lambda*)'), ...
                'Enable', 'off');
            
            % --- Status ---
            lblStatus = uilabel(lgl, 'Text', 'Ready.', ...
                'FontColor', [0.2 0.2 0.2]);
            
            txtResult = uitextarea(lgl, 'Editable', 'off', ...
                'Value', {'Click "Run LDA" to compute.'}, ...
                'FontName', 'Monospaced');
            
            % ============== Right: plots ==============
            right = uipanel(gl, 'Title', 'Results');
            right.Layout.Column = 2;
            rgl = uigridlayout(right, [3, 2]);
            rgl.ColumnWidth = {'1x', '1x'};
            rgl.RowHeight   = {'1.5x', '1x', '1x'};
            rgl.Padding     = [8 8 8 8];  rgl.RowSpacing = 8;
            rgl.ColumnSpacing = 10;
            
            % Row 1: heatmap + 1D A(tau) at fixed lambda
            axMap = uiaxes(rgl);
            axMap.Layout.Row = 1;  axMap.Layout.Column = 1;
            axMap.YScale = 'log';
            xlabel(axMap, app.xAxLabel());
            ylabel(axMap, sprintf('Lifetime \\tau (%s)', app.tUnitAx()));
            title(axMap, 'Amplitude density A(\lambda, \tau)');
            
            ax1D = uiaxes(rgl);
            ax1D.Layout.Row = 1;  ax1D.Layout.Column = 2;
            xlabel(ax1D, sprintf('Lifetime \\tau (%s)', app.tUnitAx()));
            ylabel(ax1D, 'Amplitude');
            title(ax1D, 'Click on the map to see A(\tau)');
            grid(ax1D, 'on');
            ax1D.XScale = 'log';
            
            % Row 2: spectrum slice (left) + reconstruction (right)
            axSpec = uiaxes(rgl);
            axSpec.Layout.Row = 2;  axSpec.Layout.Column = 1;
            xlabel(axSpec, app.xAxLabel());
            ylabel(axSpec, 'Amplitude');
            title(axSpec, 'Spectrum at clicked \tau');
            grid(axSpec, 'on');
            
            axKin = uiaxes(rgl);
            axKin.Layout.Row = 2;  axKin.Layout.Column = 2;
            xlabel(axKin, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axKin, '\DeltaA or Intensity');
            title(axKin, 'Reconstruction at clicked \lambda');
            grid(axKin, 'on');
            
            % Row 3: L-curve (full width)
            axL = uiaxes(rgl);
            axL.Layout.Row = 3;  axL.Layout.Column = [1 2];
            xlabel(axL, '||residual||');
            ylabel(axL, '||solution||');
            title(axL, 'L-curve (run "L-curve..." to populate)');
            grid(axL, 'on');
            axL.XScale = 'log';  axL.YScale = 'log';
            
            % ---- State ----
            ldaResult = struct('Amap', [], 'Drec', [], 'tauGrid', [], ...
                'alpha', NaN, 'regType', '', 'rms', NaN);
            curLambda = app.selWL;
            curTau    = NaN;       % set on first click
            xhV       = [];        % vertical crosshair (xline)
            xhH       = [];        % horizontal crosshair (yline)
            dragging  = false;     % LDA map drag state
            
            % ---- Wire up ----
            btnRun.ButtonPushedFcn    = @(~,~) safeCb(@doRun, 'Run LDA');
            btnLcurve.ButtonPushedFcn = @(~,~) safeCb(@doLcurve, 'L-curve');
            btnReset.ButtonPushedFcn  = @(~,~) safeCb(@doReset, 'Reset');
            btnExpAmap.ButtonPushedFcn = @(~,~) safeCb(@exportAmap, 'Export A-map');
            btnExpRec.ButtonPushedFcn  = @(~,~) safeCb(@exportRec, 'Export recon');
            btnExp1D.ButtonPushedFcn   = @(~,~) safeCb(@export1D, 'Export 1D');
            
            % ---- Restore previous analysis (if any) ----
            if ~isempty(app.ldaSaved) && isstruct(app.ldaSaved) && ...
                    isfield(app.ldaSaved, 'result') && ...
                    ~isempty(app.ldaSaved.result.Amap) && ...
                    size(app.ldaSaved.result.Amap, 1) == numel(app.wavelength)
                try
                    st = app.ldaSaved.settings;
                    edTauMin.Value = st.tauMin;  edTauMax.Value = st.tauMax;
                    edK.Value      = st.K;
                    edT0.Value     = st.t0;      edFW.Value     = st.fwhm;
                    ddReg.Value    = st.reg;     edAlpha.Value  = st.alpha;
                    if isfield(st, 'noIRF')
                        cbNoIRF.Value = st.noIRF;  onNoIRFChange();
                    end
                    if isfield(st, 'fitTmin') && isfinite(st.fitTmin)
                        edFitTmin.Value = max(min(app.delay), st.fitTmin);
                    end
                    if isfield(st, 'fitTmax') && isfinite(st.fitTmax)
                        edFitTmax.Value = min(max(app.delay), st.fitTmax);
                    end
                    ldaResult = app.ldaSaved.result;
                    curLambda = app.ldaSaved.curLambda;
                    curTau    = app.ldaSaved.curTau;
                    btnExpAmap.Enable = 'on';
                    btnExpRec.Enable  = 'on';
                    btnExp1D.Enable   = 'on';
                    plotResults();
                    lblStatus.Text  = sprintf('%s  [restored]', ...
                        app.ldaSaved.status);
                    txtResult.Value = app.ldaSaved.lines;
                catch
                    lblStatus.Text = 'Stored result could not be restored -- run again.';
                end
            end
            
            function safeCb(fcn, label)
                try
                    fcn();
                catch ME
                    uialert(f, sprintf('"%s" failed:\n%s\n\n(In %s line %d)', ...
                        label, ME.message, ...
                        ME.stack(1).name, ME.stack(1).line), ...
                        'Callback error');
                    lblStatus.Text = sprintf('Error in %s.', label);
                end
            end
            
            % ============== Nested callbacks ==============
            function tauGrid = buildGrid()
                tauGrid = logspace(log10(edTauMin.Value), ...
                                   log10(edTauMax.Value), edK.Value);
            end
            
            function doRun()
                if edTauMin.Value <= 0 || edTauMax.Value <= edTauMin.Value
                    uialert(f, 'Need 0 < tau_min < tau_max.', 'Invalid');
                    return;
                end
                if edFW.Value <= 0
                    uialert(f, 'IRF FWHM must be positive.', 'Invalid');
                    return;
                end
                tauGrid = buildGrid();
                lblStatus.Text = 'Running LDA...';  drawnow;
                try
                    if cbNoIRF.Value, irfM = 'none'; else, irfM = 'irf'; end
                    [Amap, Drec, info] = app.computeLDA( ...
                        app.deltaA, app.delay, tauGrid, ...
                        edT0.Value, edFW.Value, ...
                        edAlpha.Value, ddReg.Value, ...
                        irfM, [edFitTmin.Value, edFitTmax.Value]);
                catch ME
                    uialert(f, sprintf('LDA failed:\n%s', ME.message), ...
                        'LDA error');
                    lblStatus.Text = 'LDA failed.';
                    return;
                end
                ldaResult.Amap    = Amap;
                ldaResult.Drec    = Drec;
                ldaResult.tauGrid = tauGrid;
                ldaResult.alpha   = edAlpha.Value;
                ldaResult.regType = ddReg.Value;
                ldaResult.rms     = info.rms;
                
                btnExpAmap.Enable = 'on';
                btnExpRec.Enable  = 'on';
                btnExp1D.Enable   = 'on';
                
                lblStatus.Text = sprintf('Done. RMS = %.3g (rel. %.2g%%)', ...
                    info.rms, 100 * info.ratio);
                txtResult.Value = { ...
                    sprintf('K = %d log-spaced taus', info.K), ...
                    sprintf('tau range: %.3g .. %.3g %s', ...
                        tauGrid(1), tauGrid(end), app.tUnitTxt()), ...
                    sprintf('alpha = %.3g  (%s)', edAlpha.Value, ddReg.Value), ...
                    sprintf('RMS residual: %.4g', info.rms), ...
                    sprintf('Relative residual: %.3f%%', 100 * info.ratio), ...
                    '', ...
                    'Click on the heatmap to inspect a particular wavelength.'};
                
                plotResults();
                
                % Persist so the analysis survives closing this window.
                app.ldaSaved = struct( ...
                    'result',  ldaResult, ...
                    'curLambda', curLambda, 'curTau', curTau, ...
                    'status',  lblStatus.Text, ...
                    'lines',   {txtResult.Value}, ...
                    'settings', struct( ...
                        'tauMin', edTauMin.Value, 'tauMax', edTauMax.Value, ...
                        'K',      edK.Value, ...
                        't0',     edT0.Value,     'fwhm',   edFW.Value, ...
                        'reg',    ddReg.Value,    'alpha',  edAlpha.Value, ...
                        'noIRF',  cbNoIRF.Value, ...
                        'fitTmin', edFitTmin.Value, ...
                        'fitTmax', edFitTmax.Value));
            end
            
            function plotResults()
                if isempty(ldaResult.Amap), return; end
                cmap = app.getColormap();
                cl = max(abs(ldaResult.Amap(:)), [], 'omitnan');
                if isempty(cl) || cl == 0, cl = 1; end
                
                cla(axMap);
                wl = app.wavelength(:).';
                tauG = ldaResult.tauGrid(:).';
                [W, T] = meshgrid(wl, tauG);
                surface(axMap, W, T, zeros(size(ldaResult.Amap.')), ...
                    ldaResult.Amap.', ...
                    'EdgeColor', 'none', 'FaceColor', 'flat');
                view(axMap, 2);
                axMap.YDir = 'normal';
                axMap.YScale = 'log';
                axMap.CLim = [-cl, cl];
                colormap(axMap, cmap);
                xlim(axMap, [min(wl), max(wl)]);
                ylim(axMap, [min(tauG), max(tauG)]);
                xlabel(axMap, app.xAxLabel());
                ylabel(axMap, sprintf('Lifetime \\tau (%s)', ...
                    app.tUnitAx()));
                title(axMap, sprintf( ...
                    ['Amplitude density A(\\lambda, \\tau)   ' ...
                     '\\alpha=%.3g'], ldaResult.alpha));
                cb = colorbar(axMap);
                cb.Label.String = '\DeltaA or Intensity  amplitude';
                
                % Crosshair (initialize to map center if no click yet).
                if ~isfinite(curLambda)
                    curLambda = mean([min(wl), max(wl)]);
                end
                if ~isfinite(curTau)
                    curTau = sqrt(min(tauG) * max(tauG));   % geometric mean
                end
                xhV = xline(axMap, curLambda, 'k-', 'LineWidth', 1.0, ...
                    'Alpha', 0.7, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                xhH = yline(axMap, curTau, 'k-', 'LineWidth', 1.0, ...
                    'Alpha', 0.7, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                
                % Click handler
                ch = axMap.Children;
                for jj = 1:numel(ch)
                    if isa(ch(jj), 'matlab.graphics.primitive.Surface') || ...
                       isa(ch(jj), 'matlab.graphics.primitive.Image')
                        ch(jj).ButtonDownFcn = @(~,ev) onMapClick(ev);
                    end
                end
                axMap.ButtonDownFcn = @(~,ev) onMapClick(ev);
                
                update1DandKin();
            end
            
            function onMapClick(ev)
                pt = ev.IntersectionPoint;
                if isempty(pt), return; end
                dragging = true;
                applyLdaPoint(pt);
            end
            function onLdaMotion()
                if ~dragging, return; end
                if ~app.mouseInsideFigure(f)
                    onLdaUp();  return;
                end
                cp = axMap.CurrentPoint;
                if isempty(cp) || size(cp,1) < 1, return; end
                applyLdaPoint(cp(1, 1:2));
            end
            function onLdaUp()
                dragging = false;
            end
            function applyLdaPoint(pt)
                if numel(pt) < 2, return; end
                wlR = [min(app.wavelength), max(app.wavelength)];
                tauR = [min(ldaResult.tauGrid), max(ldaResult.tauGrid)];
                curLambda = min(max(pt(1), wlR(1)),  wlR(2));
                curTau    = min(max(pt(2), tauR(1)), tauR(2));
                if ~isempty(xhV) && isvalid(xhV), xhV.Value = curLambda; end
                if ~isempty(xhH) && isvalid(xhH), xhH.Value = curTau;    end
                update1DandKin();
                drawnow limitrate;
            end
            
            function update1DandKin()
                if isempty(ldaResult.Amap), return; end
                [~, iw] = min(abs(app.wavelength - curLambda));
                wlAct   = app.wavelength(iw);
                
                % 1D distribution
                cla(ax1D);
                hold(ax1D, 'on');
                plot(ax1D, ldaResult.tauGrid, ldaResult.Amap(iw, :), ...
                    'b-', 'LineWidth', 1.4);
                yline(ax1D, 0, 'k:');
                if isfinite(curTau)
                    xline(ax1D, curTau, 'k--', 'LineWidth', 0.8, ...
                        'Alpha', 0.6);
                end
                hold(ax1D, 'off');
                xlabel(ax1D, sprintf('Lifetime \\tau (%s)', app.tUnitAx()));
                ylabel(ax1D, 'Amplitude');
                ax1D.XScale = 'log';
                xlim(ax1D, [min(ldaResult.tauGrid), max(ldaResult.tauGrid)]);
                title(ax1D, sprintf('A(\\tau) at %s = %.2f %s', ...
                    app.xSym(), wlAct, app.xUnitAx()));
                grid(ax1D, 'on');
                
                % Spectrum slice: A(:, k_tau) for the tau closest to curTau.
                cla(axSpec);
                if isfinite(curTau)
                    [~, kTau] = min(abs(ldaResult.tauGrid - curTau));
                    tauAct = ldaResult.tauGrid(kTau);
                    hold(axSpec, 'on');
                    plot(axSpec, app.wavelength, ldaResult.Amap(:, kTau), ...
                        'b-', 'LineWidth', 1.4);
                    yline(axSpec, 0, 'k:');
                    if isfinite(curLambda)
                        xline(axSpec, curLambda, 'k--', 'LineWidth', 0.8, ...
                            'Alpha', 0.6);
                    end
                    hold(axSpec, 'off');
                    xlim(axSpec, [min(app.wavelength), max(app.wavelength)]);
                    title(axSpec, sprintf( ...
                        'A(\\lambda) at \\tau = %.3g %s', ...
                        tauAct, app.tUnitTxt()));
                else
                    title(axSpec, 'Spectrum at clicked \tau');
                end
                xlabel(axSpec, app.xAxLabel());
                ylabel(axSpec, 'Amplitude');
                grid(axSpec, 'on');
                
                % Reconstruction overlay
                app.claAll(axKin);  hold(axKin, 'on');
                plot(axKin, app.delay, app.deltaA(iw, :), 'k.', ...
                    'MarkerSize', 6, 'DisplayName', 'Data');
                plot(axKin, app.delay, ldaResult.Drec(iw, :), 'r-', ...
                    'LineWidth', 1.4, 'DisplayName', 'LDA reconstr.');
                plot(axKin, app.delay, ...
                    app.deltaA(iw, :) - ldaResult.Drec(iw, :), '-', ...
                    'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
                    'DisplayName', 'Residual');
                yline(axKin, 0, 'k:', 'HandleVisibility', 'off');
                hold(axKin, 'off');  grid(axKin, 'on');
                xlabel(axKin, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axKin, '\DeltaA or Intensity');
                if any(app.delay > 0)
                    set(axKin, 'XScale', 'log');
                    xlK = [min(app.delay(app.delay > 0)), ...
                        max(app.delay)];
                    % Zoom onto the fitted sub-range when restricted.
                    frK = getFitTRangeIfRestricted();
                    if ~isempty(frK)
                        frK(1) = max(frK(1), xlK(1));
                        if frK(1) < frK(2), xlK = frK; end
                    end
                    xlim(axKin, xlK);
                end
                title(axKin, sprintf( ...
                    'Reconstruction at %s = %.2f %s', ...
                    app.xSym(), wlAct, app.xUnitAx()));
                legend(axKin, 'Location', 'best', 'FontSize', 8);
            end
            
            function doLcurve()
                tauGrid = buildGrid();
                lblStatus.Text = 'Computing L-curve (sweeping alpha)...';
                drawnow;
                try
                    if cbNoIRF.Value, irfM = 'none'; else, irfM = 'irf'; end
                    [aList, rN, sN, aCorner] = app.computeLcurve( ...
                        app.deltaA, app.delay, tauGrid, ...
                        edT0.Value, edFW.Value, ddReg.Value, ...
                        [], irfM, [edFitTmin.Value, edFitTmax.Value]);
                catch ME
                    uialert(f, sprintf('L-curve failed:\n%s', ME.message), ...
                        'Error');
                    lblStatus.Text = 'L-curve failed.';
                    return;
                end
                cla(axL);  hold(axL, 'on');
                plot(axL, rN, sN, 'b-o', 'MarkerSize', 4, ...
                    'LineWidth', 1.0);
                [~, iC] = min(abs(aList - aCorner));
                plot(axL, rN(iC), sN(iC), 'rs', 'MarkerSize', 12, ...
                    'LineWidth', 2.0);
                hold(axL, 'off');
                set(axL, 'XScale', 'log', 'YScale', 'log');
                grid(axL, 'on');
                xlabel(axL, '||residual||');
                ylabel(axL, '||solution||');
                title(axL, sprintf( ...
                    'L-curve  (corner \\alpha = %.3g)', aCorner));
                edAlpha.Value = aCorner;
                lblStatus.Text = sprintf( ...
                    'L-curve done. Suggested alpha = %.3g (loaded).', aCorner);
            end
            
            function doReset()
                edTauMin.Value = max(0.5 * app.gaFWHM, eps);
                tauMaxR = 3 * max(app.delay);
                if ~isfinite(tauMaxR) || tauMaxR <= 0
                    tauMaxR = 5 * tspan;
                end
                edTauMax.Value = tauMaxR;
                edK.Value      = 100;
                edT0.Value     = 0;
                edFW.Value     = app.gaFWHM;
                ddReg.Value    = 'l2deriv';
                edAlpha.Value  = 1e-2;
                ldaResult = struct('Amap', [], 'Drec', [], ...
                    'tauGrid', [], 'alpha', NaN, ...
                    'regType', '', 'rms', NaN);
                curTau = NaN;
                xhV = [];  xhH = [];
                btnExpAmap.Enable = 'off';
                btnExpRec.Enable  = 'off';
                btnExp1D.Enable   = 'off';
                cla(axMap);  cla(ax1D);  cla(axSpec);
                app.claAll(axKin);  cla(axL);
                title(axMap, 'Amplitude density A(\lambda, \tau)');
                title(ax1D, 'Click on the map to see A(\tau)');
                title(axSpec, 'Spectrum at clicked \tau');
                title(axKin, 'Reconstruction at clicked \lambda');
                title(axL, 'L-curve (run "L-curve..." to populate)');
                lblStatus.Text  = 'Reset.';
                txtResult.Value = {'Click "Run LDA" to compute.'};
                cbNoIRF.Value = false;
                onNoIRFChange();
                setFitTFull();
                app.ldaSaved = [];
            end
            
            function onNoIRFChange()
                % Pure-exp mode: FWHM has no effect.
                if cbNoIRF.Value
                    edFW.Enable = 'off';
                else
                    edFW.Enable = 'on';
                end
            end
            
            function setFitTFull()
                edFitTmin.Value = min(app.delay);
                edFitTmax.Value = max(app.delay);
            end
            
            function fr = getFitTRangeIfRestricted()
                fr = [];
                lo = min(edFitTmin.Value, edFitTmax.Value);
                hi = max(edFitTmin.Value, edFitTmax.Value);
                if ~isfinite(lo) || ~isfinite(hi), return; end
                dLo = min(app.delay);  dHi = max(app.delay);
                span = dHi - dLo;
                if span <= 0, return; end
                if (lo - dLo) > 1e-9 * span || (dHi - hi) > 1e-9 * span
                    fr = [lo, hi];
                end
            end
            
            function exportAmap()
                if isempty(ldaResult.Amap), return; end
                tauG = ldaResult.tauGrid;
                hdr = [{app.xHdr()}, ...
                    arrayfun(@(x) sprintf('tau=%.4g_%s', x, app.tUnitHdr()), ...
                        tauG, 'UniformOutput', false)];
                mat = [app.wavelength(:), ldaResult.Amap];
                doExportTable(mat, hdr, 'LDA_Amap', ...
                    'Export A(lambda, tau)');
            end
            
            function exportRec()
                if isempty(ldaResult.Drec), return; end
                t = app.delay(:).';
                hdr = [{app.xHdr()}, ...
                    arrayfun(@(x) sprintf('t=%.4g_%s', x, app.tUnitHdr()), ...
                        t, 'UniformOutput', false)];
                mat = [app.wavelength(:), ldaResult.Drec];
                doExportTable(mat, hdr, 'LDA_recon', ...
                    'Export reconstruction');
            end
            
            function export1D()
                if isempty(ldaResult.Amap), return; end
                [~, iw] = min(abs(app.wavelength - curLambda));
                hdr = {sprintf('tau_%s', app.tUnitHdr()), 'amplitude'};
                mat = [ldaResult.tauGrid(:), ldaResult.Amap(iw, :).'];
                doExportTable(mat, hdr, ...
                    sprintf('LDA_Atau_%.0f%s', app.wavelength(iw), app.xUnitTxt()), ...
                    'Export A(tau)');
            end
            
            function doExportTable(mat, hdr, defaultName, dlgTitle)
                [file, path] = app.putFile( ...
                    {'*.csv','CSV (*.csv)'; '*.dat','DAT (tab, *.dat)'}, ...
                    dlgTitle, [defaultName '.csv']);
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    T = array2table(mat, 'VariableNames', hdr);
                    writetable(T, full, 'Delimiter', dlm);
                    uialert(f, sprintf('Saved: %s', full), ...
                        'Exported', 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Error');
                end
            end
            
            function onClose()
                app.ldaFig = [];
                delete(f);
            end
        end
        
        function openMcrWindow(app)
            %OPENMCRWINDOW  Multivariate Curve Resolution - ALS window.
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');
                return;
            end
            if ~isempty(app.mcrFig) && isvalid(app.mcrFig)
                figure(app.mcrFig);  return;
            end
            
            f = uifigure('Name', 'MCR-ALS', ...
                'Position', app.fitFigurePosition(1380, 880));
            app.applyLightTheme(f);
            app.mcrFig = f;
            f.CloseRequestFcn = @(~,~) onClose();
            f.WindowButtonMotionFcn = @(~,~) onMcrMotion();
            f.WindowButtonUpFcn     = @(~,~) onMcrUp();
            
            gl = uigridlayout(f, [1, 2]);
            gl.ColumnWidth = {340, '1x'};
            gl.Padding = [10 10 10 10];  gl.ColumnSpacing = 10;
            
            % ============== Left: setup ==============
            left = uipanel(gl, 'Title', 'MCR-ALS setup');
            left.Layout.Column = 1;
            lgl = uigridlayout(left, [16, 1]);
            % rows: comp, "Constraints:" label, NNc, Uni, NNs, "Init:" label,
            %       initRow, seedRow, "Convergence:" label, convRow,
            %       timeScaleRow, btnRow,
            %       "Export:" label, expRow, statusLbl, resultArea(fills rest)
            lgl.RowHeight = {30, 22, 26, 26, 26, 22, 30, 30, 22, 30, 30, 32, ...
                             22, 30, 22, '1x'};
            lgl.Padding = [8 8 8 8];  lgl.RowSpacing = 6;
            
            % Components
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Components (N):');
            edN = uieditfield(row, 'numeric', 'Value', 3, ...
                'Limits', [1 30], 'RoundFractionalValues', 'on');
            
            % Constraints
            uilabel(lgl, 'Text', 'Constraints:', 'FontWeight', 'bold');
            cbNNc = uicheckbox(lgl, 'Text', ...
                'Non-negativity on C (concentrations)', 'Value', true);
            cbUni = uicheckbox(lgl, 'Text', ...
                'Unimodality on C (single peak per profile)', ...
                'Value', false);
            cbNNs = uicheckbox(lgl, 'Text', ...
                'Non-negativity on S (spectra) - usually OFF for TA', ...
                'Value', false);
            
            % Initialization
            uilabel(lgl, 'Text', 'Initialization:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Method:');
            ddInit = uidropdown(row, ...
                'Items', {'SVD-based (default)', 'Random', ...
                          'GA result: EADS / SAS', 'GA result: DADS'}, ...
                'ItemsData', {'svd', 'random', 'gaeads', 'ga'}, ...
                'Value', 'svd', ...
                'Tooltip', ['Where the iteration starts. MCR-ALS can ' ...
                    'have several solutions of the same quality, so ' ...
                    'the start can change the result: starts that ' ...
                    'agree indicate a well-determined solution.' newline ...
                    'SVD-based: the leading singular vectors of the ' ...
                    'data. No assumptions, the same result every run.' newline ...
                    'Random: random positive profiles. Every run takes ' ...
                    'a new seed (shown in the result box) unless one ' ...
                    'is entered below.' newline ...
                    'GA result: the spectra of the stored global fit ' ...
                    'as starting spectra (components the fit does not ' ...
                    'have are taken from the data). EADS / SAS start ' ...
                    'from the sequential picture -- populations that ' ...
                    'rise and decay -- and usually give species-like ' ...
                    'results. DADS start from the parallel picture ' ...
                    '(every component decays from time zero) and stay ' ...
                    'close to it.']);
            % Seed of the random start
            row = uigridlayout(lgl, [1, 3]);
            row.ColumnWidth = {'fit', 95, '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            lblSeed = uilabel(row, 'Text', 'Random seed:');
            edSeed = uieditfield(row, 'numeric', 'Value', 0, ...
                'Limits', [0 4294967295], ...
                'RoundFractionalValues', 'on', ...
                'ValueDisplayFormat', '%d', ...
                'Tooltip', ['Seed of the random start. 0: a new seed ' ...
                    'for every run; the one used is shown in the ' ...
                    'result box. Enter that number here to repeat a ' ...
                    'run exactly.']);
            lblSeedHint = uilabel(row, 'Text', '0 = new every run', ...
                'FontColor', [0.45 0.45 0.45]);
            
            % Convergence
            uilabel(lgl, 'Text', 'Convergence:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 4]);
            row.ColumnWidth = {'fit', 80, 'fit', 60};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Tol:');
            edTol = uieditfield(row, 'numeric', 'Value', 1e-4, ...
                'Limits', [eps Inf]);
            uilabel(row, 'Text', 'Max iter:');
            edMaxIt = uieditfield(row, 'numeric', 'Value', 200, ...
                'Limits', [1 10000], 'RoundFractionalValues', 'on');
            
            % Display options
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Time scale:');
            ddTimeScale = uidropdown(row, ...
                'Items', {'Log', 'Linear'}, ...
                'ItemsData', {'log', 'linear'}, ...
                'Value', 'log', ...
                'Tooltip', ['Delay axis scale for the C(t) profiles ' ...
                    'and the reconstruction-at-lambda panel.']);
            
            % Buttons
            btnRow = uigridlayout(lgl, [1, 2]);
            btnRow.ColumnWidth = {'1x', '1x'};
            btnRow.Padding = [0 0 0 0];  btnRow.ColumnSpacing = 6;
            btnRun = uibutton(btnRow, 'Text', 'Run MCR-ALS', ...
                'BackgroundColor', [0.85 0.50 0.20], ...
                'FontColor', 'w', 'FontWeight', 'bold');
            btnReset = uibutton(btnRow, 'Text', 'Reset');
            
            % Export
            uilabel(lgl, 'Text', 'Export:', 'FontWeight', 'bold');
            expRow = uigridlayout(lgl, [1, 4]);
            expRow.ColumnWidth = {'1x', '1x', '1x', '1x'};
            expRow.Padding = [0 0 0 0];  expRow.ColumnSpacing = 6;
            btnExpC = uibutton(expRow, 'Text', 'C', 'Enable', 'off');
            btnExpS = uibutton(expRow, 'Text', 'S', 'Enable', 'off');
            btnExpRec = uibutton(expRow, 'Text', 'Recon.', 'Enable', 'off');
            btnExpConv = uibutton(expRow, 'Text', 'Conv.', 'Enable', 'off');
            
            lblStatus = uilabel(lgl, 'Text', 'Ready.', ...
                'FontColor', [0.2 0.2 0.2]);
            txtResult = uitextarea(lgl, 'Editable', 'off', ...
                'Value', {'Click "Run MCR-ALS" to compute.'}, ...
                'FontName', 'Monospaced');
            
            % ============== Right: plots ==============
            right = uipanel(gl, 'Title', 'Results');
            right.Layout.Column = 2;
            rgl = uigridlayout(right, [2, 2]);
            rgl.ColumnWidth = {'1x', '1x'};
            rgl.RowHeight   = {'1x', '1x'};
            rgl.Padding     = [8 8 8 8];  rgl.RowSpacing = 8;
            rgl.ColumnSpacing = 10;
            
            axC = uiaxes(rgl);
            axC.Layout.Row = 1;  axC.Layout.Column = 1;
            xlabel(axC, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axC, 'Concentration');
            title(axC, 'C profiles (concentration vs time)');
            grid(axC, 'on');
            
            axS = uiaxes(rgl);
            axS.Layout.Row = 1;  axS.Layout.Column = 2;
            xlabel(axS, app.xAxLabel());
            ylabel(axS, 'Pure spectrum (\DeltaA or Intensity)');
            title(axS, 'S spectra');
            grid(axS, 'on');
            
            axRec = uiaxes(rgl);
            axRec.Layout.Row = 2;  axRec.Layout.Column = 1;
            xlabel(axRec, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axRec, '\DeltaA or Intensity');
            title(axRec, 'Reconstruction at clicked \lambda');
            grid(axRec, 'on');
            
            axConv = uiaxes(rgl);
            axConv.Layout.Row = 2;  axConv.Layout.Column = 2;
            xlabel(axConv, 'ALS iteration');
            ylabel(axConv, 'LOF (%)');
            title(axConv, 'Convergence');
            grid(axConv, 'on');
            axConv.YScale = 'log';
            
            % ---- State ----
            mcrResult = struct('C', [], 'S', [], 'Drec', [], ...
                'lof', NaN, 'iter', 0, 'lofHistory', [], ...
                'converged', false);
            curLambda = app.selWL;
            xhS       = [];   % vertical wavelength marker on axS
            draggingS = false; % MCR S-plot drag state
            
            % ---- Wire up ----
            btnRun.ButtonPushedFcn    = @(~,~) safeCb(@doRun, 'Run MCR-ALS');
            btnReset.ButtonPushedFcn  = @(~,~) safeCb(@doReset, 'Reset');
            ddTimeScale.ValueChangedFcn = @(~,~) onTimeScaleChange();
            ddInit.ValueChangedFcn    = @(~,~) onInitChange();
            
            % ---- Restore previous analysis (if any) ----
            if ~isempty(app.mcrSaved) && isstruct(app.mcrSaved) && ...
                    isfield(app.mcrSaved, 'result') && ...
                    ~isempty(app.mcrSaved.result.S) && ...
                    size(app.mcrSaved.result.S, 1) == numel(app.wavelength) && ...
                    size(app.mcrSaved.result.C, 1) == numel(app.delay)
                try
                    st = app.mcrSaved.settings;
                    edN.Value     = st.N;
                    cbNNc.Value   = st.nnC;   cbUni.Value = st.uniC;
                    cbNNs.Value   = st.nnS;
                    ddInit.Value  = st.init;
                    if isfield(st, 'seed'), edSeed.Value = st.seed; end
                    edTol.Value   = st.tol;   edMaxIt.Value = st.maxIt;
                    ddTimeScale.Value = st.tscale;
                    mcrResult = app.mcrSaved.result;
                    curLambda = app.mcrSaved.curLambda;
                    btnExpC.Enable    = 'on';
                    btnExpS.Enable    = 'on';
                    btnExpRec.Enable  = 'on';
                    btnExpConv.Enable = 'on';
                    plotResults();
                    lblStatus.Text  = sprintf('%s  [restored]', ...
                        app.mcrSaved.status);
                    txtResult.Value = app.mcrSaved.lines;
                catch
                    lblStatus.Text = 'Stored result could not be restored -- run again.';
                end
            end
            onInitChange();
            btnExpC.ButtonPushedFcn   = @(~,~) safeCb(@(~) exportMat( ...
                [app.delay(:), mcrResult.C], ...
                ['delay_' app.tUnitHdr(), ...
                 arrayfun(@(j) sprintf('C%d', j), 1:size(mcrResult.C, 2), ...
                    'UniformOutput', false)], ...
                'MCR_C', 'Export C profiles'), 'Export C');
            btnExpS.ButtonPushedFcn   = @(~,~) safeCb(@(~) exportMat( ...
                [app.wavelength(:), mcrResult.S], ...
                [{app.xHdr()}, ...
                 arrayfun(@(j) sprintf('S%d', j), 1:size(mcrResult.S, 2), ...
                    'UniformOutput', false)], ...
                'MCR_S', 'Export S spectra'), 'Export S');
            btnExpRec.ButtonPushedFcn = @(~,~) safeCb(@exportRec, 'Export recon');
            btnExpConv.ButtonPushedFcn = @(~,~) safeCb(@(~) exportMat( ...
                [(1:numel(mcrResult.lofHistory)).', mcrResult.lofHistory(:)], ...
                {'iter', 'LOF_pct'}, 'MCR_conv', 'Export convergence'), ...
                'Export conv');
            
            function safeCb(fcn, label)
                try
                    fcn();
                catch ME
                    uialert(f, sprintf('"%s" failed:\n%s\n\n(In %s line %d)', ...
                        label, ME.message, ...
                        ME.stack(1).name, ME.stack(1).line), ...
                        'Callback error');
                    lblStatus.Text = sprintf('Error in %s.', label);
                end
            end
            
            % ============== Callbacks ==============
            function doRun()
                N = round(edN.Value);
                opts = struct( ...
                    'nnC',     cbNNc.Value, ...
                    'uniC',    cbUni.Value, ...
                    'nnS',     cbNNs.Value, ...
                    'tol',     edTol.Value, ...
                    'maxIter', round(edMaxIt.Value), ...
                    'initMode', ddInit.Value, ...
                    'seed',    edSeed.Value, ...
                    'Cinit',   [], ...
                    'Sinit',   []);
                
                % "GA result": the spectra of the stored global fit as
                % starting spectra -- the EADS (SAS for a fit with
                % stretched components; sequential picture, tau
                % ascending) or the DADS (parallel picture, table
                % order), as many as there are, at most N. computeMCR
                % completes a shorter set from the data (a spectrum of
                % zeros would stay zero) and falls back to the SVD start
                % when nothing usable is there.
                initNote = '';      % why the SVD start was used instead
                gaKind   = '';      % which GA spectra seed the run
                gaNote   = '';
                if any(strcmp(opts.initMode, {'ga', 'gaeads'}))
                    nWL = numel(app.wavelength);
                    Sga = [];
                    if strcmp(opts.initMode, 'gaeads') && ...
                            ~isempty(app.gaResultEADS) && ...
                            size(app.gaResultEADS, 1) == nWL
                        Sga = app.gaResultEADS;
                        if any(app.gaResultStretchOn)
                            gaKind = 'SAS';
                        else
                            gaKind = 'EADS';
                        end
                    elseif ~isempty(app.gaResultDADS) && ...
                            size(app.gaResultDADS, 1) == nWL
                        Sga = app.gaResultDADS;
                        gaKind = 'DADS';
                        if strcmp(opts.initMode, 'gaeads')
                            gaNote = '      (EADS not available)';
                        end
                    end
                    if isempty(Sga)
                        opts.initMode = 'svd';
                        initNote = ' (no GA result)';
                    else
                        opts.Sinit    = Sga(:, 1:min(N, size(Sga, 2)));
                        opts.initMode = 'custom';
                    end
                end
                
                % Data: rows = delay, cols = wavelength.
                Dmat = app.deltaA.';
                if size(Dmat, 1) < N || size(Dmat, 2) < N
                    uialert(f, sprintf( ...
                        'Need at least %d delays and %d wavelengths.', ...
                        N, N), 'Too few points');
                    return;
                end
                
                lblStatus.Text = 'Running MCR-ALS...';  drawnow;
                tic;
                [C, S, info] = app.computeMCR(Dmat, N, opts);
                elapsed = toc;
                
                % Drec from computeMCR is (delay x wavelength); transpose
                % so it matches app.deltaA orientation (wavelength x delay)
                % which is what updateRec/exportRec expect.
                Drec = (C * S.').';
                mcrResult.C          = C;
                mcrResult.S          = S;
                mcrResult.Drec       = Drec;
                mcrResult.lof        = info.lof;
                mcrResult.iter       = info.iter;
                mcrResult.lofHistory = info.lofHistory;
                mcrResult.converged  = info.converged;
                
                btnExpC.Enable    = 'on';
                btnExpS.Enable    = 'on';
                btnExpRec.Enable  = 'on';
                btnExpConv.Enable = 'on';
                
                if info.converged, tag = 'converged'; else, tag = 'NOT converged'; end
                lblStatus.Text = sprintf( ...
                    'Done. LOF = %.3f%% (%s, %d iters, %.1fs)', ...
                    info.lof, tag, info.iter, elapsed);
                lines = {};
                lines{end+1} = sprintf('Components: %d', N);
                lines{end+1} = sprintf('LOF: %.4f%%', info.lof);
                lines{end+1} = sprintf('Iterations: %d  (%s)', ...
                    info.iter, tag);
                lines{end+1} = '';
                lines{end+1} = 'Constraints:';
                lines{end+1} = sprintf('  NN(C):  %s', tf2str(opts.nnC));
                lines{end+1} = sprintf('  Uni(C): %s', tf2str(opts.uniC));
                lines{end+1} = sprintf('  NN(S):  %s', tf2str(opts.nnS));
                lines{end+1} = '';
                switch opts.initMode
                    case 'random'
                        lines{end+1} = sprintf('Init: random (seed %d)', ...
                            info.seed);
                    case 'custom'
                        if isnan(info.nSinit)
                            lines{end+1} = 'Init: SVD (GA spectra unusable)';
                        else
                            lines{end+1} = sprintf( ...
                                'Init: GA result (%d %s)', info.nSinit, ...
                                gaKind);
                            if ~isempty(gaNote)
                                lines{end+1} = gaNote;
                            end
                            if info.nSinit < N
                                lines{end+1} = sprintf( ...
                                    '      + %d from the data (SVD)', ...
                                    N - info.nSinit);
                            end
                        end
                    otherwise
                        lines{end+1} = ['Init: SVD' initNote];
                end
                lines{end+1} = '';
                lines{end+1} = 'Click on either plot to see';
                lines{end+1} = 'reconstruction at that wavelength.';
                txtResult.Value = lines;
                
                plotResults();
                
                % Persist so the analysis survives closing this window.
                app.mcrSaved = struct( ...
                    'result',  mcrResult, ...
                    'curLambda', curLambda, ...
                    'status',  lblStatus.Text, ...
                    'lines',   {txtResult.Value}, ...
                    'settings', struct( ...
                        'N',    edN.Value, ...
                        'nnC',  cbNNc.Value, 'uniC', cbUni.Value, ...
                        'nnS',  cbNNs.Value, ...
                        'init', ddInit.Value, 'seed', edSeed.Value, ...
                        'tol',  edTol.Value, 'maxIt', edMaxIt.Value, ...
                        'tscale', ddTimeScale.Value));
            end
            
            function plotResults()
                if isempty(mcrResult.C), return; end
                t = app.delay(:);  wl = app.wavelength(:);
                N = size(mcrResult.C, 2);
                cmap = lines(N);
                
                app.claAll(axC);  hold(axC, 'on');
                for jj = 1:N
                    plot(axC, t, mcrResult.C(:, jj), '-', ...
                        'Color', cmap(jj, :), 'LineWidth', 1.4, ...
                        'DisplayName', sprintf('C_{%d}', jj));
                end
                yline(axC, 0, 'k:', 'HandleVisibility', 'off');
                hold(axC, 'off');
                xlabel(axC, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axC, 'Concentration');
                title(axC, 'C profiles');
                if strcmp(ddTimeScale.Value, 'log') && any(t > 0)
                    set(axC, 'XScale', 'log');
                    xlim(axC, [min(t(t > 0)), max(t)]);
                else
                    set(axC, 'XScale', 'linear');
                    xlim(axC, [min(t), max(t)]);
                end
                legend(axC, 'Location', 'best', 'FontSize', 8);
                grid(axC, 'on');
                
                app.claAll(axS);  hold(axS, 'on');
                for jj = 1:N
                    plot(axS, wl, mcrResult.S(:, jj), '-', ...
                        'Color', cmap(jj, :), 'LineWidth', 1.4, ...
                        'DisplayName', sprintf('S_{%d}', jj));
                end
                yline(axS, 0, 'k:', 'HandleVisibility', 'off');
                if isfinite(curLambda)
                    xhS = xline(axS, curLambda, 'k-', 'LineWidth', 1.0, ...
                        'Alpha', 0.7, 'HandleVisibility', 'off', ...
                        'HitTest', 'off', 'PickableParts', 'none');
                else
                    xhS = [];
                end
                hold(axS, 'off');
                xlabel(axS, app.xAxLabel());
                ylabel(axS, 'Pure spectrum (\DeltaA or Intensity)');
                title(axS, 'S spectra');
                xlim(axS, [min(wl), max(wl)]);
                legend(axS, 'Location', 'best', 'FontSize', 8);
                grid(axS, 'on');
                
                % Convergence
                cla(axConv);
                if numel(mcrResult.lofHistory) > 1
                    plot(axConv, 1:numel(mcrResult.lofHistory), ...
                        mcrResult.lofHistory, 'b-o', ...
                        'MarkerSize', 3, 'LineWidth', 1.0);
                end
                xlabel(axConv, 'ALS iteration');
                ylabel(axConv, 'LOF (%)');
                title(axConv, sprintf('Convergence (final LOF = %.3f%%)', ...
                    mcrResult.lof));
                set(axConv, 'YScale', 'log');
                grid(axConv, 'on');
                
                % Click handlers on C and S to set the inspection wavelength.
                axS.ButtonDownFcn = @(~,ev) onPlotClick(ev, 'S');
                ch = axS.Children;
                for jj = 1:numel(ch)
                    if isa(ch(jj), 'matlab.graphics.chart.primitive.Line')
                        ch(jj).ButtonDownFcn = @(~,ev) onPlotClick(ev, 'S');
                    end
                end
                
                updateRec();
            end
            
            function onPlotClick(ev, ~)
                pt = ev.IntersectionPoint;
                if isempty(pt), return; end
                draggingS = true;
                applyMcrPoint(pt);
            end
            function onMcrMotion()
                if ~draggingS, return; end
                if ~app.mouseInsideFigure(f)
                    onMcrUp();  return;
                end
                cp = axS.CurrentPoint;
                if isempty(cp) || size(cp,1) < 1, return; end
                applyMcrPoint(cp(1, 1:2));
            end
            function onMcrUp()
                draggingS = false;
            end
            function applyMcrPoint(pt)
                if numel(pt) < 1, return; end
                wlR = [min(app.wavelength), max(app.wavelength)];
                curLambda = min(max(pt(1), wlR(1)), wlR(2));
                if ~isempty(xhS) && isvalid(xhS)
                    xhS.Value = curLambda;
                end
                updateRec();
                drawnow limitrate;
            end
            
            function updateRec()
                if isempty(mcrResult.Drec), return; end
                [~, iw] = min(abs(app.wavelength - curLambda));
                wlAct = app.wavelength(iw);
                
                app.claAll(axRec);  hold(axRec, 'on');
                plot(axRec, app.delay, app.deltaA(iw, :), 'k.', ...
                    'MarkerSize', 6, 'DisplayName', 'Data');
                plot(axRec, app.delay, mcrResult.Drec(iw, :), 'r-', ...
                    'LineWidth', 1.4, 'DisplayName', 'MCR reconstr.');
                plot(axRec, app.delay, ...
                    app.deltaA(iw, :) - mcrResult.Drec(iw, :), '-', ...
                    'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
                    'DisplayName', 'Residual');
                yline(axRec, 0, 'k:', 'HandleVisibility', 'off');
                hold(axRec, 'off');
                xlabel(axRec, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axRec, '\DeltaA or Intensity');
                title(axRec, sprintf( ...
                    'Reconstruction at %s = %.2f %s', ...
                    app.xSym(), wlAct, app.xUnitAx()));
                if strcmp(ddTimeScale.Value, 'log') && any(app.delay > 0)
                    set(axRec, 'XScale', 'log');
                    xlim(axRec, [min(app.delay(app.delay > 0)), ...
                        max(app.delay)]);
                else
                    set(axRec, 'XScale', 'linear');
                    xlim(axRec, [min(app.delay), max(app.delay)]);
                end
                legend(axRec, 'Location', 'best', 'FontSize', 8);
                grid(axRec, 'on');
            end
            
            function onTimeScaleChange()
                % Apply the chosen time scale to both panels that show
                % delay on x. Only redraw when we already have results.
                if isempty(mcrResult.C), return; end
                plotResults();
                if ~isempty(mcrResult.Drec), updateRec(); end
            end
            
            function onInitChange()
                % The seed belongs to the random start only.
                if strcmp(ddInit.Value, 'random')
                    edSeed.Enable         = 'on';
                    lblSeed.FontColor     = [0.10 0.10 0.10];
                    lblSeedHint.FontColor = [0.45 0.45 0.45];
                else
                    edSeed.Enable         = 'off';
                    lblSeed.FontColor     = [0.60 0.60 0.60];
                    lblSeedHint.FontColor = [0.70 0.70 0.70];
                end
            end
            
            function doReset()
                edN.Value      = 3;
                cbNNc.Value    = true;
                cbUni.Value    = false;
                cbNNs.Value    = false;
                ddInit.Value   = 'svd';
                edSeed.Value   = 0;
                onInitChange();
                edTol.Value    = 1e-4;
                edMaxIt.Value  = 200;
                ddTimeScale.Value = 'log';
                mcrResult = struct('C', [], 'S', [], 'Drec', [], ...
                    'lof', NaN, 'iter', 0, 'lofHistory', [], ...
                    'converged', false);
                xhS = [];
                btnExpC.Enable    = 'off';
                btnExpS.Enable    = 'off';
                btnExpRec.Enable  = 'off';
                btnExpConv.Enable = 'off';
                app.claAll(axC);  app.claAll(axS);  app.claAll(axRec);  cla(axConv);
                title(axC, 'C profiles (concentration vs time)');
                title(axS, 'S spectra');
                title(axRec, 'Reconstruction at clicked \lambda');
                title(axConv, 'Convergence');
                lblStatus.Text  = 'Reset.';
                txtResult.Value = {'Click "Run MCR-ALS" to compute.'};
                app.mcrSaved = [];
            end
            
            function exportRec()
                if isempty(mcrResult.Drec), return; end
                t = app.delay(:).';
                hdr = [{app.xHdr()}, ...
                    arrayfun(@(x) sprintf('t=%.4g_%s', x, app.tUnitHdr()), ...
                        t, 'UniformOutput', false)];
                mat = [app.wavelength(:), mcrResult.Drec];
                exportMat(mat, hdr, 'MCR_recon', 'Export reconstruction');
            end
            
            function exportMat(mat, hdr, defName, dlgTitle)
                [file, path] = app.putFile( ...
                    {'*.csv','CSV (*.csv)'; '*.dat','DAT (tab, *.dat)'}, ...
                    dlgTitle, [defName '.csv']);
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    T = array2table(mat, 'VariableNames', hdr);
                    writetable(T, full, 'Delimiter', dlm);
                    uialert(f, sprintf('Saved: %s', full), ...
                        'Exported', 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Error');
                end
            end
            
            function s = tf2str(b)
                if b, s = 'ON'; else, s = 'off'; end
            end
            
            function onClose()
                app.mcrFig = [];
                delete(f);
            end
        end
        
        function openCoherenceWindow(app)
            %OPENCOHERENCEWINDOW  Vibrational coherence (oscillation map).
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');
                return;
            end
            if ~isempty(app.coherenceFig) && isvalid(app.coherenceFig)
                figure(app.coherenceFig);  return;
            end
            
            f = uifigure('Name', 'Vibrational Coherence Analysis', ...
                'Position', app.fitFigurePosition(1380, 900));
            app.applyLightTheme(f);
            app.coherenceFig = f;
            f.CloseRequestFcn = @(~,~) onClose();
            f.WindowButtonMotionFcn = @(~,~) onCohMotion();
            f.WindowButtonUpFcn     = @(~,~) onCohUp();
            
            gl = uigridlayout(f, [1, 2]);
            gl.ColumnWidth = {360, '1x'};
            gl.Padding = [10 10 10 10];  gl.ColumnSpacing = 10;
            
            % ============== Left: setup ==============
            left = uipanel(gl, 'Title', 'Coherence setup');
            left.Layout.Column = 1;
            lgl = uigridlayout(left, [26, 1]);
            % rows: residualLbl, residRow, NexpRow, twindowLbl, tMinRow,
            %       tMaxRow, fftLbl, apodRow, zpadRow, detrendRow,
            %       freqLbl, freqUnitRow, normRow, yScaleRow,
            %       avgLbl, avgMinRow, avgMaxRow,
            %       dispLbl, dispMinRow, dispMaxRow,
            %       btnRow, tfBtnRow, expLbl, expRow, statusLbl, resultArea
            lgl.RowHeight = {22, 30, 30, 22, 30, 30, 22, 30, 30, 30, ...
                             22, 30, 30, 30, 22, 30, 30, ...
                             22, 30, 30, ...
                             32, 30, 22, 30, 22, '1x'};
            lgl.Padding = [8 8 8 8];  lgl.RowSpacing = 6;
            
            % --- Residual source ---
            uilabel(lgl, 'Text', 'Residual source:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Method:');
            ddResid = uidropdown(row, ...
                'Items', {'Use GA residual', 'Per-pixel exp fit'}, ...
                'ItemsData', {'ga', 'pixfit'}, ...
                'Value', 'ga');
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            lblNexp = uilabel(row, 'Text', 'Per-pixel N exp:');
            edNexp = uieditfield(row, 'numeric', 'Value', 2, ...
                'Limits', [1 5], 'RoundFractionalValues', 'on', ...
                'Tooltip', ['Number of exponentials fitted at every ' ...
                    'wavelength. Used by "Per-pixel exp fit" only.']);
            
            % --- Time window ---
            uilabel(lgl, 'Text', 'Time window:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('t_min (%s):', app.tUnitTxt()));
            edTMin = uieditfield(row, 'numeric', ...
                'Value', defaultTMin());
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('t_max (%s):', app.tUnitTxt()));
            edTMax = uieditfield(row, 'numeric', ...
                'Value', defaultTMax());
            
            % --- FFT options ---
            uilabel(lgl, 'Text', 'FFT options:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Apodization:');
            ddApod = uidropdown(row, ...
                'Items', {'Rectangular', 'Hann', 'Hamming', 'Blackman'}, ...
                'ItemsData', {'rect', 'hann', 'hamming', 'blackman'}, ...
                'Value', 'hann');
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Zero-pad factor:');
            edZpad = uieditfield(row, 'numeric', 'Value', 2, ...
                'Limits', [1 16], 'RoundFractionalValues', 'on');
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Detrend:');
            ddDet = uidropdown(row, ...
                'Items', {'None', 'Mean', 'Linear'}, ...
                'ItemsData', {'none', 'mean', 'linear'}, ...
                'Value', 'linear');
            
            % --- Output ---
            uilabel(lgl, 'Text', 'Output:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Frequency unit:');
            ddFreqU = uidropdown(row, ...
                'Items', app.uiStr({'Wavenumber (cm^{-1})', 'THz'}), ...
                'ItemsData', {'cm-1', 'thz'}, ...
                'Value', 'cm-1');
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Normalize:');
            ddNorm = uidropdown(row, ...
                'Items', {'Peak (global)', 'Per-wavelength', 'None'}, ...
                'ItemsData', {'peak', 'perwl', 'none'}, ...
                'Value', 'peak');
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Heatmap z / 1D y scale:');
            ddYScale = uidropdown(row, ...
                'Items', {'Linear', 'Log'}, ...
                'ItemsData', {'linear', 'log'}, ...
                'Value', 'linear', ...
                'Tooltip', app.uiStr(['Apply log scale to both top plots: ' ...
                    'the |FFT|^2 colormap on the heatmap (z-axis) ' ...
                    'and the y-axis of the power spectrum at \lambda. ' ...
                    'Log highlights weak modes.']), ...
                'ValueChangedFcn', @(~,~) onYScaleChange());
            
            % --- Averaging window ---
            uilabel(lgl, 'Text', app.uiStr('Average power over \lambda range:'), ...
                'FontWeight', 'bold');
            wlLo0 = min(app.wavelength);  wlHi0 = max(app.wavelength);
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('%s_min (%s):', app.xSymTxt(), app.xUnitTxt()));
            edAvgWLmin = uieditfield(row, 'numeric', 'Value', wlLo0, ...
                'Limits', [wlLo0, wlHi0], ...
                'ValueChangedFcn', @(~,~) onAvgRangeChange());
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('%s_max (%s):', app.xSymTxt(), app.xUnitTxt()));
            edAvgWLmax = uieditfield(row, 'numeric', 'Value', wlHi0, ...
                'Limits', [wlLo0, wlHi0], ...
                'ValueChangedFcn', @(~,~) onAvgRangeChange());
            
            % --- Display nu range (visual zoom on the three frequency
            %     plots: heatmap, P(nu) at lambda, and the averaged
            %     P(nu). Does NOT recompute the FFT -- only the x-axis
            %     ranges are clipped.)
            uilabel(lgl, 'Text', app.uiStr('Display \nu range:'), ...
                'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', app.uiStr('\nu_min:'));
            edDispVmin = uieditfield(row, 'numeric', 'Value', 0, ...
                'ValueChangedFcn', @(~,~) onDispRangeChange());
            
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', app.uiStr('\nu_max:'));
            edDispVmax = uieditfield(row, 'numeric', 'Value', 0, ...
                'Tooltip', app.uiStr(['Visualization window on the wavenumber ' ...
                    'axis. \nu_max = 0 means "auto" (full Nyquist ' ...
                    'range of the current FFT).']), ...
                'ValueChangedFcn', @(~,~) onDispRangeChange());
            
            % --- Buttons ---
            btnRow = uigridlayout(lgl, [1, 2]);
            btnRow.ColumnWidth = {'1x', '1x'};
            btnRow.Padding = [0 0 0 0];  btnRow.ColumnSpacing = 6;
            btnRun = uibutton(btnRow, 'Text', 'Run Coherence', ...
                'BackgroundColor', [0.30 0.65 0.75], ...
                'FontColor', 'w', 'FontWeight', 'bold');
            btnReset = uibutton(btnRow, 'Text', 'Reset');
            
            btnTF = uibutton(lgl, 'Text', 'Time-frequency (STFT / wavelet)...', ...
                'Enable', 'off', ...
                'Tooltip', ['Open the nu-t analysis of this residual: ' ...
                    'when does each mode live and how fast does it ' ...
                    'dephase. Available after a run.'], ...
                'ButtonPushedFcn', @(~,~) app.openTimeFrequencyWindow());
            
            % --- Export ---
            uilabel(lgl, 'Text', 'Export:', 'FontWeight', 'bold');
            expRow = uigridlayout(lgl, [1, 4]);
            expRow.ColumnWidth = {'1x', '1x', '1x', '1x'};
            expRow.Padding = [0 0 0 0];  expRow.ColumnSpacing = 6;
            btnExpMap   = uibutton(expRow, 'Text', app.uiStr('|FFT|^2 map'), ...
                'Enable', 'off');
            btnExpSpec  = uibutton(expRow, 'Text', app.uiStr('P(\nu)@\lambda'), ...
                'Enable', 'off');
            btnExpAvg   = uibutton(expRow, 'Text', app.uiStr('Avg P(\nu)'), ...
                'Enable', 'off', ...
                'Tooltip', app.uiStr(['Export the power spectrum averaged ' ...
                    'over the current "Average over \lambda range" ' ...
                    'window.']));
            btnExpProf  = uibutton(expRow, 'Text', app.uiStr('A(\lambda)@\nu'), ...
                'Enable', 'off');
            
            lblStatus = uilabel(lgl, 'Text', 'Ready.', ...
                'FontColor', [0.2 0.2 0.2]);
            txtResult = uitextarea(lgl, 'Editable', 'off', ...
                'Value', {'Click "Run Coherence" to compute the FFT map.'}, ...
                'FontName', 'Monospaced');
            
            % ============== Right: plots ==============
            right = uipanel(gl, 'Title', 'Results');
            right.Layout.Column = 2;
            rgl = uigridlayout(right, [2, 2]);
            rgl.ColumnWidth = {'1x', '1x'};
            rgl.RowHeight   = {'1.3x', '1x'};
            rgl.Padding     = [8 8 8 8];  rgl.RowSpacing = 8;
            rgl.ColumnSpacing = 10;
            
            % Heatmap (large) + 1D power spectrum at lambda
            axMap = uiaxes(rgl);
            axMap.Layout.Row = 1;  axMap.Layout.Column = 1;
            xlabel(axMap, freqAxLabel('cm-1'));
            ylabel(axMap, app.xAxLabel());
            title(axMap, '|FFT|^2 map');
            % UIAxes default interactions (datatip / rotate3d / zoom /
            % pan) can intercept mouse-down/up sequences on a 3D surface
            % such as the heatmap, which prevents our figure-level
            % WindowButtonUpFcn from firing and leaves dragging stuck.
            % Strip them so only our ButtonDownFcn drives the drag.
            try, disableDefaultInteractivity(axMap); catch, end
            try, axMap.Toolbar.Visible = 'off'; catch, end
            try, axMap.Interactions = []; catch, end
            
            % Power spectrum panel: split top (at clicked lambda) and
            % bottom (averaged over a user-set wavelength range).
            pspecPnl = uigridlayout(rgl, [2, 1]);
            pspecPnl.Layout.Row = 1;  pspecPnl.Layout.Column = 2;
            pspecPnl.RowHeight = {'1x', '1x'};
            pspecPnl.Padding = [0 0 0 0];  pspecPnl.RowSpacing = 4;
            
            axPSpec = uiaxes(pspecPnl);
            axPSpec.Layout.Row = 1;
            xlabel(axPSpec, '');                  % shared x with avg below
            ylabel(axPSpec, '|FFT|^2');
            title(axPSpec, 'Power spectrum at clicked \lambda');
            grid(axPSpec, 'on');
            
            axPSpecAvg = uiaxes(pspecPnl);
            axPSpecAvg.Layout.Row = 2;
            xlabel(axPSpecAvg, freqAxLabel('cm-1'));
            ylabel(axPSpecAvg, '|FFT|^2');
            title(axPSpecAvg, 'Average over \lambda range');
            grid(axPSpecAvg, 'on');
            
            % Profile at fixed wavenumber: split top (|FFT|^2) and
            % bottom (phase). The sub-grid takes the (Row 2, Col 1) cell
            % of the parent rgl that originally held a single axProf.
            profPnl = uigridlayout(rgl, [2, 1]);
            profPnl.Layout.Row = 2;  profPnl.Layout.Column = 1;
            profPnl.RowHeight = {'1x', '1x'};
            profPnl.Padding = [0 0 0 0];  profPnl.RowSpacing = 4;
            
            axProf = uiaxes(profPnl);
            axProf.Layout.Row = 1;
            xlabel(axProf, '');                  % shared x with phase below
            ylabel(axProf, '|FFT|^2');
            title(axProf, 'Spectral profile at clicked \nu');
            grid(axProf, 'on');
            
            axProfPh = uiaxes(profPnl);
            axProfPh.Layout.Row = 2;
            xlabel(axProfPh, app.xAxLabel());
            ylabel(axProfPh, 'Phase (rad)');
            grid(axProfPh, 'on');
            
            axRes = uiaxes(rgl);
            axRes.Layout.Row = 2;  axRes.Layout.Column = 2;
            xlabel(axRes, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axRes, 'Residual');
            title(axRes, 'Residual time trace at clicked \lambda');
            grid(axRes, 'on');
            
            % ---- State ----
            cohResult = struct('Pmap', [], 'freq', [], ...
                'Rresid', [], 'tWin', [], ...
                'info', struct());
            curLambda = app.selWL;
            curFreq   = NaN;
            xhV       = [];   % vertical crosshair on heatmap (freq)
            xhH       = [];   % horizontal crosshair on heatmap (lambda)
            yhAvgLo   = [];   % heatmap horizontal line: avg lambda_min
            yhAvgHi   = [];   % heatmap horizontal line: avg lambda_max
            draggingC = false; % coherence heatmap drag state
            cohDragTimer = []; % polling timer used as drag safety net
            
            % ---- Wire up ----
            btnRun.ButtonPushedFcn   = @(~,~) safeCb(@doRun, 'Run Coherence');
            btnReset.ButtonPushedFcn = @(~,~) safeCb(@doReset, 'Reset');
            ddResid.ValueChangedFcn  = @(~,~) safeCb(@syncNexpEnable, 'Method');
            btnExpMap.ButtonPushedFcn  = @(~,~) safeCb(@exportMap, 'Export map');
            btnExpSpec.ButtonPushedFcn = @(~,~) safeCb(@exportSpec, 'Export spec');
            btnExpAvg.ButtonPushedFcn  = @(~,~) safeCb(@exportAvg,  'Export avg spec');
            btnExpProf.ButtonPushedFcn = @(~,~) safeCb(@exportProf, 'Export profile');
            
            % ---- Restore previous analysis (if any) ----
            % A completed run is snapshotted into app.cohSaved, so a
            % closed / reopened window comes back exactly where the
            % user left it. Restore only when the saved map is still
            % dimensionally consistent with the current wavelength
            % axis (crop or re-load invalidates it).
            if ~isempty(app.cohSaved) && ...
                    isstruct(app.cohSaved) && ...
                    isfield(app.cohSaved, 'result') && ...
                    ~isempty(app.cohSaved.result.Pmap) && ...
                    size(app.cohSaved.result.Pmap, 1) == numel(app.wavelength)
                try
                    st = app.cohSaved.settings;
                    ddResid.Value  = st.resid;   edNexp.Value = st.nexp;
                    edTMin.Value   = st.tmin;    edTMax.Value = st.tmax;
                    ddApod.Value   = st.apod;    edZpad.Value = st.zpad;
                    ddDet.Value    = st.detrend; ddFreqU.Value = st.frequ;
                    ddNorm.Value   = st.norm;    ddYScale.Value = st.yscale;
                    edAvgWLmin.Value = st.avgLo; edAvgWLmax.Value = st.avgHi;
                    edDispVmin.Value = st.dispLo; edDispVmax.Value = st.dispHi;
                    cohResult = app.cohSaved.result;
                    curLambda = app.cohSaved.curLambda;
                    curFreq   = app.cohSaved.curFreq;
                    btnExpMap.Enable  = 'on';
                    btnExpSpec.Enable = 'on';
                    btnExpAvg.Enable  = 'on';
                    btnExpProf.Enable = 'on';
                    btnTF.Enable      = 'on';
                    plotResults();
                    lblStatus.Text  = sprintf('%s  [restored]', ...
                        app.cohSaved.status);
                    txtResult.Value = app.cohSaved.lines;
                catch
                    % Restoration is best-effort; say so instead of
                    % leaving a half-restored window unexplained.
                    lblStatus.Text = 'Stored result could not be restored -- run again.';
                end
            end
            syncNexpEnable();
            
            function syncNexpEnable()
                % "Per-pixel N exp" belongs to the per-pixel fit only;
                % with the GA residual it has no effect, so grey it out.
                if strcmp(ddResid.Value, 'pixfit')
                    onOff = 'on';
                else
                    onOff = 'off';
                end
                edNexp.Enable  = onOff;
                lblNexp.Enable = onOff;
            end
            
            function onYScaleChange()
                % Apply log scale to the two top plots:
                %   * heatmap : color z-axis (CData = log10(Pmap))
                %   * power spectrum at lambda : YScale
                % The heatmap's CData has to be recomputed, so we
                % call plotResults() (which then calls update1Ds).
                if ~isempty(cohResult.Pmap)
                    plotResults();
                end
            end
            
            function onAvgRangeChange()
                % Refresh just the heatmap range-marker lines and the
                % averaged-power plot (no full plotResults() needed).
                if isempty(cohResult.Pmap), return; end
                if ~isempty(yhAvgLo) && isvalid(yhAvgLo)
                    yhAvgLo.Value = edAvgWLmin.Value;
                end
                if ~isempty(yhAvgHi) && isvalid(yhAvgHi)
                    yhAvgHi.Value = edAvgWLmax.Value;
                end
                update1Ds();
            end
            
            function onDispRangeChange()
                % Apply the new visualization window to the three
                % wavenumber plots without recomputing the FFT.
                % Additionally renormalize the heatmap's color scale to
                % the data inside the visible wavenumber band, so that
                % zooming into a quieter region brings its small modes
                % into the full color range automatically.
                if isempty(cohResult.Pmap), return; end
                lims = getDispLims();
                if isempty(lims), return; end
                try, xlim(axMap,      lims); catch, end
                try, xlim(axPSpec,    lims); catch, end
                try, xlim(axPSpecAvg, lims); catch, end
                
                freq = cohResult.freq;
                fMask = freq >= lims(1) & freq <= lims(2);
                if any(fMask)
                    sub = cohResult.Pmap(:, fMask);
                    if strcmp(ddYScale.Value, 'log')
                        pkSub = max(sub(:), [], 'omitnan');
                        if ~isfinite(pkSub) || pkSub <= 0, pkSub = 1; end
                        floorVal = pkSub * 1e-4;   % 4-decade range
                        try
                            axMap.CLim = [log10(floorVal), log10(pkSub)];
                        catch
                        end
                    else
                        cl = max(sub(:), [], 'omitnan');
                        if ~isfinite(cl) || cl <= 0, cl = 1; end
                        try, axMap.CLim = [0, cl]; catch, end
                    end
                end
            end
            
            function lims = getDispLims()
                % Resolve the user-set display window to a valid
                % [lo hi] pair. nu_max <= 0 (the default) means
                % "auto" -- use the full Nyquist range.
                lims = [];
                if isempty(cohResult.freq), return; end
                fLo = min(cohResult.freq);
                fHi = max(cohResult.freq);
                lo = edDispVmin.Value;
                hi = edDispVmax.Value;
                if ~isfinite(hi) || hi <= 0, hi = fHi; end
                if ~isfinite(lo), lo = fLo;            end
                if lo > hi, tmp = lo; lo = hi; hi = tmp; end
                lo = max(lo, fLo);
                hi = min(hi, fHi);
                if lo >= hi, return; end
                lims = [lo, hi];
            end
            
            function lims = getXlimForFreq(freqAxis)
                % Convenience wrapper used inside plotResults /
                % update1Ds where the local "freq" vector is in
                % scope. Falls back to the full range when no display
                % window is set.
                lims = getDispLims();
                if isempty(lims)
                    lims = [min(freqAxis), max(freqAxis)];
                end
            end
            
            function safeCb(fcn, label)
                try
                    fcn();
                catch ME
                    uialert(f, sprintf('"%s" failed:\n%s\n\n(In %s line %d)', ...
                        label, ME.message, ...
                        ME.stack(1).name, ME.stack(1).line), ...
                        'Callback error');
                    lblStatus.Text = sprintf('Error in %s.', label);
                end
            end
            
            % ============== Helpers ==============
            function v = defaultTMin()
                % Reasonable starting time after IRF: 3*sigma, in
                % the active time unit. If GA FWHM is set, use that.
                if ~isempty(app.gaResultFWHM) && isfinite(app.gaResultFWHM)
                    sigma = app.gaResultFWHM / (2 * sqrt(2 * log(2)));
                elseif ~isempty(app.gaFWHM) && isfinite(app.gaFWHM)
                    sigma = app.gaFWHM / (2 * sqrt(2 * log(2)));
                else
                    sigma = 0;
                end
                switch app.timeUnit
                    case 'fs', v = max(100, 3 * sigma);
                    case 'ps', v = max(0.3, 3 * sigma);
                    case 'ns', v = max(3e-4, 3 * sigma);
                    case 'us', v = max(3e-7, 3 * sigma);
                    otherwise, v = max(0.3, 3 * sigma);
                end
            end
            function v = defaultTMax()
                tspan = max(app.delay) - min(app.delay);
                if tspan <= 0, tspan = 1000; end
                switch app.timeUnit
                    case 'fs', v = min(5000, max(app.delay));
                    case 'ps', v = min(5, max(app.delay));
                    case 'ns', v = min(0.005, max(app.delay));
                    case 'us', v = min(5e-6, max(app.delay));
                    otherwise, v = min(5, max(app.delay));
                end
            end
            function lab = freqAxLabel(unit)
                switch lower(unit)
                    case 'cm-1', lab = 'Wavenumber (cm^{-1})';
                    case 'thz',  lab = 'Frequency (THz)';
                    otherwise,   lab = 'Frequency (Hz)';
                end
            end
            
            function [R, note] = computeResidual(method, Nexp)
                % Returns R (Nwl x Nt), same orientation as app.deltaA.
                % note: non-empty when some rows could not be fitted.
                note = '';
                switch method
                    case 'ga'
                        if isempty(app.gaResultFit)
                            error(['No GA fit available. Run Global ' ...
                                'Analysis first, or switch to ' ...
                                '"Per-pixel exp fit".']);
                        end
                        R = app.deltaA - app.gaResultFit;
                    case 'pixfit'
                        % Per-row multi-exp fit. Use the wavelength-
                        % independent IRF parameters from GA setup.
                        [R, nFail, nTried] = perPixelResidual(Nexp);
                        if nFail > 0
                            note = sprintf(['%d of %d rows: exp fit ' ...
                                'failed (mean-subtracted)'], nFail, nTried);
                        end
                    otherwise
                        error('Unknown residual method: %s', method);
                end
            end
            
            function [R, nFail, nTried] = perPixelResidual(Nexp)
                % Fit each row (wavelength) with N independent
                % exponentials sharing a Gaussian IRF (t0, FWHM from
                % GA setup).  Slow but model-light.
                Nwl = size(app.deltaA, 1);
                Nt  = size(app.deltaA, 2);
                R   = zeros(Nwl, Nt);
                t   = app.delay;
                if ~isempty(app.gaResultT0),   t0in = app.gaResultT0; else, t0in = app.gaT0;   end
                if ~isempty(app.gaResultFWHM), fwin = app.gaResultFWHM; else, fwin = app.gaFWHM; end
                tauInit  = logspace(log10(max(0.1, fwin)), ...
                                    log10(max(app.delay)/2), Nexp);
                tauFixed = false(1, Nexp);
                betaInit = ones(1, Nexp);
                betaFixed = true(1, Nexp);
                stretchOn = false(1, Nexp);   % plain exponentials
                lblStatus.Text = sprintf( ...
                    'Per-pixel exp fit on %d rows...', Nwl);
                drawnow;
                nFail = 0;  nTried = 0;  lastMsg = '';
                for ii = 1:Nwl
                    y = app.deltaA(ii, :);
                    if all(~isfinite(y))
                        R(ii, :) = 0;  continue;
                    end
                    nTried = nTried + 1;
                    if mod(nTried, 20) == 0
                        % One fit per wavelength takes a while: show
                        % progress (without running other callbacks).
                        lblStatus.Text = sprintf( ...
                            'Per-pixel exp fit: row %d of %d...', ii, Nwl);
                        drawnow limitrate nocallbacks;
                    end
                    try
                        [~, ~, ~, ~, ~, ~, resVec, ~] = ...
                            app.fitSingleTrace(t, y, tauInit, tauFixed, ...
                            betaInit, betaFixed, stretchOn, ...
                            t0in, true, fwin, true, false, 'numerical');
                        R(ii, :) = resVec;
                    catch ME
                        % Fall back: residual = data - mean
                        R(ii, :) = y - mean(y, 'omitnan');
                        nFail = nFail + 1;  lastMsg = ME.message;
                    end
                end
                % A fallback on every row means the fit itself is
                % broken for this data set -- report it instead of
                % silently returning mean-subtracted data.
                if nTried > 0 && nFail == nTried
                    error(['Per-pixel exp fit failed on every row ' ...
                        '(%s). Use the GA residual instead, or change ' ...
                        'the number of exponentials.'], lastMsg);
                end
            end
            
            % ============== Callbacks ==============
            function doRun()
                if edTMin.Value >= edTMax.Value
                    uialert(f, 'Need t_min < t_max.', 'Invalid');
                    return;
                end
                method = ddResid.Value;
                lblStatus.Text = 'Computing residual...';  drawnow;
                [Rresid, residNote] = computeResidual(method, ...
                    round(edNexp.Value));
                
                opts = struct( ...
                    'tMin', edTMin.Value, ...
                    'tMax', edTMax.Value, ...
                    'apod', ddApod.Value, ...
                    'zeroPad', round(edZpad.Value), ...
                    'detrend', ddDet.Value, ...
                    'freqUnit', ddFreqU.Value, ...
                    'normMode', ddNorm.Value);
                
                lblStatus.Text = 'Computing FFT map...';  drawnow;
                tic;
                [Pmap, freq, info] = app.computeCoherence( ...
                    Rresid, app.delay, opts);
                elapsed = toc;
                
                % Time-window mask for the residual that gets shown in
                % the time-domain panel
                mask = (app.delay >= edTMin.Value) & ...
                       (app.delay <= edTMax.Value);
                
                cohResult.Pmap   = Pmap;
                cohResult.freq   = freq;
                cohResult.Rresid = Rresid;
                cohResult.tWin   = app.delay(mask);
                cohResult.RinWin = Rresid(:, mask);
                cohResult.info   = info;
                
                btnExpMap.Enable  = 'on';
                btnExpSpec.Enable = 'on';
                btnExpAvg.Enable  = 'on';
                btnExpProf.Enable = 'on';
                btnTF.Enable      = 'on';
                % The time-frequency result was computed from the
                % previous residual -- drop it (and close its window).
                app.tfSaved = [];
                if ~isempty(app.tfFig) && isvalid(app.tfFig)
                    try, close(app.tfFig); catch, end
                end
                
                lblStatus.Text = sprintf( ...
                    'Done (%d freqs, %d wl, %.2fs).', ...
                    numel(freq), size(Pmap, 1), elapsed);
                lines = {};
                lines{end+1} = sprintf('Residual: %s', method);
                if ~isempty(residNote)
                    lines{end+1} = residNote;
                end
                lines{end+1} = sprintf('Window: %.4g .. %.4g %s', ...
                    edTMin.Value, edTMax.Value, app.tUnitTxt());
                lines{end+1} = sprintf('dt (uniform): %.4g %s', ...
                    info.dtUni, app.tUnitTxt());
                lines{end+1} = sprintf('Nfft: %d  (input %d)', ...
                    info.Nfft, info.Nwin);
                lines{end+1} = sprintf('Apod: %s   Zero-pad: %dx', ...
                    info.apod, info.zeroPad);
                lines{end+1} = sprintf('Detrend: %s', info.detrend);
                lines{end+1} = sprintf('Freq range: %.3g .. %.3g %s', ...
                    freq(1), freq(end), info.freqUnit);
                lines{end+1} = '';
                lines{end+1} = 'Click on the heatmap to inspect';
                lines{end+1} = 'a particular (lambda, nu).';
                txtResult.Value = lines;
                
                plotResults();
                
                % Persist to app so the analysis survives closing this
                % window. Settings snapshot lets reopen restore the UI.
                app.cohSaved = struct( ...
                    'result',   cohResult, ...
                    'curLambda', curLambda, 'curFreq', curFreq, ...
                    'status',   lblStatus.Text, ...
                    'lines',    {txtResult.Value}, ...
                    'settings', struct( ...
                        'resid',   ddResid.Value, 'nexp', edNexp.Value, ...
                        'tmin',    edTMin.Value,  'tmax', edTMax.Value, ...
                        'apod',    ddApod.Value,  'zpad', edZpad.Value, ...
                        'detrend', ddDet.Value,   'frequ', ddFreqU.Value, ...
                        'norm',    ddNorm.Value,  'yscale', ddYScale.Value, ...
                        'avgLo',   edAvgWLmin.Value, ...
                        'avgHi',   edAvgWLmax.Value, ...
                        'dispLo',  edDispVmin.Value, ...
                        'dispHi',  edDispVmax.Value));
            end
            
            function plotResults()
                if isempty(cohResult.Pmap), return; end
                wl   = app.wavelength(:).';
                freq = cohResult.freq(:).';
                Pmap = cohResult.Pmap;
                
                cmap = turboMap(256);
                
                % Choose the data subset used for the color-axis
                % normalization: if the user has set a Display nu range
                % we use only the visible band, otherwise the full FFT.
                dispLimsTmp = getDispLims();
                if ~isempty(dispLimsTmp)
                    fMaskDisp = freq >= dispLimsTmp(1) & ...
                                freq <= dispLimsTmp(2);
                    Psub = Pmap(:, fMaskDisp);
                    if isempty(Psub), Psub = Pmap; end
                else
                    Psub = Pmap;
                end
                
                % Build the surface CData. In 'log' mode we display
                % log10(P) so the colormap lays a logarithmic z-scale on
                % the heatmap; very small values are clipped to a floor
                % to avoid -Inf.
                if strcmp(ddYScale.Value, 'log')
                    pkAll = max(Psub(:), [], 'omitnan');
                    if ~isfinite(pkAll) || pkAll <= 0, pkAll = 1; end
                    floorVal = pkAll * 1e-4;       % 4 decades dyn. range
                    Pdisp    = max(Pmap, floorVal);
                    Cdat     = log10(Pdisp);
                    cLo = log10(floorVal);
                    cHi = log10(pkAll);
                    cbLabel = 'log_{10}|FFT|^2';
                else
                    Cdat   = Pmap;
                    cbLabel = '|FFT|^2';
                    cl = max(Psub(:), [], 'omitnan');
                    if ~isfinite(cl) || cl <= 0, cl = 1; end
                    cLo = 0;  cHi = cl;
                end
                
                cla(axMap);
                [F, W] = meshgrid(freq, wl);
                hSurf = surface(axMap, F, W, zeros(size(Cdat)), Cdat, ...
                    'EdgeColor', 'none', 'FaceColor', 'flat', ...
                    'HitTest', 'off', ...
                    'PickableParts', 'none');  %#ok<NASGU>
                view(axMap, 2);
                axMap.YDir = 'normal';
                axMap.CLim = [cLo, cHi];
                colormap(axMap, cmap);
                xlim(axMap, getXlimForFreq(freq));
                ylim(axMap, [min(wl),   max(wl)]);
                xlabel(axMap, freqAxLabel(cohResult.info.freqUnit));
                ylabel(axMap, app.xAxLabel());
                title(axMap, '|FFT|^2 map');
                cb = colorbar(axMap);
                cb.Label.String = cbLabel;
                % Re-disable interactivity AFTER cla / surface, since
                % MATLAB may quietly re-attach default interactions when
                % the axes content is rebuilt.
                try, disableDefaultInteractivity(axMap); catch, end
                try, axMap.Toolbar.Visible = 'off'; catch, end
                try, axMap.Interactions = []; catch, end
                
                % Initialize crosshair if first run
                if ~isfinite(curLambda)
                    curLambda = mean([min(wl), max(wl)]);
                end
                if ~isfinite(curFreq)
                    [~, iC]  = max(max(Pmap, [], 1));
                    curFreq  = freq(iC);
                end
                xhV = xline(axMap, curFreq, 'k-', 'LineWidth', 1.0, ...
                    'Alpha', 0.7, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                xhH = yline(axMap, curLambda, 'k-', 'LineWidth', 1.0, ...
                    'Alpha', 0.7, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                % Avg-range bounds
                yhAvgLo = yline(axMap, edAvgWLmin.Value, ':', ...
                    'Color', [0.30 0.30 0.30], 'LineWidth', 1.2, ...
                    'Alpha', 0.8, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                yhAvgHi = yline(axMap, edAvgWLmax.Value, ':', ...
                    'Color', [0.30 0.30 0.30], 'LineWidth', 1.2, ...
                    'Alpha', 0.8, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                
                % Click handler: surface is now non-interactive
                % (HitTest = off, PickableParts = none), so clicks pass
                % through to the axes itself which fires onMapClick.
                axMap.ButtonDownFcn = @(~,ev) onMapClick(ev);
                
                update1Ds();
            end
            
            function onMapClick(ev)
                pt = ev.IntersectionPoint;
                if isempty(pt), return; end
                % IMPORTANT: set the drag flag and start the safety
                % timer BEFORE applyCohPoint runs. applyCohPoint ends
                % with drawnow which yields the callback queue; if the
                % user has already released the mouse by then, the
                % WindowButtonUpFcn callback would fire there and set
                % draggingC = false. If we set draggingC = true AFTER
                % that yield we'd silently overwrite the release and
                % the drag would stay stuck. Setting the flag first
                % lets the up callback win when it should.
                draggingC = true;
                startCohDragTimer();
                applyCohPoint(pt);
            end
            function onCohMotion()
                if ~draggingC, return; end
                cp = axMap.CurrentPoint;
                if isempty(cp) || size(cp,1) < 1, return; end
                applyCohPoint(cp(1, 1:2));
            end
            function onCohUp()
                draggingC = false;
                stopCohDragTimer();
            end
            function startCohDragTimer()
                if isempty(cohDragTimer) || ~isvalid(cohDragTimer)
                    try
                        cohDragTimer = timer( ...
                            'ExecutionMode', 'fixedRate', ...
                            'Period', 0.1, 'BusyMode', 'drop', ...
                            'TimerFcn', @(~,~) checkCohDrag());
                    catch
                        cohDragTimer = [];  return;
                    end
                end
                try
                    if strcmp(cohDragTimer.Running, 'off')
                        start(cohDragTimer);
                    end
                catch
                end
            end
            function stopCohDragTimer()
                if ~isempty(cohDragTimer) && isvalid(cohDragTimer)
                    try
                        if strcmp(cohDragTimer.Running, 'on')
                            stop(cohDragTimer);
                        end
                    catch
                    end
                end
            end
            function checkCohDrag()
                % Runs every 100 ms while a drag is active. If the OS
                % mouse cursor has left the figure window we end the
                % drag here, even if no WindowButtonMotion / Up event
                % ever arrived.
                if ~draggingC
                    stopCohDragTimer();  return;
                end
                if ~app.mouseInsideFigure(f)
                    onCohUp();
                end
            end
            function applyCohPoint(pt)
                if numel(pt) < 2, return; end
                wl  = app.wavelength;
                fr  = cohResult.freq;
                if isempty(fr), return; end
                curFreq   = min(max(pt(1), min(fr)), max(fr));
                curLambda = min(max(pt(2), min(wl)), max(wl));
                if ~isempty(xhV) && isvalid(xhV), xhV.Value = curFreq;   end
                if ~isempty(xhH) && isvalid(xhH), xhH.Value = curLambda; end
                update1Ds();
                drawnow limitrate;
            end
            
            function update1Ds()
                if isempty(cohResult.Pmap), return; end
                wl   = app.wavelength(:);
                freq = cohResult.freq(:);
                Pmap = cohResult.Pmap;
                
                [~, iw] = min(abs(wl - curLambda));
                [~, ifq] = min(abs(freq - curFreq));
                wlAct  = wl(iw);
                fqAct  = freq(ifq);
                
                % Power spectrum at clicked lambda (top of right panel)
                cla(axPSpec);
                hold(axPSpec, 'on');
                plot(axPSpec, freq, Pmap(iw, :), 'b-', 'LineWidth', 1.4);
                xline(axPSpec, fqAct, 'k--', 'LineWidth', 0.8, 'Alpha', 0.6);
                hold(axPSpec, 'off');
                xlim(axPSpec, getXlimForFreq(freq));
                xlabel(axPSpec, '');                    % shared with avg below
                ylabel(axPSpec, '|FFT|^2');
                title(axPSpec, sprintf( ...
                    'Power spectrum at %s = %.2f %s', ...
                    app.xSym(), wlAct, app.xUnitAx()));
                set(axPSpec, 'YScale', ddYScale.Value);
                set(axPSpec, 'XTickLabel', []);
                grid(axPSpec, 'on');
                
                % Average over wavelength range (bottom of right panel)
                cla(axPSpecAvg);
                wlLo = min(edAvgWLmin.Value, edAvgWLmax.Value);
                wlHi = max(edAvgWLmin.Value, edAvgWLmax.Value);
                inRange = (wl >= wlLo) & (wl <= wlHi);
                if any(inRange)
                    Pavg = mean(Pmap(inRange, :), 1, 'omitnan');
                else
                    Pavg = nan(1, numel(freq));
                end
                hold(axPSpecAvg, 'on');
                plot(axPSpecAvg, freq, Pavg, '-', ...
                    'Color', [0.20 0.55 0.40], 'LineWidth', 1.4);
                xline(axPSpecAvg, fqAct, 'k--', 'LineWidth', 0.8, ...
                    'Alpha', 0.6);
                hold(axPSpecAvg, 'off');
                xlim(axPSpecAvg, getXlimForFreq(freq));
                xlabel(axPSpecAvg, freqAxLabel(cohResult.info.freqUnit));
                ylabel(axPSpecAvg, '|FFT|^2');
                title(axPSpecAvg, sprintf( ...
                    'Average over %s = %.0f .. %.0f %s (N=%d)', ...
                    app.xSym(), wlLo, wlHi, app.xUnitAx(), nnz(inRange)));
                set(axPSpecAvg, 'YScale', ddYScale.Value);
                grid(axPSpecAvg, 'on');
                
                % Spectral profile at wavenumber: top |FFT|^2, bottom phase
                cla(axProf);
                hold(axProf, 'on');
                plot(axProf, wl, Pmap(:, ifq), 'b-', 'LineWidth', 1.4);
                xline(axProf, wlAct, 'k--', 'LineWidth', 0.8, 'Alpha', 0.6);
                hold(axProf, 'off');
                xlim(axProf, [min(wl), max(wl)]);
                xlabel(axProf, '');
                ylabel(axProf, '|FFT|^2');
                title(axProf, sprintf( ...
                    'Spectral profile at \\nu = %.1f %s', fqAct, ...
                    cohResult.info.freqUnit));
                set(axProf, 'XTickLabel', []);
                grid(axProf, 'on');
                
                % Phase profile (lower half)
                cla(axProfPh);
                if isfield(cohResult.info, 'Phimap') && ...
                        ~isempty(cohResult.info.Phimap)
                    Phimap = cohResult.info.Phimap;
                    phaseCol = Phimap(:, ifq);
                    
                    % Mask phase where amplitude is below 5% of the
                    % peak in the same wavenumber column. The masked
                    % entries become NaN and so do not draw, while
                    % the high-amplitude phase reads cleanly.
                    pkCol = max(Pmap(:, ifq), [], 'omitnan');
                    if isfinite(pkCol) && pkCol > 0
                        mask = Pmap(:, ifq) < 0.05 * pkCol;
                        phaseCol(mask) = NaN;
                    end
                    
                    hold(axProfPh, 'on');
                    plot(axProfPh, wl, phaseCol, '-', ...
                        'Color', [0.55 0.30 0.70], 'LineWidth', 1.4);
                    yline(axProfPh, 0, 'k:');
                    xline(axProfPh, wlAct, 'k--', 'LineWidth', 0.8, ...
                        'Alpha', 0.6);
                    hold(axProfPh, 'off');
                    xlim(axProfPh, [min(wl), max(wl)]);
                    ylim(axProfPh, [-pi, pi]);
                    set(axProfPh, 'YTick', [-pi, -pi/2, 0, pi/2, pi]);
                    set(axProfPh, 'YTickLabel', ...
                        {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
                end
                xlabel(axProfPh, app.xAxLabel());
                ylabel(axProfPh, 'Phase (rad)');
                grid(axProfPh, 'on');
                
                % Residual time trace at lambda (within window)
                cla(axRes);
                hold(axRes, 'on');
                plot(axRes, cohResult.tWin, ...
                    cohResult.RinWin(iw, :), 'k-', 'LineWidth', 1.0);
                yline(axRes, 0, 'k:');
                hold(axRes, 'off');
                xlabel(axRes, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axRes, 'Residual');
                title(axRes, sprintf( ...
                    'Residual at %s = %.2f %s', ...
                    app.xSym(), wlAct, app.xUnitAx()));
                grid(axRes, 'on');
            end
            
            function doReset()
                edTMin.Value  = defaultTMin();
                edTMax.Value  = defaultTMax();
                ddResid.Value = 'ga';
                edNexp.Value  = 2;
                syncNexpEnable();
                ddApod.Value  = 'hann';
                edZpad.Value  = 2;
                ddDet.Value   = 'linear';
                ddFreqU.Value = 'cm-1';
                ddNorm.Value  = 'peak';
                ddYScale.Value = 'linear';
                edAvgWLmin.Value = min(app.wavelength);
                edAvgWLmax.Value = max(app.wavelength);
                edDispVmin.Value = 0;
                edDispVmax.Value = 0;
                cohResult = struct('Pmap', [], 'freq', [], ...
                    'Rresid', [], 'tWin', [], 'info', struct());
                curFreq = NaN;
                xhV = [];  xhH = [];
                yhAvgLo = [];  yhAvgHi = [];
                btnExpMap.Enable  = 'off';
                btnExpSpec.Enable = 'off';
                btnExpAvg.Enable  = 'off';
                btnExpProf.Enable = 'off';
                btnTF.Enable      = 'off';
                app.tfSaved = [];
                if ~isempty(app.tfFig) && isvalid(app.tfFig)
                    try, close(app.tfFig); catch, end
                end
                cla(axMap);  cla(axPSpec);  cla(axPSpecAvg);
                cla(axProf);  cla(axProfPh);  cla(axRes);
                title(axMap, '|FFT|^2 map');
                title(axPSpec, 'Power spectrum at clicked \lambda');
                title(axPSpecAvg, 'Average over \lambda range');
                title(axProf, 'Spectral profile at clicked \nu');
                title(axRes, 'Residual time trace at clicked \lambda');
                lblStatus.Text  = 'Reset.';
                txtResult.Value = {'Click "Run Coherence" to compute the FFT map.'};
                app.cohSaved = [];
            end
            
            function exportMap()
                if isempty(cohResult.Pmap), return; end
                [file, path] = app.putFile( ...
                    {'*.csv','CSV (*.csv)'; '*.dat','DAT (tab, *.dat)'}, ...
                    app.uiStr('Export |FFT|^2 map'), 'coherence_map.csv');
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    freq = cohResult.freq(:);          % Nfreq x 1
                    wl   = app.wavelength(:).';        % 1 x Nwl
                    % Transposed layout: rows index frequency, columns
                    % index wavelength.  First column = frequency axis,
                    % header row = wavelengths.
                    matT = [freq, cohResult.Pmap.'];   % Nfreq x (1 + Nwl)
                    fid = fopen(full, 'w');
                    if fid < 0
                        error('Could not open file for writing: %s', full);
                    end
                    cleanup = onCleanup(@() fclose(fid));
                    fprintf(fid, 'wavenumber_%s', cohResult.info.freqUnit);
                    for ii = 1:numel(wl)
                        fprintf(fid, '%s%.6g', dlm, wl(ii));
                    end
                    fprintf(fid, '\n');
                    clear cleanup;     % flush + close before writematrix
                    writematrix(matT, full, 'Delimiter', dlm, ...
                        'WriteMode', 'append');
                    uialert(f, sprintf('Saved: %s', full), ...
                        'Exported', 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Error');
                end
            end
            function exportSpec()
                if isempty(cohResult.Pmap), return; end
                wl = app.wavelength(:);
                [~, iw] = min(abs(wl - curLambda));
                hdr = {sprintf('nu_%s', cohResult.info.freqUnit), ...
                       'power'};
                mat = [cohResult.freq(:), cohResult.Pmap(iw, :).'];
                exportTbl(mat, hdr, ...
                    sprintf('coherence_Pnu_%.0f%s', wl(iw), app.xUnitTxt()), ...
                    'Export power spectrum');
            end
            function exportAvg()
                % Power spectrum averaged over the current lambda
                % range (the dashed band on the heatmap and the
                % lower-right panel).
                if isempty(cohResult.Pmap), return; end
                wl = app.wavelength(:);
                wlLo = min(edAvgWLmin.Value, edAvgWLmax.Value);
                wlHi = max(edAvgWLmin.Value, edAvgWLmax.Value);
                inRange = (wl >= wlLo) & (wl <= wlHi);
                if ~any(inRange)
                    uialert(f, ['No wavelengths fall inside the ' ...
                        'current averaging range.'], 'Empty range');
                    return;
                end
                Pavg = mean(cohResult.Pmap(inRange, :), 1, 'omitnan').';
                hdr = {sprintf('nu_%s', cohResult.info.freqUnit), ...
                       'power_avg'};
                mat = [cohResult.freq(:), Pavg];
                exportTbl(mat, hdr, ...
                    sprintf('coherence_avgPnu_%.0f-%.0f%s', ...
                        wlLo, wlHi, app.xUnitTxt()), ...
                    'Export averaged power spectrum');
            end
            function exportProf()
                if isempty(cohResult.Pmap), return; end
                freq = cohResult.freq(:);
                [~, ifq] = min(abs(freq - curFreq));
                if isfield(cohResult.info, 'Phimap') && ...
                        ~isempty(cohResult.info.Phimap)
                    hdr = {app.xHdr(), 'power', 'phase_rad'};
                    mat = [app.wavelength(:), ...
                        cohResult.Pmap(:, ifq), ...
                        cohResult.info.Phimap(:, ifq)];
                else
                    hdr = {app.xHdr(), 'power'};
                    mat = [app.wavelength(:), cohResult.Pmap(:, ifq)];
                end
                exportTbl(mat, hdr, ...
                    sprintf('coherence_profile_%.0f%s', ...
                        freq(ifq), cohResult.info.freqUnit), ...
                    'Export spectral profile');
            end
            
            function exportTbl(mat, hdr, defName, dlgTitle)
                [file, path] = app.putFile( ...
                    {'*.csv','CSV (*.csv)'; '*.dat','DAT (tab, *.dat)'}, ...
                    dlgTitle, [defName '.csv']);
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    % Variable names must be valid identifiers; sanitize.
                    hdrSan = matlab.lang.makeValidName(hdr, ...
                        'ReplacementStyle', 'underscore');
                    T = array2table(mat, 'VariableNames', hdrSan);
                    writetable(T, full, 'Delimiter', dlm);
                    uialert(f, sprintf('Saved: %s', full), ...
                        'Exported', 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Error');
                end
            end
            
            function cm = turboMap(n)
                % Approximate turbo colormap (no Image Processing Toolbox).
                % If newer MATLAB has built-in turbo(), prefer it.
                try
                    cm = turbo(n);  return;
                catch
                    cm = parula(n);   % graceful fallback
                end
            end
            
            function onClose()
                % Stop and delete the drag-safety timer if it is alive.
                if ~isempty(cohDragTimer) && isvalid(cohDragTimer)
                    try, stop(cohDragTimer); catch, end
                    try, delete(cohDragTimer); catch, end
                    cohDragTimer = [];
                end
                app.coherenceFig = [];
                delete(f);
            end
        end
        
        function openTimeFrequencyWindow(app)
            %OPENTIMEFREQUENCYWINDOW  nu-t map (STFT or Morlet CWT) of the
            %   coherence residual at one lambda or a lambda-range average.
            %   Input comes from the last Coherence run (app.cohSaved), so
            %   the residual method and analysis window are inherited.
            if isempty(app.cohSaved) || ~isstruct(app.cohSaved) || ...
                    ~isfield(app.cohSaved, 'result') || ...
                    ~isfield(app.cohSaved.result, 'Rresid') || ...
                    isempty(app.cohSaved.result.Rresid)
                uialert(app.fig, ['Run the Coherence analysis first -- ' ...
                    'its residual is the input of the time-frequency ' ...
                    'analysis.'], 'No residual');
                return;
            end
            R  = app.cohSaved.result.Rresid;
            cs = app.cohSaved.settings;
            if size(R, 1) ~= numel(app.wavelength) || ...
                    size(R, 2) ~= numel(app.delay)
                uialert(app.fig, ['The stored coherence residual no ' ...
                    'longer matches the current data (crop / reload). ' ...
                    'Re-run the Coherence analysis.'], 'Stale residual');
                return;
            end
            if ~isempty(app.tfFig) && isvalid(app.tfFig)
                figure(app.tfFig);  return;
            end

            f = uifigure('Name', 'Time-Frequency Analysis (STFT / Wavelet)', ...
                'Position', app.fitFigurePosition(1380, 880));
            app.applyLightTheme(f);
            app.tfFig = f;
            f.CloseRequestFcn       = @(~,~) onClose();
            f.WindowButtonMotionFcn = @(~,~) onTfMotion();
            f.WindowButtonUpFcn     = @(~,~) onTfUp();

            gl = uigridlayout(f, [1, 2]);
            gl.ColumnWidth = {345, '1x'};
            gl.Padding = [10 10 10 10];  gl.ColumnSpacing = 10;

            % ============== Left: setup ==============
            left = uipanel(gl, 'Title', 'Time-frequency setup');
            left.Layout.Column = 1;
            lgl = uigridlayout(left, [24, 1]);
            % rows: srcLbl, srcRow, lambdaRow, avgRow,
            %       winLbl, tMinRow, tMaxRow, detrendRow,
            %       methLbl, methRow, stftWinRow, stftHopRow, stftApodRow,
            %       stftZpadRow, cwtOm0Row, cwtNfRow,
            %       dispLbl, dispRow, zRow,
            %       btnRow, expLbl, expRow, statusLbl, resultArea
            lgl.RowHeight = {22, 30, 30, 30, ...
                             22, 30, 30, 30, ...
                             22, 30, 30, 30, 30, 30, 30, 30, ...
                             22, 30, 30, ...
                             32, 22, 30, 22, '1x'};
            lgl.Padding = [8 8 8 8];  lgl.RowSpacing = 6;

            wlLo0 = min(app.wavelength);  wlHi0 = max(app.wavelength);

            % --- Input trace ---
            uilabel(lgl, 'Text', 'Input trace (coherence residual):', ...
                'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Source:');
            ddSource = uidropdown(row, ...
                'Items',     app.uiStr({'Residual at \lambda', ...
                              '\lambda-range average'}), ...
                'ItemsData', {'single', 'avg'}, ...
                'Value', 'single', ...
                'ValueChangedFcn', @(~,~) onSourceChange());

            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('%s (%s):', app.xSymTxt(), app.xUnitTxt()));
            lam0 = app.cohSaved.curLambda;
            if isempty(lam0) || ~isfinite(lam0), lam0 = app.wavelength(round(end/2)); end
            edLambda = uieditfield(row, 'numeric', 'Value', lam0, ...
                'Limits', [wlLo0, wlHi0], ...
                'Tooltip', 'Wavelength whose residual trace is analysed.');

            row = uigridlayout(lgl, [1, 4]);
            row.ColumnWidth = {'fit', '1x', 'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 4;
            uilabel(row, 'Text', app.uiStr('Avg \lambda:'));
            edAvgLo = uieditfield(row, 'numeric', 'Value', cs.avgLo, ...
                'Limits', [wlLo0, wlHi0]);
            uilabel(row, 'Text', '~', 'HorizontalAlignment', 'center');
            edAvgHi = uieditfield(row, 'numeric', 'Value', cs.avgHi, ...
                'Limits', [wlLo0, wlHi0]);

            % --- Analysis window ---
            uilabel(lgl, 'Text', 'Analysis window:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('t_min (%s):', app.tUnitTxt()));
            edTMin = uieditfield(row, 'numeric', 'Value', cs.tmin);
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('t_max (%s):', app.tUnitTxt()));
            edTMax = uieditfield(row, 'numeric', 'Value', cs.tmax);
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Detrend:');
            ddDet = uidropdown(row, ...
                'Items',     {'None', 'Mean', 'Linear'}, ...
                'ItemsData', {'none', 'mean', 'linear'}, ...
                'Value', cs.detrend);

            % --- Method ---
            uilabel(lgl, 'Text', 'Method:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Transform:');
            ddMethod = uidropdown(row, ...
                'Items',     {'STFT (sliding-window FFT)', ...
                              'Wavelet (Morlet CWT)'}, ...
                'ItemsData', {'stft', 'cwt'}, ...
                'Value', 'stft', ...
                'Tooltip', ['STFT: fixed window -> uniform time / ' ...
                    'frequency resolution.' newline ...
                    'Wavelet: window scales with 1/nu -> frequency ' ...
                    'resolution at low nu, time resolution at high nu.'], ...
                'ValueChangedFcn', @(~,~) onMethodChange());

            winDef = (cs.tmax - cs.tmin) / 4;
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            lblWin = uilabel(row, 'Text', sprintf('STFT window (%s):', app.tUnitTxt()));
            edWin = uieditfield(row, 'numeric', 'Value', winDef, ...
                'Tooltip', ['Sliding-window length. Frequency resolution ' ...
                    '~ 1 / window; time resolution ~ window. Longer = ' ...
                    'sharper peaks, blurrier in time.']);
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            lblHop = uilabel(row, 'Text', sprintf('STFT hop (%s):', app.tUnitTxt()));
            edHop = uieditfield(row, 'numeric', 'Value', 0, ...
                'Tooltip', 'Step between window centres. 0 = every sample.');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            lblApod = uilabel(row, 'Text', 'STFT apodization:');
            ddApod = uidropdown(row, ...
                'Items',     {'Hann', 'Hamming', 'Blackman', 'Rectangular'}, ...
                'ItemsData', {'hann', 'hamming', 'blackman', 'rect'}, ...
                'Value', cs.apod);
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            lblZpad = uilabel(row, 'Text', 'STFT zero-pad factor:');
            edZpad = uieditfield(row, 'numeric', 'Value', cs.zpad, ...
                'Limits', [1 32], 'RoundFractionalValues', 'on');

            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            lblOm0 = uilabel(row, 'Text', app.uiStr('Wavelet \omega_0:'));
            edOm0 = uieditfield(row, 'numeric', 'Value', 6, ...
                'Limits', [2 40], ...
                'Tooltip', ['Morlet time-frequency trade-off. Higher = ' ...
                    'sharper frequency, blurrier time (6 is the usual ' ...
                    'default; 4-5 for fast dephasing, 10-12 for close ' ...
                    'modes).']);
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            lblNf = uilabel(row, 'Text', app.uiStr('Wavelet \nu grid points:'));
            edNf = uieditfield(row, 'numeric', 'Value', 200, ...
                'Limits', [16 2000], 'RoundFractionalValues', 'on');

            % --- Display ---
            uilabel(lgl, 'Text', sprintf('Display (%s):', app.uiStr(freqUnitTxt(cs.frequ))), ...
                'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 4]);
            row.ColumnWidth = {'fit', '1x', 'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 4;
            uilabel(row, 'Text', app.uiStr('\nu range:'));
            edDispVmin = uieditfield(row, 'numeric', 'Value', cs.dispLo, ...
                'ValueChangedFcn', @(~,~) onDispRangeChange());
            uilabel(row, 'Text', '~', 'HorizontalAlignment', 'center');
            edDispVmax = uieditfield(row, 'numeric', 'Value', cs.dispHi, ...
                'Tooltip', ['Visualization window on the frequency axis; ' ...
                    'also bounds the wavelet frequency grid at Run. ' ...
                    '0 = auto (up to Nyquist).'], ...
                'ValueChangedFcn', @(~,~) onDispRangeChange());
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Map z / 1D y scale:');
            ddZScale = uidropdown(row, ...
                'Items',     {'Linear', 'Log'}, ...
                'ItemsData', {'linear', 'log'}, ...
                'Value', 'linear', ...
                'ValueChangedFcn', @(~,~) onZScaleChange());

            % --- Buttons ---
            btnRow = uigridlayout(lgl, [1, 2]);
            btnRow.ColumnWidth = {'1x', '1x'};
            btnRow.Padding = [0 0 0 0];  btnRow.ColumnSpacing = 6;
            btnRun = uibutton(btnRow, 'Text', 'Run', ...
                'BackgroundColor', [0.25 0.55 0.75], 'FontColor', 'w', ...
                'FontWeight', 'bold');
            btnReset = uibutton(btnRow, 'Text', 'Reset');

            % --- Export ---
            uilabel(lgl, 'Text', 'Export:', 'FontWeight', 'bold');
            expRow = uigridlayout(lgl, [1, 3]);
            expRow.ColumnWidth = {'1x', '1x', '1x'};
            expRow.Padding = [0 0 0 0];  expRow.ColumnSpacing = 6;
            btnExpMap = uibutton(expRow, 'Text', app.uiStr('\nu-t map'), 'Enable', 'off');
            btnExpDec = uibutton(expRow, 'Text', app.uiStr('P(t)@\nu'), 'Enable', 'off', ...
                'Tooltip', 'Time decay of the power at the selected frequency.');
            btnExpSpc = uibutton(expRow, 'Text', app.uiStr('P(\nu)@t'), 'Enable', 'off', ...
                'Tooltip', 'Instantaneous spectrum at the selected time.');

            lblStatus = uilabel(lgl, 'Text', 'Ready.', 'FontColor', [0.2 0.2 0.2]);
            txtResult = uitextarea(lgl, 'Editable', 'off', ...
                'Value', {app.uiStr('Click "Run" to compute the \nu-t map.')}, ...
                'FontName', 'Courier New', 'FontSize', 10);

            % ============== Right: results ==============
            right = uipanel(gl, 'Title', 'Results');
            right.Layout.Column = 2;
            rgl = uigridlayout(right, [2, 2]);
            rgl.RowHeight   = {'3x', '2x'};
            rgl.ColumnWidth = {'3x', '2x'};
            rgl.Padding = [8 8 8 8];  rgl.RowSpacing = 8;  rgl.ColumnSpacing = 8;

            axMap = uiaxes(rgl);  axMap.Layout.Row = 1;  axMap.Layout.Column = 1;
            disableDefaultInteractivity(axMap);
            axMap.Toolbar.Visible = 'off';
            axMap.Interactions = [];
            title(axMap, '\nu-t map');

            axDecay = uiaxes(rgl);  axDecay.Layout.Row = 1;  axDecay.Layout.Column = 2;
            title(axDecay, 'Power vs time at selected \nu');
            axSpecT = uiaxes(rgl);  axSpecT.Layout.Row = 2;  axSpecT.Layout.Column = 1;
            title(axSpecT, 'Spectrum at selected t');
            axTrace = uiaxes(rgl);  axTrace.Layout.Row = 2;  axTrace.Layout.Column = 2;
            title(axTrace, 'Input trace (analysis window)');

            % ---- State ----
            tfResult = struct('Ptf', [], 'freq', [], 'tc', [], ...
                'info', struct(), 'trace', [], 'srcLabel', '');
            freqUnit = cs.frequ;
            curT  = NaN;  curNu = NaN;
            xhV = [];  xhH = [];
            draggingT   = false;
            tfDragTimer = [];

            % ---- Wire up ----
            btnRun.ButtonPushedFcn    = @(~,~) safeCb(@doRun,   'Run');
            btnReset.ButtonPushedFcn  = @(~,~) safeCb(@doReset, 'Reset');
            btnExpMap.ButtonPushedFcn = @(~,~) safeCb(@exportMap,   'Export map');
            btnExpDec.ButtonPushedFcn = @(~,~) safeCb(@exportDecay, 'Export decay');
            btnExpSpc.ButtonPushedFcn = @(~,~) safeCb(@exportSpecT, 'Export spectrum');
            axMap.ButtonDownFcn = @(~,ev) onMapClick(ev);
            onSourceChange();
            onMethodChange();

            % ---- Restore a previous analysis (if any) ----
            if ~isempty(app.tfSaved) && isstruct(app.tfSaved) && ...
                    isfield(app.tfSaved, 'result') && ...
                    ~isempty(app.tfSaved.result.Ptf) && ...
                    numel(app.tfSaved.result.trace) == numel(app.delay)
                try
                    st = app.tfSaved.settings;
                    ddSource.Value = st.source;  edLambda.Value = st.lambda;
                    edAvgLo.Value = st.avgLo;    edAvgHi.Value = st.avgHi;
                    edTMin.Value = st.tmin;      edTMax.Value = st.tmax;
                    ddDet.Value = st.detrend;    ddMethod.Value = st.method;
                    edWin.Value = st.winLen;     edHop.Value = st.hop;
                    ddApod.Value = st.apod;      edZpad.Value = st.zpad;
                    edOm0.Value = st.omega0;     edNf.Value = st.nFreq;
                    edDispVmin.Value = st.dispLo; edDispVmax.Value = st.dispHi;
                    ddZScale.Value = st.zscale;
                    tfResult = app.tfSaved.result;
                    curT  = app.tfSaved.curT;
                    curNu = app.tfSaved.curNu;
                    onSourceChange();  onMethodChange();
                    setExports('on');
                    plotResults();
                    lblStatus.Text  = sprintf('%s  [restored]', app.tfSaved.status);
                    txtResult.Value = app.tfSaved.lines;
                catch
                    lblStatus.Text = 'Stored result could not be restored -- run again.';
                end
            end

            % =====================  nested helpers  =====================

            function s = freqUnitTxt(u)
                switch lower(u)
                    case 'cm-1', s = 'cm^{-1}';
                    case 'thz',  s = 'THz';
                    otherwise,   s = 'Hz';
                end
            end
            function s = freqAxLabel(u)
                switch lower(u)
                    case 'cm-1', s = 'Wavenumber (cm^{-1})';
                    case 'thz',  s = 'Frequency (THz)';
                    otherwise,   s = 'Frequency (Hz)';
                end
            end

            function safeCb(fcn, label)
                try
                    fcn();
                catch ME
                    uialert(f, sprintf('"%s" failed:\n%s\n\n(In %s line %d)', ...
                        label, ME.message, ME.stack(1).name, ME.stack(1).line), ...
                        'Callback error');
                end
            end

            function onSourceChange()
                isSingle = strcmp(ddSource.Value, 'single');
                edLambda.Enable = onoff(isSingle);
                edAvgLo.Enable  = onoff(~isSingle);
                edAvgHi.Enable  = onoff(~isSingle);
            end
            function onMethodChange()
                isStft = strcmp(ddMethod.Value, 'stft');
                hs = {edWin, edHop, ddApod, edZpad};
                for k = 1:numel(hs), hs{k}.Enable = onoff(isStft); end
                hs = {lblWin, lblHop, lblApod, lblZpad};
                for k = 1:numel(hs), hs{k}.FontColor = pickColor(isStft); end
                hs = {edOm0, edNf};
                for k = 1:numel(hs), hs{k}.Enable = onoff(~isStft); end
                hs = {lblOm0, lblNf};
                for k = 1:numel(hs), hs{k}.FontColor = pickColor(~isStft); end
            end
            function s = onoff(b)
                if b, s = 'on'; else, s = 'off'; end
            end
            function c = pickColor(b)
                if b, c = [0.10 0.10 0.10]; else, c = [0.60 0.60 0.60]; end
            end
            function setExports(state)
                btnExpMap.Enable = state;
                btnExpDec.Enable = state;
                btnExpSpc.Enable = state;
            end

            function [lims, ok] = getDispLims()
                % [lo hi] in freqUnit; nu_max <= 0 -> up to the map's max.
                ok = false;  lims = [];
                if isempty(tfResult.freq), return; end
                fLo = min(tfResult.freq);  fHi = max(tfResult.freq);
                lo = edDispVmin.Value;  hi = edDispVmax.Value;
                if ~isfinite(hi) || hi <= 0, hi = fHi; end
                if ~isfinite(lo) || lo < 0, lo = fLo; end
                if lo > hi, tmp = lo; lo = hi; hi = tmp; end
                lo = max(lo, fLo);  hi = min(hi, fHi);
                if lo >= hi, return; end
                lims = [lo, hi];  ok = true;
            end

            % ---------- Run / Reset ----------
            function doRun()
                wl = app.wavelength(:);
                if strcmp(ddSource.Value, 'single')
                    [~, iw] = min(abs(wl - edLambda.Value));
                    y = R(iw, :);
                    edLambda.Value = wl(iw);
                    srcLabel = sprintf('%s = %.1f %s', app.xSym(), wl(iw), app.xUnitAx());
                else
                    lo = min(edAvgLo.Value, edAvgHi.Value);
                    hi = max(edAvgLo.Value, edAvgHi.Value);
                    m  = wl >= lo & wl <= hi;
                    if ~any(m)
                        uialert(f, 'No wavelengths inside the averaging range.', ...
                            'Empty range');
                        return;
                    end
                    y = mean(R(m, :), 1, 'omitnan');
                    srcLabel = sprintf('avg %s = %.0f..%.0f %s (N=%d)', ...
                        app.xSym(), lo, hi, app.xUnitAx(), nnz(m));
                end
                if edTMax.Value <= edTMin.Value
                    uialert(f, 't_max must exceed t_min.', 'Invalid window');
                    return;
                end

                opts = struct( ...
                    'tMin', edTMin.Value, 'tMax', edTMax.Value, ...
                    'detrend', ddDet.Value, 'method', ddMethod.Value, ...
                    'freqUnit', freqUnit, 'normMode', 'peak', ...
                    'winLen', edWin.Value, 'hop', edHop.Value, ...
                    'apod', ddApod.Value, 'zeroPad', edZpad.Value, ...
                    'omega0', edOm0.Value, 'nFreq', edNf.Value, ...
                    'fLo', edDispVmin.Value, 'fHi', edDispVmax.Value);

                lblStatus.Text = 'Computing...';  drawnow;
                tStart = tic;
                try
                    [Ptf, freq, tc, info] = app.computeTimeFrequency( ...
                        y, app.delay, opts);
                catch ME
                    uialert(f, sprintf('Time-frequency analysis failed:\n%s', ...
                        ME.message), 'Error');
                    lblStatus.Text = 'Failed.';
                    return;
                end
                elapsed = toc(tStart);

                tfResult = struct('Ptf', Ptf, 'freq', freq, 'tc', tc, ...
                    'info', info, 'trace', y, 'srcLabel', srcLabel);
                % Default cursor: strongest point of the map inside the
                % reliable region.
                Pm = Ptf;
                Pm(bsxfun(@lt, freq(:), info.coiFreq(:).')) = -Inf;
                [~, imax] = max(Pm(:));
                if isfinite(Pm(imax))
                    [iF, iT] = ind2sub(size(Pm), imax);
                else
                    [~, imax] = max(Ptf(:));  [iF, iT] = ind2sub(size(Ptf), imax);
                end
                curNu = freq(iF);  curT = tc(iT);

                setExports('on');
                plotResults();

                lines = {};
                lines{end+1} = sprintf('Source  : %s', app.uiStr(srcLabel));
                lines{end+1} = sprintf('Window  : %.4g .. %.4g %s', ...
                    edTMin.Value, edTMax.Value, app.tUnitTxt());
                lines{end+1} = sprintf('dt      : %.4g %s (%d pts)', ...
                    info.dtUni, app.tUnitTxt(), numel(info.tUni));
                lines{end+1} = sprintf('Detrend : %s', ddDet.Value);
                if strcmp(info.method, 'stft')
                    lines{end+1} = 'Method  : STFT';
                    lines{end+1} = sprintf('  window = %.4g %s (%d pts, %s)', ...
                        info.winLen, app.tUnitTxt(), info.nWin, info.apod);
                    lines{end+1} = sprintf('  hop    = %d pt, zero-pad x%d', ...
                        info.hopN, info.zeroPad);
                    lines{end+1} = sprintf('  d(nu)  = %.3g %s (FFT bin)', ...
                        freq(2) - freq(1), freqUnit);
                    lines{end+1} = sprintf('  edge-affected: t < %.4g, t > %.4g', ...
                        info.coiT(1), info.coiT(2));
                else
                    lines{end+1} = 'Method  : Morlet CWT';
                    lines{end+1} = sprintf('  omega0 = %.3g, %d freq pts', ...
                        info.omega0, info.nFreq);
                    lines{end+1} = sprintf('  nu grid: %.4g .. %.4g %s', ...
                        min(freq), max(freq), freqUnit);
                    lines{end+1} = '  Dashed line on map = cone of influence';
                    lines{end+1} = '  (below it: edge-affected, unreliable).';
                end
                lines{end+1} = '';
                lines{end+1} = 'Click + drag on the map to move the';
                lines{end+1} = 'cursor (nu, t).';
                txtResult.Value = lines;
                lblStatus.Text = sprintf('Done (%d x %d, %.2fs).', ...
                    size(Ptf, 1), size(Ptf, 2), elapsed);

                % Persist so the analysis survives closing this window.
                app.tfSaved = struct( ...
                    'result', tfResult, 'curT', curT, 'curNu', curNu, ...
                    'status', lblStatus.Text, 'lines', {txtResult.Value}, ...
                    'settings', struct( ...
                        'source', ddSource.Value, 'lambda', edLambda.Value, ...
                        'avgLo', edAvgLo.Value, 'avgHi', edAvgHi.Value, ...
                        'tmin', edTMin.Value, 'tmax', edTMax.Value, ...
                        'detrend', ddDet.Value, 'method', ddMethod.Value, ...
                        'winLen', edWin.Value, 'hop', edHop.Value, ...
                        'apod', ddApod.Value, 'zpad', edZpad.Value, ...
                        'omega0', edOm0.Value, 'nFreq', edNf.Value, ...
                        'dispLo', edDispVmin.Value, 'dispHi', edDispVmax.Value, ...
                        'zscale', ddZScale.Value));
            end

            function doReset()
                ddSource.Value = 'single';
                lam = app.cohSaved.curLambda;
                if isempty(lam) || ~isfinite(lam), lam = app.wavelength(round(end/2)); end
                edLambda.Value = lam;
                edAvgLo.Value = cs.avgLo;  edAvgHi.Value = cs.avgHi;
                edTMin.Value = cs.tmin;    edTMax.Value = cs.tmax;
                ddDet.Value = cs.detrend;  ddMethod.Value = 'stft';
                edWin.Value = winDef;      edHop.Value = 0;
                ddApod.Value = cs.apod;    edZpad.Value = cs.zpad;
                edOm0.Value = 6;           edNf.Value = 200;
                edDispVmin.Value = cs.dispLo;  edDispVmax.Value = cs.dispHi;
                ddZScale.Value = 'linear';
                onSourceChange();  onMethodChange();
                tfResult = struct('Ptf', [], 'freq', [], 'tc', [], ...
                    'info', struct(), 'trace', [], 'srcLabel', '');
                curT = NaN;  curNu = NaN;  xhV = [];  xhH = [];
                setExports('off');
                cla(axMap);  app.claAll(axDecay);  app.claAll(axSpecT);  app.claAll(axTrace);
                title(axMap, '\nu-t map');
                title(axDecay, 'Power vs time at selected \nu');
                title(axSpecT, 'Spectrum at selected t');
                title(axTrace, 'Input trace (analysis window)');
                lblStatus.Text  = 'Reset.';
                txtResult.Value = {app.uiStr('Click "Run" to compute the \nu-t map.')};
                app.tfSaved = [];
            end

            % ---------- Plotting ----------
            function plotResults()
                if isempty(tfResult.Ptf), return; end
                P    = tfResult.Ptf;  freq = tfResult.freq;  tc = tfResult.tc;
                info = tfResult.info;
                cmap = turboMap(256);

                [lims, ok] = getDispLims();
                if ~ok, lims = [min(freq), max(freq)]; end
                fMask = freq >= lims(1) & freq <= lims(2);
                if ~any(fMask), fMask = true(size(freq)); end
                Psub  = P(fMask, :);

                if strcmp(ddZScale.Value, 'log')
                    pk = max(Psub(:), [], 'omitnan');
                    if ~isfinite(pk) || pk <= 0, pk = 1; end
                    floorVal = pk * 1e-4;
                    Cdat = log10(max(P, floorVal));
                    cLo = log10(floorVal);  cHi = log10(pk);
                    cbLabel = 'log_{10} power';
                else
                    Cdat = P;
                    cHi = max(Psub(:), [], 'omitnan');
                    if ~isfinite(cHi) || cHi <= 0, cHi = 1; end
                    cLo = 0;  cbLabel = 'Power (norm.)';
                end

                cla(axMap);
                hold(axMap, 'on');
                hs = surface(axMap, tc, freq, zeros(numel(freq), numel(tc)), Cdat, ...
                    'EdgeColor', 'none', 'FaceColor', 'flat');
                hs.HitTest = 'off';  hs.PickableParts = 'none';
                view(axMap, 2);
                axMap.YDir = 'normal';
                axMap.CLim = [cLo, cHi];
                colormap(axMap, cmap);
                % Cone of influence
                if strcmp(info.method, 'cwt')
                    cf = info.coiFreq;  cf(cf > lims(2)) = NaN;
                    plot(axMap, tc, cf, 'w--', 'LineWidth', 1.2, ...
                        'HitTest', 'off', 'PickableParts', 'none');
                else
                    for e = info.coiT
                        xline(axMap, e, 'w--', 'LineWidth', 1.2, ...
                            'HitTest', 'off', 'PickableParts', 'none');
                    end
                end
                xhV = xline(axMap, curT,  'k-', 'LineWidth', 1.0, ...
                    'HitTest', 'off', 'PickableParts', 'none');
                xhH = yline(axMap, curNu, 'k-', 'LineWidth', 1.0, ...
                    'HitTest', 'off', 'PickableParts', 'none');
                hold(axMap, 'off');
                xlim(axMap, [min(tc), max(tc)]);
                ylim(axMap, lims);
                xlabel(axMap, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axMap, freqAxLabel(freqUnit));
                if strcmp(info.method, 'stft')
                    ttl = sprintf('STFT  (window %.3g %s)  --  %s', ...
                        info.winLen, app.tUnitAx(), tfResult.srcLabel);
                else
                    ttl = sprintf('Morlet CWT  (\\omega_0 = %.3g)  --  %s', ...
                        info.omega0, tfResult.srcLabel);
                end
                title(axMap, ttl);
                cb = colorbar(axMap);  cb.Label.String = cbLabel;
                % Re-apply interaction lock after cla (uiaxes resets it).
                disableDefaultInteractivity(axMap);
                axMap.Toolbar.Visible = 'off';
                axMap.Interactions = [];
                axMap.ButtonDownFcn = @(~,ev) onMapClick(ev);

                update1Ds();
            end

            function update1Ds()
                if isempty(tfResult.Ptf), return; end
                P = tfResult.Ptf;  freq = tfResult.freq;  tc = tfResult.tc;
                info = tfResult.info;
                useLog = strcmp(ddZScale.Value, 'log');
                [lims, ok] = getDispLims();
                if ~ok, lims = [min(freq), max(freq)]; end

                % --- Decay at nu ---
                [~, iF] = min(abs(freq - curNu));
                rowP = P(iF, :);
                reliable = info.coiFreq <= freq(iF);
                rOK = rowP;  rOK(~reliable) = NaN;
                rNO = rowP;  rNO(reliable)  = NaN;
                app.claAll(axDecay);  hold(axDecay, 'on');
                plot(axDecay, tc, rNO, '-', 'Color', [0.65 0.65 0.65], ...
                    'LineWidth', 1.2);
                plot(axDecay, tc, rOK, '-', 'Color', [0.10 0.30 0.80], ...
                    'LineWidth', 1.6);
                xline(axDecay, curT, 'k--', 'HandleVisibility', 'off');
                hold(axDecay, 'off');  grid(axDecay, 'on');
                xlim(axDecay, [min(tc), max(tc)]);
                if useLog, axDecay.YScale = 'log'; else, axDecay.YScale = 'linear'; end
                xlabel(axDecay, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axDecay, 'Power');
                title(axDecay, sprintf('P(t) at \\nu = %.1f %s   (gray = edge-affected)', ...
                    freq(iF), freqUnitTxt(freqUnit)));

                % --- Spectrum at t ---
                [~, iT] = min(abs(tc - curT));
                colP = P(:, iT);
                app.claAll(axSpecT);  hold(axSpecT, 'on');
                plot(axSpecT, freq, colP, '-', 'Color', [0.10 0.55 0.35], ...
                    'LineWidth', 1.6);
                xline(axSpecT, curNu, 'k--', 'HandleVisibility', 'off');
                cfT = info.coiFreq(iT);
                if isfinite(cfT) && cfT > lims(1) && cfT < lims(2)
                    xline(axSpecT, cfT, '--', 'Color', [0.5 0.5 0.5], ...
                        'Label', 'COI', 'LabelOrientation', 'horizontal', ...
                        'HandleVisibility', 'off');
                end
                hold(axSpecT, 'off');  grid(axSpecT, 'on');
                xlim(axSpecT, lims);
                if useLog, axSpecT.YScale = 'log'; else, axSpecT.YScale = 'linear'; end
                xlabel(axSpecT, freqAxLabel(freqUnit));
                ylabel(axSpecT, 'Power');
                if isfinite(cfT)
                    title(axSpecT, sprintf('P(\\nu) at t = %.4g %s', ...
                        tc(iT), app.tUnitAx()));
                else
                    title(axSpecT, sprintf('P(\\nu) at t = %.4g %s  (edge-affected)', ...
                        tc(iT), app.tUnitAx()));
                end

                % --- Input trace with the effective window ---
                tU = info.tUni;  yu = info.yu;
                app.claAll(axTrace);  hold(axTrace, 'on');
                yr = [min(yu), max(yu)];
                if ~(yr(2) > yr(1)), yr = yr + [-1 1]; end
                pad = 0.08 * (yr(2) - yr(1));  yr = yr + [-pad, pad];
                if strcmp(info.method, 'stft')
                    halfSpan = info.winLen / 2;
                else
                    sC = interp1(freq, info.scales, curNu, 'linear', 'extrap');
                    halfSpan = sqrt(2) * sC;          % e-folding half-width
                end
                tA = max(tU(1), curT - halfSpan);  tB = min(tU(end), curT + halfSpan);
                if tB > tA
                    patch(axTrace, [tA tB tB tA], [yr(1) yr(1) yr(2) yr(2)], ...
                        [0.95 0.75 0.30], 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
                        'HitTest', 'off', 'PickableParts', 'none');
                end
                plot(axTrace, tU, yu, 'k-', 'LineWidth', 1.0);
                yline(axTrace, 0, 'k:', 'HandleVisibility', 'off');
                xline(axTrace, curT, 'k--', 'HandleVisibility', 'off');
                hold(axTrace, 'off');  grid(axTrace, 'on');
                xlim(axTrace, [tU(1), tU(end)]);  ylim(axTrace, yr);
                xlabel(axTrace, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axTrace, 'Residual');
                if strcmp(info.method, 'stft')
                    title(axTrace, 'Input trace  (shaded = STFT window at cursor)');
                else
                    title(axTrace, 'Input trace  (shaded = wavelet e-folding width at cursor \nu)');
                end
            end

            function onZScaleChange()
                if isempty(tfResult.Ptf), return; end
                plotResults();
            end
            function onDispRangeChange()
                if isempty(tfResult.Ptf), return; end
                plotResults();
            end

            % ---------- Click + drag ----------
            function onMapClick(ev)
                pt = ev.IntersectionPoint;
                if isempty(pt), return; end
                % Flag + timer BEFORE the apply (which yields via drawnow):
                % a fast release must be allowed to win.
                draggingT = true;
                startTfDragTimer();
                applyTfPoint(pt);
            end
            function onTfMotion()
                if ~draggingT, return; end
                cp = axMap.CurrentPoint;
                if isempty(cp) || size(cp, 1) < 1, return; end
                applyTfPoint(cp(1, 1:2));
            end
            function onTfUp()
                draggingT = false;
                stopTfDragTimer();
            end
            function startTfDragTimer()
                if isempty(tfDragTimer) || ~isvalid(tfDragTimer)
                    try
                        tfDragTimer = timer('ExecutionMode', 'fixedRate', ...
                            'Period', 0.1, 'BusyMode', 'drop', ...
                            'TimerFcn', @(~,~) checkTfDrag());
                    catch
                        tfDragTimer = [];  return;
                    end
                end
                try
                    if strcmp(tfDragTimer.Running, 'off'), start(tfDragTimer); end
                catch
                end
            end
            function stopTfDragTimer()
                if ~isempty(tfDragTimer) && isvalid(tfDragTimer)
                    try
                        if strcmp(tfDragTimer.Running, 'on'), stop(tfDragTimer); end
                    catch
                    end
                end
            end
            function checkTfDrag()
                if ~draggingT, stopTfDragTimer();  return; end
                if ~app.mouseInsideFigure(f), onTfUp(); end
            end
            function applyTfPoint(pt)
                if numel(pt) < 2 || isempty(tfResult.Ptf), return; end
                tc = tfResult.tc;  fr = tfResult.freq;
                curT  = min(max(pt(1), min(tc)), max(tc));
                curNu = min(max(pt(2), min(fr)), max(fr));
                if ~isempty(xhV) && isvalid(xhV), xhV.Value = curT;  end
                if ~isempty(xhH) && isvalid(xhH), xhH.Value = curNu; end
                update1Ds();
                if ~isempty(app.tfSaved) && isstruct(app.tfSaved)
                    app.tfSaved.curT = curT;  app.tfSaved.curNu = curNu;
                end
                drawnow limitrate;
            end

            % ---------- Export ----------
            function exportMap()
                if isempty(tfResult.Ptf), return; end
                [file, path] = app.putFile({'*.csv', 'CSV (*.csv)'}, ...
                    'Export nu-t map', 'timefreq_map.csv');
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                try
                    % Header row = time axis; first column = frequency.
                    fid = fopen(full, 'w');
                    if fid < 0, error('Could not open file for writing: %s', full); end
                    fprintf(fid, 'nu_%s', freqUnit);
                    fprintf(fid, ',%.10g', tfResult.tc);
                    fprintf(fid, '\n');
                    fclose(fid);
                    writematrix([tfResult.freq(:), tfResult.Ptf], full, ...
                        'WriteMode', 'append');
                    uialert(f, sprintf('Saved: %s', full), 'Exported', ...
                        'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), 'Error');
                end
            end
            function exportDecay()
                if isempty(tfResult.Ptf), return; end
                [~, iF] = min(abs(tfResult.freq - curNu));
                hdr = {sprintf('t_%s', app.tUnitTxt()), 'power', 'reliable'};
                rel = double(tfResult.info.coiFreq(:) <= tfResult.freq(iF));
                mat = [tfResult.tc(:), tfResult.Ptf(iF, :).', rel];
                exportTbl(mat, hdr, sprintf('timefreq_Pt_%.0f', tfResult.freq(iF)), ...
                    'Export decay at nu');
            end
            function exportSpecT()
                if isempty(tfResult.Ptf), return; end
                [~, iT] = min(abs(tfResult.tc - curT));
                hdr = {sprintf('nu_%s', freqUnit), 'power'};
                mat = [tfResult.freq(:), tfResult.Ptf(:, iT)];
                exportTbl(mat, hdr, sprintf('timefreq_Pnu_t%.4g', tfResult.tc(iT)), ...
                    'Export spectrum at t');
            end
            function exportTbl(mat, hdr, defName, dlgTitle)
                [file, path] = app.putFile( ...
                    {'*.csv', 'CSV (*.csv)'; '*.dat', 'DAT (tab, *.dat)'}, ...
                    dlgTitle, [defName '.csv']);
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    hdrSan = matlab.lang.makeValidName(hdr, ...
                        'ReplacementStyle', 'underscore');
                    T = array2table(mat, 'VariableNames', hdrSan);
                    writetable(T, full, 'Delimiter', dlm);
                    uialert(f, sprintf('Saved: %s', full), 'Exported', ...
                        'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), 'Error');
                end
            end

            function cm = turboMap(n)
                try
                    cm = turbo(n);
                catch
                    cm = parula(n);
                end
            end

            function onClose()
                if ~isempty(tfDragTimer) && isvalid(tfDragTimer)
                    try, stop(tfDragTimer); catch, end
                    try, delete(tfDragTimer); catch, end
                    tfDragTimer = [];
                end
                app.tfFig = [];
                delete(f);
            end
        end

        function openCropWindow(app)
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');  return;
            end
            if ~isempty(app.cropFig) && isvalid(app.cropFig)
                figure(app.cropFig);  return;
            end
            
            f = uifigure('Name', 'Crop Data', 'Position', app.fitFigurePosition(820, 670));
            app.applyLightTheme(f);
            app.cropFig = f;
            
            gl = uigridlayout(f, [3, 1]);
            gl.RowHeight = {95, '1x', 44};
            gl.Padding = [10 10 10 10];  gl.RowSpacing = 8;
            
            % ---- Row 1: controls ----
            cgl = uigridlayout(gl, [2, 7]);
            cgl.Layout.Row = 1;
            cgl.RowHeight   = {30, 30};
            cgl.ColumnWidth = {'fit', 90, 20, 90, 'fit', 130, '1x'};
            cgl.Padding = [0 0 0 0];
            cgl.ColumnSpacing = 5;  cgl.RowSpacing = 3;
            
            % Use the CURRENT working data as both the default crop
            % range and the "Full" range. The raw originalDelay /
            % originalWavelength are kept only for "Reset all".
            wlMinVal = min(app.wavelength);   wlMaxVal = max(app.wavelength);
            tMinVal  = min(app.delay);        tMaxVal  = max(app.delay);
            wlMinFull = wlMinVal;             wlMaxFull = wlMaxVal;
            tMinFull  = tMinVal;              tMaxFull  = tMaxVal;
            
            uilabel(cgl, 'Text', sprintf('%s range (%s):', app.xSymTxt(), app.xUnitTxt()), ...
                'HorizontalAlignment', 'right');
            edWLmin = uieditfield(cgl, 'numeric', 'Value', wlMinVal);
            uilabel(cgl, 'Text', 'to', 'HorizontalAlignment', 'center');
            edWLmax = uieditfield(cgl, 'numeric', 'Value', wlMaxVal);
            uilabel(cgl, 'Text', '');   % spacer
            btnFullWL = uibutton(cgl, 'Text', app.uiStr('Full \lambda range'));
            uilabel(cgl, 'Text', '');   % spacer
            
            uilabel(cgl, 'Text', sprintf('t range (%s):', app.tUnitTxt()), ...
                'HorizontalAlignment', 'right');
            edTmin = uieditfield(cgl, 'numeric', 'Value', tMinVal);
            uilabel(cgl, 'Text', 'to', 'HorizontalAlignment', 'center');
            edTmax = uieditfield(cgl, 'numeric', 'Value', tMaxVal);
            uilabel(cgl, 'Text', '');
            btnFullT = uibutton(cgl, 'Text', 'Full t range');
            uilabel(cgl, 'Text', '');
            
            % ---- Row 2: preview ----
            ax = uiaxes(gl);
            ax.Layout.Row = 2;
            ax.Color = [0.93 0.93 0.93];
            xlabel(ax, app.xAxLabel());
            ylabel(ax, sprintf('Delay time (%s)', app.tUnitAx()));
            
            % ---- Row 3: action buttons ----
            bgl = uigridlayout(gl, [1, 4]);
            bgl.Layout.Row = 3;
            bgl.ColumnWidth = {'1x', 'fit', 'fit'};
            uilabel(bgl, 'Text', '');
            btnCancel = uibutton(bgl, 'Text', 'Cancel', ...
                'ButtonPushedFcn', @(~,~) close(f));   %#ok<NASGU>
            btnApply = uibutton(bgl, 'Text', 'Apply Crop', ...
                'BackgroundColor', [0.30 0.65 0.35], 'FontColor', 'w');
            
            edWLmin.ValueChangedFcn   = @(~,~) drawPreview();
            edWLmax.ValueChangedFcn   = @(~,~) drawPreview();
            edTmin.ValueChangedFcn    = @(~,~) drawPreview();
            edTmax.ValueChangedFcn    = @(~,~) drawPreview();
            btnFullWL.ButtonPushedFcn = @(~,~) setFullWL();
            btnFullT.ButtonPushedFcn  = @(~,~) setFullT();
            applyBusy = false;   % true while the button below is being processed
            btnApply.ButtonPushedFcn  = @(~,~) doApply();
            
            drawPreview();
            
            % ------- nested -------
            function drawPreview()
                wl1 = edWLmin.Value;  wl2 = edWLmax.Value;
                t1  = edTmin.Value;   t2  = edTmax.Value;
                if wl1 > wl2, tmp = wl1; wl1 = wl2; wl2 = tmp; end
                if t1  > t2,  tmp = t1;  t1  = t2;  t2  = tmp;  end
                
                cla(ax);  hold(ax, 'off');
                wlAx = app.wavelength(:).';
                tAx  = app.delay(:).';
                [W, T] = meshgrid(wlAx, tAx);
                C = app.deltaA.';
                h = surface(ax, W, T, zeros(size(C)), C, ...
                    'EdgeColor', 'none', 'FaceColor', 'flat');
                view(ax, 2);  axis(ax, 'tight');
                set(ax, 'YDir', 'normal');
                cl = max(abs(C(:)), [], 'omitnan');
                if isempty(cl) || cl == 0 || isnan(cl), cl = 1; end
                ax.CLim = [-cl cl];  colormap(ax, app.getColormap());
                cb = colorbar(ax);  cb.Label.String = '\DeltaA or Intensity';
                
                % alpha mask: 1 inside crop box, 0.2 outside
                mask = (W >= wl1) & (W <= wl2) & (T >= t1) & (T <= t2);
                alpha = 0.2 * ones(size(C));
                alpha(mask) = 1;
                h.AlphaData = alpha;
                h.AlphaDataMapping = 'none';
                h.FaceAlpha = 'flat';
                
                hold(ax, 'on');
                xline(ax, wl1, 'w-', 'LineWidth', 1.2);
                xline(ax, wl2, 'w-', 'LineWidth', 1.2);
                yline(ax, t1,  'w-', 'LineWidth', 1.2);
                yline(ax, t2,  'w-', 'LineWidth', 1.2);
                hold(ax, 'off');
                
                nWL = nnz(app.originalWavelength >= wl1 & ...
                          app.originalWavelength <= wl2);
                nT  = nnz(app.originalDelay >= t1 & ...
                          app.originalDelay <= t2);
                title(ax, sprintf(['Kept: %s = [%.1f, %.1f] %s  ' ...
                    '(%d pts)  |  t = [%.3g, %.3g] %s  (%d pts)'], ...
                    app.xSym(), wl1, wl2, app.xUnitAx(), nWL, t1, t2, app.tUnitAx(), nT));
                xlabel(ax, app.xAxLabel());
                ylabel(ax, sprintf('Delay time (%s)', app.tUnitAx()));
            end
            
            function setFullWL()
                edWLmin.Value = wlMinFull;
                edWLmax.Value = wlMaxFull;
                drawPreview();
            end
            function setFullT()
                edTmin.Value = tMinFull;
                edTmax.Value = tMaxFull;
                drawPreview();
            end
            function doApply()
                % The main window is redrawn while the crop is applied,
                % and a second click that arrives in the meantime would
                % be processed in the middle of the first: it would crop
                % again and close this window, and the first run would
                % then fail on the deleted window. Ignore that click.
                if applyBusy, return; end
                applyBusy = true;
                try
                    app.applyCropByRange( ...
                        edWLmin.Value, edWLmax.Value, ...
                        edTmin.Value,  edTmax.Value);
                catch errA
                    applyBusy = false;
                    rethrow(errA);
                end
                applyBusy = false;
                if isvalid(f), close(f); end
            end
        end
        
        % =================================================================
        % Load & Average window
        % =================================================================
        function openLoadAverageWindow(app, pathBase, fileNames)
            % ------ parse all files ------
            files = struct('name', {}, 'path', {}, 'wavelength', {}, ...
                'delay', {}, 'deltaA', {}, 'include', {}, 'maxAbs', {}, ...
                'shift', {}, 'isRef', {});
            for k = 1:numel(fileNames)
                try
                    [wl, t, A] = app.parseDataFile( ...
                        fullfile(pathBase, fileNames{k}));
                catch ME
                    uialert(app.fig, sprintf('Failed to read %s:\n%s', ...
                        fileNames{k}, ME.message), 'Load Error');
                    return;
                end
                files(end+1).name = fileNames{k};           %#ok<AGROW>
                files(end).path = fullfile(pathBase, fileNames{k});
                files(end).wavelength = wl;
                files(end).delay = t;
                files(end).deltaA = A;
                files(end).include = true;
                files(end).maxAbs = max(abs(A(:)), [], 'omitnan');
                files(end).shift = 0;
                files(end).isRef = (k == 1);   % first file is the default reference
            end
            
            if isempty(files)
                uialert(app.fig, 'No files were loaded.', 'Load Error');
                return;
            end
            
            szRef = size(files(1).deltaA);
            for k = 2:numel(files)
                if ~isequal(size(files(k).deltaA), szRef)
                    uialert(app.fig, sprintf( ...
                        ['Shape mismatch: %s is %dx%d, %s is %dx%d.\n' ...
                         'All files must share the same dimensions.'], ...
                        files(1).name, szRef(1), szRef(2), ...
                        files(k).name, size(files(k).deltaA,1), ...
                        size(files(k).deltaA,2)), 'Load Error');
                    return;
                end
            end
            
            % Keep the raw parse results so the transpose toggle can flip
            % back to the original orientation without re-reading disk.
            filesOrig = files;
            swapSuggested = app.detectAxisSwap( ...
                files(1).wavelength, files(1).delay);
            
            % ------ create window ------
            f = uifigure('Name', 'Load & Average Datasets', ...
                'Position', app.fitFigurePosition(1080, 780));
            app.applyLightTheme(f);
            app.loadFig = f;
            
            gl = uigridlayout(f, [6, 1]);
            gl.RowHeight = {'fit', 36, 170, 40, '1x', 44};
            gl.Padding   = [10 10 10 10];  gl.RowSpacing = 8;
            
            % Row 1: header label
            uilabel(gl, 'Text', ...
                app.uiStr(['Compare files visually, then check those you want ' ...
                 'to include in the average.  Unchecked files are ' ...
                 'shown dashed / semi-transparent.  Click + drag on ' ...
                 'the preview plots to move Ref \lambda / Ref t.  ' ...
                 'If the axes look swapped, toggle "Transpose axes".  ' ...
                 'To time-align scans, mark one file as Ref, ' ...
                 'enter Shift values manually or click ' ...
                 '"Auto-align to Ref".']), ...
                 'FontColor', [0.15 0.15 0.15]);
            
            % Row 2: file buttons + transpose
            btnPnl = uigridlayout(gl, [1, 6]);
            btnPnl.ColumnWidth = {'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
            btnPnl.Padding = [0 0 0 0];
            btnAdd = uibutton(btnPnl, 'Text', 'Add more files...');
            btnRem = uibutton(btnPnl, 'Text', 'Remove selected');
            btnAll = uibutton(btnPnl, 'Text', 'Check all');
            btnNone = uibutton(btnPnl, 'Text', 'Uncheck all');
            btnTrans = uibutton(btnPnl, 'state', ...
                'Text', 'Transpose axes', ...
                'Value', swapSuggested, ...
                'Tooltip', ['Swap wavelength and delay for all files. ' ...
                    'Use this when the source files were saved with ' ...
                    'rows = delay and columns = wavelength instead of ' ...
                    'the expected layout.']);
            uilabel(btnPnl, 'Text', '');
            
            % Row 3: file table
            tbl = uitable(gl);
            tbl.ColumnName     = {'File', app.uiStr('Max |\DeltaA|'), 'Include', ...
                sprintf('Shift (%s)', app.tUnitTxt()), 'Ref'};
            tbl.ColumnFormat   = {'char', 'char', 'logical', ...
                'numeric', 'logical'};
            tbl.ColumnEditable = [false false true true true];
            tbl.ColumnWidth    = {'auto', 110, 80, 110, 60};
            tbl.Data = buildTableData();
            tbl.CellEditCallback      = @(~,ev) onEdit(ev);
            tbl.CellSelectionCallback = @(~,ev) onSelect(ev);
            tbl.UserData.sel = [];
            
            % Row 4: ref lambda / t + alignment buttons
            refGl = uigridlayout(gl, [1, 9]);
            refGl.ColumnWidth = {'fit', 100, 'fit', 100, 'fit', 'fit', ...
                                 'fit', 'fit', 'fit'};
            refGl.Padding = [0 0 0 0];  refGl.ColumnSpacing = 6;
            
            wlAx0 = files(1).wavelength;
            tAx0  = files(1).delay;
            refWLinit = wlAx0(round(end/2));
            posIdx = find(tAx0 >= 1, 1);
            if isempty(posIdx), posIdx = round(numel(tAx0)/2); end
            refTinit = tAx0(posIdx);
            
            uilabel(refGl, 'Text', sprintf('Ref %s (%s):', app.xSymTxt(), app.xUnitTxt()), ...
                'HorizontalAlignment', 'right');
            edRefWL = uieditfield(refGl, 'numeric', 'Value', refWLinit, ...
                'Limits', [min(wlAx0), max(wlAx0)]);
            uilabel(refGl, 'Text', sprintf('Ref t (%s):', app.tUnitTxt()), ...
                'HorizontalAlignment', 'right');
            edRefT = uieditfield(refGl, 'numeric', 'Value', refTinit, ...
                'Limits', [min(tAx0), max(tAx0)]);
            uilabel(refGl, 'Text', 'Kinetics scale:', ...
                'HorizontalAlignment', 'right');
            ddKinScale = uidropdown(refGl, ...
                'Items', {'Linear', 'Log'}, 'Value', 'Log');
            btnAutoAlign = uibutton(refGl, 'Text', 'Auto-align to Ref', ...
                'Tooltip', app.uiStr(['Estimate the time shift of every other ' ...
                    'included file by cross-correlating the kinetic ' ...
                    'trace at "Ref \lambda" against the file marked ' ...
                    'as Ref. Results land in the Shift column.']));
            btnResetShift = uibutton(refGl, 'Text', 'Reset shifts', ...
                'Tooltip', 'Set all shifts back to 0.');
            cbShowLegend = uicheckbox(refGl, ...
                'Text', 'Show legend', 'Value', false, ...
                'Tooltip', ['Toggle the legend on the preview plots. ' ...
                    'When on, the legend shows short IDs ' ...
                    '(Scan 1, Scan 2, ...); full filenames remain in ' ...
                    'the table above.']);
            
            % Row 5: preview axes (2 side by side)
            pvGl = uigridlayout(gl, [1, 2]);
            pvGl.Padding = [0 0 0 0];  pvGl.ColumnSpacing = 8;
            axPrevS = uiaxes(pvGl);   axPrevS.Layout.Column = 1;
            axPrevK = uiaxes(pvGl);   axPrevK.Layout.Column = 2;
            
            % Click + drag inside the previews moves Ref lambda / Ref t.
            % Motion / up handlers are permanently attached and no-op
            % until a flag is set, matching the pattern used in other
            % drag-aware windows.
            draggingPrevS = false;
            draggingPrevK = false;
            axPrevS.ButtonDownFcn = @(~,ev) onPrevSClick(ev);
            axPrevK.ButtonDownFcn = @(~,ev) onPrevKClick(ev);
            f.WindowButtonMotionFcn = @(~,~) onPrevMotion();
            f.WindowButtonUpFcn     = @(~,~) onPrevUp();
            
            % Row 6: action buttons
            actGl = uigridlayout(gl, [1, 4]);
            actGl.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
            btnSave = uibutton(actGl, 'Text', 'Save Avg to CSV...');
            uilabel(actGl, 'Text', '');
            btnCancel = uibutton(actGl, 'Text', 'Cancel', ...
                'ButtonPushedFcn', @(~,~) close(f));  %#ok<NASGU>
            btnLoad = uibutton(actGl, 'Text', 'Load Average of Checked', ...
                'BackgroundColor', [0.30 0.65 0.35], 'FontColor', 'w');
            
            btnAdd.ButtonPushedFcn  = @(~,~) onAddMore();
            btnRem.ButtonPushedFcn  = @(~,~) onRemoveSelected();
            btnAll.ButtonPushedFcn  = @(~,~) onSetAll(true);
            btnNone.ButtonPushedFcn = @(~,~) onSetAll(false);
            btnTrans.ValueChangedFcn = @(src,~) onTransposeToggle(src.Value);
            edRefWL.ValueChangedFcn = @(~,~) drawPreview();
            edRefT.ValueChangedFcn  = @(~,~) drawPreview();
            ddKinScale.ValueChangedFcn = @(~,~) drawPreview();
            btnAutoAlign.ButtonPushedFcn  = @(~,~) onAutoAlign();
            btnResetShift.ButtonPushedFcn = @(~,~) onResetShifts();
            cbShowLegend.ValueChangedFcn  = @(~,~) drawPreview();
            btnSave.ButtonPushedFcn = @(~,~) onSaveAvg();
            applyBusy = false;   % true while the button below is being processed
            btnLoad.ButtonPushedFcn = @(~,~) onLoadAvg();
            
            % If heuristic suggested a swap, apply it now. This must run
            % after the reference fields and preview axes exist, since
            % onTransposeToggle updates them.
            if swapSuggested
                onTransposeToggle(true);
            else
                drawPreview();
            end
            
            % ============ nested helpers ============
            function onTransposeToggle(doSwap)
                % Swap or restore wavelength/delay axes of every file in
                % the preview list. filesOrig holds the raw parse results
                % so we always rebuild from a known baseline.
                for i = 1:numel(files)
                    wlOrig = filesOrig(i).wavelength;
                    tOrig  = filesOrig(i).delay;
                    Aorig  = filesOrig(i).deltaA;
                    if doSwap
                        wlNew = tOrig(:);
                        tNew  = wlOrig(:).';
                        Anew  = Aorig.';
                        % After transpose, what were columns become rows
                        % (wavelength) and vice versa. Sort both axes to
                        % preserve the monotonic convention used
                        % throughout the rest of the program.
                        [wlNew, iw] = sort(wlNew);  Anew = Anew(iw, :);
                        [tNew,  it] = sort(tNew);   Anew = Anew(:,  it);
                        files(i).wavelength = wlNew;
                        files(i).delay      = tNew;
                        files(i).deltaA     = Anew;
                    else
                        files(i).wavelength = wlOrig(:);
                        files(i).delay      = tOrig(:).';
                        files(i).deltaA     = Aorig;
                    end
                    files(i).maxAbs = max(abs(files(i).deltaA(:)), ...
                        [], 'omitnan');
                end
                % Refresh the shared reference axis fields.
                wlCur = files(1).wavelength;
                tCur  = files(1).delay;
                if isvalid(edRefWL)
                    edRefWL.Limits = [min(wlCur), max(wlCur)];
                    edRefWL.Value  = wlCur(round(end/2));
                end
                if isvalid(edRefT)
                    edRefT.Limits = [min(tCur), max(tCur)];
                    posIdx = find(tCur >= 1, 1);
                    if isempty(posIdx), posIdx = round(numel(tCur)/2); end
                    edRefT.Value = tCur(posIdx);
                end
                tbl.Data = buildTableData();
                drawPreview();
            end
            
            function D = buildTableData()
                D = cell(numel(files), 5);
                for i = 1:numel(files)
                    D{i,1} = files(i).name;
                    D{i,2} = sprintf('%.3e', files(i).maxAbs);
                    D{i,3} = files(i).include;
                    D{i,4} = files(i).shift;
                    D{i,5} = files(i).isRef;
                end
            end
            function refreshTable()
                tbl.Data = buildTableData();
            end
            function onEdit(ev)
                col = ev.Indices(2);
                row = ev.Indices(1);
                switch col
                    case 3   % Include
                        files(row).include = ev.NewData;
                        drawPreview();
                    case 4   % Shift
                        if isnumeric(ev.NewData) && isfinite(ev.NewData)
                            files(row).shift = ev.NewData;
                        else
                            files(row).shift = 0;
                            tbl.Data{row, 4} = 0;
                        end
                        drawPreview();
                    case 5   % Ref (single-select; checking one unchecks others)
                        if ev.NewData
                            for jj = 1:numel(files)
                                files(jj).isRef = (jj == row);
                            end
                        else
                            files(row).isRef = false;
                        end
                        refreshTable();
                end
            end
            function onSetAll(val)
                for i = 1:numel(files), files(i).include = val; end
                refreshTable();  drawPreview();
            end
            function onAddMore()
                [newFiles, newPath] = uigetfile( ...
                    {'*.csv;*.dat;*.txt', 'Data (*.csv, *.dat, *.txt)'; ...
                     '*.*', 'All Files'}, ...
                    'Select more TA data file(s)', pathBase, ...
                    'MultiSelect', 'on');
                if isequal(newFiles, 0), return; end
                app.rememberDir(newPath);
                if ~iscell(newFiles), newFiles = {newFiles}; end
                curTranspose = logical(btnTrans.Value);
                for i = 1:numel(newFiles)
                    try
                        [wl, t, A] = app.parseDataFile( ...
                            fullfile(newPath, newFiles{i}));
                    catch ME
                        uialert(f, sprintf('Failed to read %s:\n%s', ...
                            newFiles{i}, ME.message), 'Load Error');
                        continue;
                    end
                    if ~isequal(size(A), szRef)
                        uialert(f, sprintf( ...
                            'Shape mismatch with %s; skipping %s.', ...
                            files(1).name, newFiles{i}), 'Shape Error');
                        continue;
                    end
                    % Append raw parse to filesOrig (untransposed baseline).
                    filesOrig(end+1).name = newFiles{i};     %#ok<AGROW>
                    filesOrig(end).path = fullfile(newPath, newFiles{i});
                    filesOrig(end).wavelength = wl;
                    filesOrig(end).delay = t;
                    filesOrig(end).deltaA = A;
                    filesOrig(end).include = true;
                    filesOrig(end).maxAbs = max(abs(A(:)), [], 'omitnan');
                    filesOrig(end).shift = 0;
                    filesOrig(end).isRef = false;
                    % Append to files with current transpose state applied.
                    files(end+1).name = newFiles{i};         %#ok<AGROW>
                    files(end).path = fullfile(newPath, newFiles{i});
                    if curTranspose
                        wlN = t(:);  tN = wl(:).';  An = A.';
                        [wlN, iw] = sort(wlN);  An = An(iw, :);
                        [tN,  it] = sort(tN);   An = An(:,  it);
                        files(end).wavelength = wlN;
                        files(end).delay      = tN;
                        files(end).deltaA     = An;
                    else
                        files(end).wavelength = wl(:);
                        files(end).delay      = t(:).';
                        files(end).deltaA     = A;
                    end
                    files(end).include = true;
                    files(end).maxAbs  = max(abs(files(end).deltaA(:)), ...
                        [], 'omitnan');
                    files(end).shift = 0;
                    files(end).isRef = false;
                end
                refreshTable();  drawPreview();
            end
            function onRemoveSelected()
                sel = tbl.UserData.sel;
                if isempty(sel), return; end
                rows = unique(sel(:, 1));
                rows = rows(rows <= numel(files));
                if isempty(rows), return; end
                files(rows) = [];
                filesOrig(rows) = [];
                if isempty(files), close(f); return; end
                refreshTable();  drawPreview();
            end
            
            % ---- Preview click + drag (Ref lambda / Ref t control) ----
            function onPrevSClick(ev)
                % Set the drag flag first so a fast release that fires
                % during the apply step does not get overwritten.
                draggingPrevS = true;
                applyPrevSPoint(ev.IntersectionPoint);
            end
            function onPrevKClick(ev)
                draggingPrevK = true;
                applyPrevKPoint(ev.IntersectionPoint);
            end
            function onPrevMotion()
                if ~(draggingPrevS || draggingPrevK), return; end
                if ~app.mouseInsideFigure(f), onPrevUp(); return; end
                if draggingPrevS
                    cp = axPrevS.CurrentPoint;
                    if ~isempty(cp) && size(cp,1) >= 1
                        applyPrevSPoint(cp(1, 1:2));
                    end
                elseif draggingPrevK
                    cp = axPrevK.CurrentPoint;
                    if ~isempty(cp) && size(cp,1) >= 1
                        applyPrevKPoint(cp(1, 1:2));
                    end
                end
            end
            function onPrevUp()
                draggingPrevS = false;
                draggingPrevK = false;
            end
            function applyPrevSPoint(pt)
                % Spectrum panel: x = wavelength -> updates Ref lambda
                if numel(pt) < 1, return; end
                wlAx = files(1).wavelength;
                wl = min(max(pt(1), min(wlAx)), max(wlAx));
                edRefWL.Value = wl;
                drawPreview();
                drawnow limitrate;
            end
            function applyPrevKPoint(pt)
                % Kinetics panel: x = delay -> updates Ref t
                if numel(pt) < 1, return; end
                tAx = files(1).delay;
                t = min(max(pt(1), min(tAx)), max(tAx));
                edRefT.Value = t;
                drawPreview();
                drawnow limitrate;
            end
            
            function drawPreview()
                refWL = edRefWL.Value;  refT = edRefT.Value;
                wlAx = files(1).wavelength;
                tAx  = files(1).delay;
                [~, iW] = min(abs(wlAx - refWL));
                cmap = lines(max(8, numel(files)));
                
                % --- Spectrum @ t_ref ---
                % After shifting file i's time axis by +shift_i, the
                % column at refT corresponds to original time
                % (refT - shift_i) on that file. So we interpolate
                % each file's wavelength-row at that wavelength is not
                % needed; we just take the column on the shifted grid.
                cla(axPrevS);  hold(axPrevS, 'on');
                legH = gobjects(0);  legL = {};
                for i = 1:numel(files)
                    if files(i).include
                        style = '-';  width = 1.3;  alpha = 1;
                    else
                        style = '--'; width = 0.9;  alpha = 0.45;
                    end
                    col = [cmap(i,:), alpha];
                    % Spectrum at refT on the *shifted* time axis: that
                    % is the column at original time refT - shift.
                    tOrig = refT - files(i).shift;
                    [~, iT_i] = min(abs(files(i).delay - tOrig));
                    spec_i = files(i).deltaA(:, iT_i);
                    h = plot(axPrevS, wlAx, spec_i, style, ...
                        'Color', col, 'LineWidth', width);
                    legH(end+1) = h;                         %#ok<AGROW>
                    shortID = sprintf('Scan %d', i);
                    if files(i).isRef
                        legL{end+1} = sprintf('%s [REF]', shortID); %#ok<AGROW>
                    elseif files(i).shift ~= 0
                        legL{end+1} = sprintf('%s (shift=%+.3g)', ...
                            shortID, files(i).shift); %#ok<AGROW>
                    else
                        legL{end+1} = shortID;                   %#ok<AGROW>
                    end
                end
                yline(axPrevS, 0, 'k:');
                % Vertical line marking the current Ref lambda (which is
                % the wavelength tracked by the right kinetics panel).
                % Click + drag inside the spectrum to move it.
                xline(axPrevS, refWL, 'k-', 'LineWidth', 1.0, ...
                    'Alpha', 0.7, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                hold(axPrevS, 'off');
                xlabel(axPrevS, app.xAxLabel());
                ylabel(axPrevS, '\DeltaA or Intensity');
                grid(axPrevS, 'on');
                xlim(axPrevS, [min(wlAx), max(wlAx)]);
                title(axPrevS, sprintf('Spectrum at t = %.3g %s', ...
                    refT, app.tUnitAx()));
                if cbShowLegend.Value
                    legend(axPrevS, legH, legL, 'Location', 'best', ...
                        'FontSize', 8, 'Interpreter', 'none');
                else
                    legend(axPrevS, 'off');
                end
                
                % --- Kinetics @ lambda_ref ---
                % Plot each file's kinetic on its *shifted* time axis,
                % so visually-aligned curves indicate consistent shifts.
                cla(axPrevK);  hold(axPrevK, 'on');
                legH = gobjects(0);  legL = {};
                for i = 1:numel(files)
                    if files(i).include
                        style = '-';  width = 1.3;  alpha = 1;
                    else
                        style = '--'; width = 0.9;  alpha = 0.45;
                    end
                    col = [cmap(i,:), alpha];
                    [~, iW_i] = min(abs(files(i).wavelength - refWL));
                    tShift = files(i).delay + files(i).shift;
                    h = plot(axPrevK, tShift, ...
                        files(i).deltaA(iW_i, :), style, ...
                        'Color', col, 'LineWidth', width);
                    legH(end+1) = h;                         %#ok<AGROW>
                    shortID = sprintf('Scan %d', i);
                    if files(i).isRef
                        legL{end+1} = sprintf('%s [REF]', shortID); %#ok<AGROW>
                    elseif files(i).shift ~= 0
                        legL{end+1} = sprintf('%s (shift=%+.3g)', ...
                            shortID, files(i).shift); %#ok<AGROW>
                    else
                        legL{end+1} = shortID;                   %#ok<AGROW>
                    end
                end
                yline(axPrevK, 0, 'k:');
                % Vertical line marking the current Ref t. In log mode
                % only draw it when refT > 0 (the log axis hides
                % non-positive values).
                if ~strcmp(ddKinScale.Value, 'Log') || refT > 0
                    xline(axPrevK, refT, 'k-', 'LineWidth', 1.0, ...
                        'Alpha', 0.7, 'HitTest', 'off', ...
                        'PickableParts', 'none');
                end
                hold(axPrevK, 'off');
                xlabel(axPrevK, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axPrevK, '\DeltaA or Intensity');
                grid(axPrevK, 'on');
                if strcmp(ddKinScale.Value, 'Log')
                    set(axPrevK, 'XScale', 'log');
                    pos = tAx > 0;
                    if any(pos)
                        xlim(axPrevK, [min(tAx(pos)), max(tAx)]);
                    end
                else
                    set(axPrevK, 'XScale', 'linear');
                    xlim(axPrevK, [min(tAx), max(tAx)]);
                end
                title(axPrevK, sprintf('Kinetics at %s = %.1f %s', ...
                    app.xSym(), wlAx(iW), app.xUnitAx()));
                if cbShowLegend.Value
                    legend(axPrevK, legH, legL, 'Location', 'best', ...
                        'FontSize', 8, 'Interpreter', 'none');
                else
                    legend(axPrevK, 'off');
                end
            end
            
            function [wlAv, tAv, avgD, nUsed] = computeAverage()
                inclIdx = find([files.include]);
                nUsed = numel(inclIdx);
                if nUsed == 0
                    wlAv = [];  tAv = [];  avgD = [];  return;
                end
                
                % Use file 1's grid as the common output grid.
                wlAv = files(1).wavelength;
                tAv  = files(1).delay;
                tAv  = tAv(:).';
                Nwl  = numel(wlAv);
                Nt   = numel(tAv);
                
                % Stack each included file resampled onto (wlAv, tAv)
                % after applying its time shift. NaN where extrapolation
                % would be required, so the mean is taken over the valid
                % files only (omitnan).
                stack = nan(Nwl, Nt, nUsed);
                anyShift = false;
                for ii = 1:nUsed
                    fi  = files(inclIdx(ii));
                    shi = fi.shift;
                    if shi ~= 0, anyShift = true; end
                    
                    % After the shift, sample j of file i corresponds to
                    % time fi.delay(j) + shi on the common grid.
                    % Equivalently, to obtain the value of file i at
                    % time tAv(k) on the COMMON grid, we read file i at
                    % its ORIGINAL time tAv(k) - shi.
                    tQuery = tAv - shi;
                    
                    if isequal(fi.wavelength(:), wlAv(:)) ...
                            && shi == 0 ...
                            && isequal(fi.delay(:).', tAv)
                        % No interpolation needed -- exact match.
                        stack(:, :, ii) = fi.deltaA;
                        continue;
                    end
                    
                    % Interpolate per-wavelength row in time.
                    Atmp = nan(Nwl, Nt);
                    for jj = 1:Nwl
                        % Find the matching wavelength index on this
                        % file (most files share a grid, so this is
                        % typically just jj).
                        if isequal(fi.wavelength(:), wlAv(:))
                            yrow = fi.deltaA(jj, :);
                        else
                            % Different wavelength axes -- interpolate
                            % the full 2D field on demand.
                            [WLG, TG] = meshgrid(fi.wavelength, fi.delay);
                            % Note: this branch is rare; we keep it
                            % simple and slow.
                            Aint = interp2(WLG, TG, fi.deltaA.', ...
                                wlAv(jj) * ones(1, Nt), tQuery, ...
                                'linear', NaN);
                            Atmp(jj, :) = Aint;
                            continue;
                        end
                        Atmp(jj, :) = interp1(fi.delay, yrow, ...
                            tQuery, 'linear', NaN);
                    end
                    stack(:, :, ii) = Atmp;
                end
                
                avgD = mean(stack, 3, 'omitnan');
            end
            
            function onAutoAlign()
                % Cross-correlate each file's kinetic at the reference
                % wavelength against the file marked Ref, find the lag
                % that maximizes correlation, and store it as the file's
                % shift.  Sub-sample resolution is obtained by parabolic
                % interpolation of the correlation peak.
                refIdx = find([files.isRef], 1);
                if isempty(refIdx)
                    refIdx = 1;
                    files(refIdx).isRef = true;
                    refreshTable();
                end
                if ~files(refIdx).include
                    uialert(f, ['The Ref file is unchecked. Check it ' ...
                        'before running auto-align.'], 'Auto-align');
                    return;
                end
                refWL = edRefWL.Value;
                
                tAxC = files(refIdx).delay(:).';
                wlR  = files(refIdx).wavelength;
                [~, iWr] = min(abs(wlR - refWL));
                yref = files(refIdx).deltaA(iWr, :);
                
                % Resample reference to a uniform grid (xcorr needs
                % uniform sampling for the lag-to-time conversion).
                dts   = diff(tAxC);
                dtUni = median(dts);
                if dtUni <= 0 || ~isfinite(dtUni), dtUni = mean(dts); end
                tUni  = tAxC(1):dtUni:tAxC(end);
                yrefU = interp1(tAxC, yref, tUni, 'linear', 0);
                yrefU = detrendRow(yrefU);
                
                changed = 0;
                for ii = 1:numel(files)
                    if ii == refIdx, continue; end
                    if ~files(ii).include, continue; end
                    [~, iWi] = min(abs(files(ii).wavelength - refWL));
                    yi = files(ii).deltaA(iWi, :);
                    yiU = interp1(files(ii).delay, yi, tUni, ...
                        'linear', 0);
                    yiU = detrendRow(yiU);
                    
                    % Cross-correlation r(m) = sum_n yref(n + m) yi(n)
                    % for the lags m = -(n-1) ... (n-1), by FFT (xcorr
                    % would need the Signal Processing Toolbox).
                    nU   = numel(yrefU);
                    nfft = 2^nextpow2(2 * nU - 1);
                    r    = real(ifft(fft(yrefU(:), nfft) .* ...
                                     conj(fft(yiU(:), nfft))));
                    r    = [r(nfft - nU + 2 : nfft); r(1:nU)];
                    lags = -(nU - 1) : (nU - 1);
                    [~, iMax] = max(r);
                    iMax = max(2, min(numel(r) - 1, iMax));
                    % Parabolic interpolation around the peak for
                    % sub-sample lag accuracy.
                    rL = r(iMax-1);  r0 = r(iMax);  rR = r(iMax+1);
                    denom = (rL - 2*r0 + rR);
                    if abs(denom) > 0
                        delta = 0.5 * (rL - rR) / denom;
                    else
                        delta = 0;
                    end
                    lagSamples = lags(iMax) + delta;
                    shift_i = lagSamples * dtUni;
                    files(ii).shift = shift_i;
                    changed = changed + 1;
                end
                refreshTable();  drawPreview();
                lastMsg = sprintf( ...
                    'Auto-aligned %d file(s) to "%s".', changed, ...
                    files(refIdx).name);
                uialert(f, lastMsg, 'Auto-align', 'Icon', 'success');
            end
            function y = detrendRow(y)
                y = y - mean(y, 'omitnan');
                y(~isfinite(y)) = 0;
            end
            
            function onResetShifts()
                for ii = 1:numel(files)
                    files(ii).shift = 0;
                end
                refreshTable();  drawPreview();
            end
            
            function onSaveAvg()
                [wlAv, tAv, avgD, N] = computeAverage();
                if isempty(avgD)
                    uialert(f, 'Check at least one file.', 'Nothing to save');
                    return;
                end
                [fn, fp] = app.putFile({'*.csv', 'CSV (*.csv)'}, ...
                    'Save averaged data', 'averaged.csv');
                if isequal(fn, 0), return; end
                M = zeros(numel(wlAv)+1, numel(tAv)+1);
                M(1, 1) = 0;
                M(1, 2:end) = tAv;
                M(2:end, 1) = wlAv;
                M(2:end, 2:end) = avgD;
                try
                    writematrix(M, fullfile(fp, fn));
                    uialert(f, sprintf(['Saved average of %d files ' ...
                        'to:\n%s'], N, fn), 'Saved', 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Save failed:\n%s', ME.message), ...
                        'Save Error');
                end
            end
            
            function onSelect(ev)
                tbl.UserData.sel = ev.Indices;
            end
            
            function onLoadAvg()
                % A second click while the data is being loaded and the
                % main window redrawn is ignored (see doApply of the
                % crop window).
                if applyBusy, return; end
                [wlAv, tAv, avgD, N] = computeAverage();
                if isempty(avgD)
                    uialert(f, 'Check at least one file.', 'None selected');
                    return;
                end
                inclNames = {files([files.include]).name};
                if N == 1
                    desc = sprintf('Loaded: %s', inclNames{1});
                else
                    desc = sprintf('Avg of %d files (%s...)', ...
                        N, inclNames{1});
                end
                applyBusy = true;
                try
                    app.setLoadedData(wlAv, tAv, avgD, desc);
                catch errA
                    applyBusy = false;
                    rethrow(errA);
                end
                applyBusy = false;
                if isvalid(f), close(f); end
            end
        end
        
        % =================================================================
        % Background correction window
        % =================================================================
        function openBackgroundWindow(app)
            if isempty(app.deltaA_raw)
                uialert(app.fig, 'Load data first.', 'No data');  return;
            end
            if ~isempty(app.bgFig) && isvalid(app.bgFig)
                figure(app.bgFig);  return;
            end
            
            f = uifigure('Name', 'Background Correction', ...
                'Position', app.fitFigurePosition(720, 520));
            app.applyLightTheme(f);
            app.bgFig = f;
            
            gl = uigridlayout(f, [3, 3]);
            gl.RowHeight   = {40, '1x', 44};
            gl.ColumnWidth = {'fit', 120, '1x'};
            gl.Padding     = [10 10 10 10];
            
            uilabel(gl, 'Text', '# of initial delay points to average:', ...
                'HorizontalAlignment', 'right');
            defaultN = min(5, numel(app.delay));
            if app.bgApplied, defaultN = app.bgN; end
            edN = uieditfield(gl, 'numeric', ...
                'Limits', [1 numel(app.delay)], ...
                'RoundFractionalValues', 'on', ...
                'Value', defaultN);
            lblInfo = uilabel(gl, 'Text', '', ...
                'FontColor', [0.25 0.25 0.25]);
            
            ax = uiaxes(gl);
            ax.Layout.Row = 2;  ax.Layout.Column = [1 3];
            xlabel(ax, app.xAxLabel());  ylabel(ax, '\DeltaA or Intensity');
            title(ax, 'Early-time spectra and their mean (background)');
            
            ctrl = uigridlayout(gl, [1, 3]);
            ctrl.Layout.Row = 3;  ctrl.Layout.Column = [1 3];
            ctrl.ColumnWidth = {'1x', 'fit', 'fit'};
            uilabel(ctrl, 'Text', '');
            btnApply = uibutton(ctrl, 'Text', 'Apply & Close', ...
                'BackgroundColor', [0.30 0.65 0.35], 'FontColor', 'w');
            btnCancel = uibutton(ctrl, 'Text', 'Cancel', ...
                'ButtonPushedFcn', @(~,~) close(f));  %#ok<NASGU>
            
            edN.ValueChangedFcn      = @(~,~) drawPreview();
            applyBusy = false;   % true while the button below is being processed
            btnApply.ButtonPushedFcn = @(~,~) applyBG();
            drawPreview();
            
            function drawPreview()
                N = round(edN.Value);
                if N < 1 || N > numel(app.delay), return; end
                app.claAll(ax);  hold(ax, 'on');   % preview lines are hidden-handle
                cmap = parula(max(N, 2));
                for k = 1:N
                    plot(ax, app.wavelength, app.deltaA_raw(:, k), '-', ...
                        'Color', [cmap(k,:), 0.6], 'LineWidth', 0.7, ...
                        'HandleVisibility', 'off');
                end
                bg = mean(app.deltaA_raw(:, 1:N), 2, 'omitnan');
                hMean = plot(ax, app.wavelength, bg, 'k-', 'LineWidth', 2, ...
                    'DisplayName', sprintf('Mean of first %d point(s)', N));
                hold(ax, 'off');
                grid(ax, 'on');  yline(ax, 0, 'k:');
                legend(ax, hMean, 'Location', 'best');
                lblInfo.Text = sprintf('Using t = %.4g ... %.4g %s', ...
                    app.delay(1), app.delay(N), app.tUnitTxt());
            end
            
            function applyBG()
                % A second click while the main window is being redrawn
                % is ignored (see doApply of the crop window).
                if applyBusy, return; end
                N = round(edN.Value);
                if N < 1 || N > numel(app.delay), return; end
                applyBusy = true;
                try
                    app.bgSpectrum = mean(app.deltaA_raw(:, 1:N), 2, 'omitnan');
                    app.bgN        = N;
                    app.bgApplied  = true;
                    app.updateAll();
                catch errA
                    applyBusy = false;
                    rethrow(errA);
                end
                applyBusy = false;
                if isvalid(f), close(f); end
            end
        end
        
        % =================================================================
        % Ridge computation (chirp helper)
        % =================================================================
        function computeRidge(app, method, tMin, tMax, smoothN)
            base = app.getChirpBaseData();
            m = app.delay >= tMin & app.delay <= tMax;
            tSub = app.delay(m);
            if nnz(m) < 3
                app.ridgeWL = [];  app.ridgeT = [];  return;
            end
            sub = base(:, m);
            switch method
                case 'max|dA/dt|'
                    dA = diff(sub, 1, 2);
                    dt = diff(tSub);
                    rate = dA ./ dt;
                    tAxis = (tSub(1:end-1) + tSub(2:end)) / 2;
                    [~, idx] = max(abs(rate), [], 2);
                case 'max|deltaA|'
                    tAxis = tSub;
                    [~, idx] = max(abs(sub), [], 2);
                otherwise
                    app.ridgeWL = [];  app.ridgeT = [];  return;
            end
            app.ridgeWL = app.wavelength(:);
            app.ridgeT  = tAxis(idx).';
            if smoothN > 1 && numel(app.ridgeT) > smoothN
                app.ridgeT = smoothdata(app.ridgeT, 'movmean', smoothN);
            end
        end
        
        % =================================================================
        % Chirp correction window
        % =================================================================
        function openChirpWindow(app)
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');  return;
            end
            if ~isempty(app.chirpFig) && isvalid(app.chirpFig)
                figure(app.chirpFig);  return;
            end
            
            f = uifigure('Name', 'Chirp Correction', ...
                'Position', app.fitFigurePosition(1180, 780));
            app.applyLightTheme(f);
            app.chirpFig = f;
            f.WindowButtonMotionFcn = @(~,~) onChirpDrag();
            f.WindowButtonUpFcn     = @(~,~) onChirpDragUp();
            
            gl = uigridlayout(f, [4, 2]);
            gl.RowHeight   = {35, '1x', 95, 44};
            gl.ColumnWidth = {'3x', '2x'};
            gl.RowSpacing  = 8;  gl.ColumnSpacing = 8;
            gl.Padding     = [10 10 10 10];
            
            lblHdr = uilabel(gl, 'Text', ...
                app.uiStr(['Click + drag on the 2D map to position the (\lambda, t_0) ' ...
                 'crosshair, then press "Save point" to commit it ' ...
                 '(>= 4 points needed for fit).  ' ...
                 'Model:  t_0 = a \cdot sqrt( (b w^2 - 1)/(c w^2 - 1) ) + d.']), ...
                 'FontColor', [0.15 0.15 0.15]);
            lblHdr.Layout.Row = 1;  lblHdr.Layout.Column = [1 2];
            
            axMap = uiaxes(gl);
            axMap.Layout.Row = 2;  axMap.Layout.Column = 1;
            axMap.Color = [0.93 0.93 0.93];
            xlabel(axMap, app.xAxLabel());
            ylabel(axMap, sprintf('Delay time (%s)', app.tUnitAx()));
            
            rightPnl = uipanel(gl, 'BorderType', 'none');
            rightPnl.Layout.Row = 2;  rightPnl.Layout.Column = 2;
            rgl = uigridlayout(rightPnl, [3, 1]);
            rgl.RowHeight   = {'1x', '1x', '1x'};
            rgl.Padding     = [0 0 0 0];  rgl.RowSpacing = 8;
            
            axFit = uiaxes(rgl);
            xlabel(axFit, app.xAxLabel());  ylabel(axFit, sprintf('t_0 (%s)', app.tUnitAx()));
            title(axFit, 'Polynomial fit preview');
            grid(axFit, 'on');
            
            tbl = uitable(rgl, ...
                'ColumnName', {sprintf('%s (%s)', app.xSymTxt(), app.xUnitTxt()), ...
                               sprintf('t (%s)', app.tUnitTxt())}, ...
                'ColumnEditable', [false false]);
            
            axKinP = uiaxes(rgl);
            xlabel(axKinP, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axKinP, '\DeltaA or Intensity');
            title(axKinP, 'Kinetics at last-clicked \lambda');
            grid(axKinP, 'on');
            
            ridgePnl = uipanel(gl, 'Title', 'Ridge Helper (reference t_0 curve)');
            ridgePnl.Layout.Row = 3;  ridgePnl.Layout.Column = [1 2];
            rhgl = uigridlayout(ridgePnl, [2, 10]);
            rhgl.RowHeight   = {30, 30};
            rhgl.ColumnWidth = {150, 'fit', 130, 'fit', 70, 20, 70, 'fit', 70, 100};
            rhgl.Padding     = [5 5 5 5];  rhgl.ColumnSpacing = 4;  rhgl.RowSpacing = 3;
            
            cbShowRidge = uicheckbox(rhgl, 'Text', 'Show ridge overlay', ...
                'Value', false);
            uilabel(rhgl, 'Text', 'Method:', 'HorizontalAlignment', 'right');
            ddMethod = uidropdown(rhgl, ...
                'Items', {'max|dA/dt|', 'max|deltaA|'}, ...
                'Value', app.ridgeMethod);
            uilabel(rhgl, 'Text', sprintf('Search t (%s):', app.tUnitTxt()), 'HorizontalAlignment', 'right');
            edTmin = uieditfield(rhgl, 'numeric', 'Value', app.ridgeTmin);
            uilabel(rhgl, 'Text', 'to', 'HorizontalAlignment', 'center');
            edTmax = uieditfield(rhgl, 'numeric', 'Value', app.ridgeTmax);
            uilabel(rhgl, 'Text', app.uiStr('\lambda-smooth:'), 'HorizontalAlignment', 'right');
            edSm = uieditfield(rhgl, 'numeric', 'Value', app.ridgeSmoothN, ...
                'Limits', [1 200], 'RoundFractionalValues', 'on');
            btnRecRidge = uibutton(rhgl, 'Text', 'Recompute');
            
            uilabel(rhgl, 'Text', '');
            uilabel(rhgl, 'Text', 'Auto-suggest N points:', ...
                'HorizontalAlignment', 'right');
            edAutoN = uieditfield(rhgl, 'numeric', 'Value', 8, ...
                'Limits', [4 100], 'RoundFractionalValues', 'on');
            uilabel(rhgl, 'Text', '');
            btnAuto = uibutton(rhgl, 'Text', 'Auto-place');
            uilabel(rhgl, 'Text', '');
            uilabel(rhgl, 'Text', '');
            uilabel(rhgl, 'Text', '');
            uilabel(rhgl, 'Text', '');
            uilabel(rhgl, 'Text', '');
            
            ctrl = uigridlayout(gl, [1, 11]);
            ctrl.Layout.Row = 4;  ctrl.Layout.Column = [1 2];
            ctrl.ColumnWidth = {'fit', 70, 'fit', 70, 'fit', 'fit', 'fit', 'fit', '1x', 'fit', 'fit'};
            uilabel(ctrl, 'Text', sprintf('View t_min (%s):', app.tUnitTxt()), ...
                'HorizontalAlignment', 'right');
            edViewTmin = uieditfield(ctrl, 'numeric', ...
                'Value', app.chirpViewTMin);
            uilabel(ctrl, 'Text', sprintf('t_max (%s):', app.tUnitTxt()), ...
                'HorizontalAlignment', 'right');
            edViewTmax = uieditfield(ctrl, 'numeric', ...
                'Value', app.chirpViewTMax);
            btnViewReset = uibutton(ctrl, 'Text', 'Reset view');
            btnSavePt = uibutton(ctrl, 'Text', 'Save point', ...
                'BackgroundColor', [0.30 0.55 0.80], 'FontColor', 'w', ...
                'Tooltip', app.uiStr(['Commit the current crosshair position ' ...
                    '(\lambda, t_0) to the points table. Use after ' ...
                    'click-and-drag to position the crosshair on the ' ...
                    'kinetic rising edge.']));
            btnUndo = uibutton(ctrl, 'Text', 'Remove Last');
            btnClr  = uibutton(ctrl, 'Text', 'Clear Points');
            uilabel(ctrl, 'Text', '');
            btnCancel = uibutton(ctrl, 'Text', 'Cancel', ...
                'ButtonPushedFcn', @(~,~) close(f));  %#ok<NASGU>
            btnApply = uibutton(ctrl, 'Text', 'Fit & Apply', ...
                'BackgroundColor', [0.30 0.65 0.35], 'FontColor', 'w');
            
            hPts        = gobjects(0);
            hFitOnMap   = gobjects(0);
            hRidgeOnMap = gobjects(0);
            % Candidate-point state: a single (lambda, t) the user is
            % positioning by click + drag. It becomes a saved point only
            % when "Save point" is pressed.
            candWL      = NaN;
            candT       = NaN;
            hCandV      = [];      % vertical line on the heatmap
            hCandH      = [];      % horizontal line on the heatmap
            % lastClick* remain as the most recent click/drag position,
            % which kin preview also tracks.
            lastClickWL = [];
            lastClickT  = [];
            draggingChirp = false;
            
            drawMapOnce();
            redraw();
            applyViewRange();
            
            cbShowRidge.ValueChangedFcn = @(~,~) toggleRidge();
            ddMethod.ValueChangedFcn    = @(src,~) setAndRecompute('method', src.Value);
            edTmin.ValueChangedFcn      = @(src,~) setAndRecompute('tmin',   src.Value);
            edTmax.ValueChangedFcn      = @(src,~) setAndRecompute('tmax',   src.Value);
            edSm.ValueChangedFcn        = @(src,~) setAndRecompute('smooth', round(src.Value));
            btnRecRidge.ButtonPushedFcn = @(~,~) recomputeRidge();
            btnAuto.ButtonPushedFcn     = @(~,~) autoPlace(round(edAutoN.Value));
            btnUndo.ButtonPushedFcn     = @(~,~) undoPt();
            btnClr.ButtonPushedFcn      = @(~,~) clearPts();
            btnSavePt.ButtonPushedFcn   = @(~,~) savePt();
            applyBusy = false;   % true while the button below is being processed
            btnApply.ButtonPushedFcn    = @(~,~) applyChirp();
            edViewTmin.ValueChangedFcn   = @(~,~) onViewChange();
            edViewTmax.ValueChangedFcn   = @(~,~) onViewChange();
            btnViewReset.ButtonPushedFcn = @(~,~) resetView();
            
            function drawMapOnce()
                base = app.getChirpBaseData();
                cla(axMap);
                wl = app.wavelength(:).';  t = app.delay(:).';
                [W, T] = meshgrid(wl, t);
                C = base.';
                h = surface(axMap, W, T, zeros(size(C)), C, ...
                    'EdgeColor', 'none', 'FaceColor', 'flat');
                view(axMap, 2);
                axMap.YDir   = 'normal';
                axMap.YScale = 'linear';
                cl = max(abs(base(:)), [], 'omitnan');
                if isempty(cl) || cl == 0 || isnan(cl), cl = 1; end
                axMap.CLim = [-cl cl];  colormap(axMap, app.getColormap());
                cb = colorbar(axMap);  cb.Label.String = '\DeltaA or Intensity';
                h.ButtonDownFcn     = @(~,ev) addPt(ev);
                axMap.ButtonDownFcn = @(~,ev) addPt(ev);
                if app.bgApplied
                    title(axMap, 'Click + drag, then Save point  -  BG-corrected base');
                else
                    title(axMap, 'Click + drag, then Save point  -  raw base');
                end
                xlabel(axMap, app.xAxLabel());
                ylabel(axMap, sprintf('Delay time (%s)', app.tUnitAx()));
                % Lock limits explicitly via property assignment.
                axMap.XLim = [min(wl), max(wl)];
                tmin = app.chirpViewTMin;  tmax = app.chirpViewTMax;
                if isempty(tmin) || ~isfinite(tmin), tmin = t(1); end
                if isempty(tmax) || ~isfinite(tmax), tmax = t(end); end
                if tmax <= tmin, tmax = tmin + 0.1; end
                axMap.YLim = [tmin tmax];
                axMap.XLimMode = 'manual';
                axMap.YLimMode = 'manual';
                axMap.NextPlot = 'add';   % permanent "hold on" equivalent
            end
            
            function addPt(ev)
                % Click on the heatmap places (or moves) the candidate
                % crosshair. The candidate is NOT yet a saved point --
                % press "Save point" to commit it to the table.
                pt = ev.IntersectionPoint;
                draggingChirp = true;
                placeCandidate(pt(1), pt(2));
            end
            
            function onChirpDrag()
                if ~draggingChirp, return; end
                if ~app.mouseInsideFigure(f)
                    onChirpDragUp();  return;
                end
                cp = axMap.CurrentPoint;
                if isempty(cp) || size(cp, 1) < 1, return; end
                placeCandidate(cp(1, 1), cp(1, 2));
            end
            
            function onChirpDragUp()
                draggingChirp = false;
            end
            
            function placeCandidate(wl, t)
                % Clamp to data + current view ranges and update the
                % candidate crosshair + kinetics preview.
                wlR = [min(app.wavelength), max(app.wavelength)];
                tR  = [app.chirpViewTMin, app.chirpViewTMax];
                if isempty(tR(1)) || ~isfinite(tR(1)), tR(1) = min(app.delay); end
                if isempty(tR(2)) || ~isfinite(tR(2)), tR(2) = max(app.delay); end
                candWL = min(max(wl, wlR(1)), wlR(2));
                candT  = min(max(t,  tR(1)),  tR(2));
                lastClickWL = candWL;  lastClickT = candT;
                drawCandidate();
                updateKinPreview();
                drawnow limitrate;
            end
            
            function drawCandidate()
                % (Re)draw the candidate crosshair on the heatmap.
                if ~isempty(hCandV) && isgraphics(hCandV), delete(hCandV); end
                if ~isempty(hCandH) && isgraphics(hCandH), delete(hCandH); end
                hCandV = [];  hCandH = [];
                if ~isfinite(candWL) || ~isfinite(candT), return; end
                % Use line() in the same axes -- xline/yline do not
                % stretch across uiaxes when NextPlot='add' / manual
                % limits the way the heatmap is set up here, so we draw
                % plain lines spanning the current axis ranges.
                xl = axMap.XLim;  yl = axMap.YLim;
                hCandH = line(axMap, xl, [candT candT], ...
                    'Color', [0 0 0 0.7], 'LineStyle', '-', ...
                    'LineWidth', 1.0, 'HitTest', 'off', ...
                    'PickableParts', 'none');
                hCandV = line(axMap, [candWL candWL], yl, ...
                    'Color', [0 0 0 0.7], 'LineStyle', '-', ...
                    'LineWidth', 1.0, 'HitTest', 'off', ...
                    'PickableParts', 'none');
            end
            
            function savePt()
                % Commit the current candidate to the points list.
                if ~isfinite(candWL) || ~isfinite(candT)
                    uialert(f, ['No candidate point yet -- click on ' ...
                        'the heatmap first.'], 'No point');
                    return;
                end
                app.chirpPts(end+1, :) = [candWL, candT];
                redraw();
            end
            
            function undoPt()
                if ~isempty(app.chirpPts)
                    app.chirpPts(end, :) = [];
                    redraw();
                end
            end
            
            function clearPts()
                app.chirpPts = zeros(0, 2);
                redraw();
            end
            
            function redraw()
                % Use line() primitives only -- they do NOT interact with hold
                % state, NextPlot, or axis limit modes, so the heatmap beneath
                % them does not re-render or shift.
                delete(hPts(isgraphics(hPts)));        hPts      = gobjects(0);
                delete(hFitOnMap(isgraphics(hFitOnMap))); hFitOnMap = gobjects(0);
                tbl.Data = app.chirpPts;
                
                if ~isempty(app.chirpPts)
                    hPts = line(axMap, app.chirpPts(:,1), app.chirpPts(:,2), ...
                        'LineStyle', 'none', 'Marker', '+', ...
                        'MarkerSize', 14, 'Color', 'w', 'LineWidth', 1.8, ...
                        'HitTest', 'off', 'PickableParts', 'none');
                end
                
                cla(axFit);
                nPts = size(app.chirpPts, 1);
                if nPts >= 4
                    try
                        [pp, rmsFit] = app.fitChirpParams(app.chirpPts);
                        wlFine = linspace(min(app.wavelength), ...
                                          max(app.wavelength), 400);
                        tFine  = app.chirpModel(pp, wlFine);
                        hFitOnMap = line(axMap, wlFine, tFine, ...
                            'Color', 'k', 'LineWidth', 1.6, ...
                            'HitTest', 'off', 'PickableParts', 'none');
                        hold(axFit, 'on');
                        plot(axFit, wlFine, tFine, 'b-', 'LineWidth', 1.4);
                        scatter(axFit, app.chirpPts(:,1), app.chirpPts(:,2), ...
                            50, 'ro', 'filled');
                        hold(axFit, 'off');  grid(axFit, 'on');
                        xlabel(axFit, app.xAxLabel());
                        ylabel(axFit, sprintf('t_0 (%s)', app.tUnitAx()));
                        title(axFit, sprintf(...
                            'N=%d   RMS=%.3g %s\na=%.3g  b=%.3g  c=%.3g  d=%.3g', ...
                            nPts, rmsFit, app.tUnitAx(), ...
                            pp(1), pp(2), pp(3), pp(4)));
                    catch ME
                        title(axFit, sprintf('Fit failed: %s', ME.message));
                    end
                else
                    title(axFit, sprintf('Need %d more point(s)  (>= 4 required)', ...
                        max(0, 4 - nPts)));
                end
                % Candidate crosshair sits on top of everything else.
                drawCandidate();
            end
            
            function updateKinPreview()
                if isempty(lastClickWL), return; end
                base = app.getChirpBaseData();
                [~, iw] = min(abs(app.wavelength - lastClickWL));
                cla(axKinP);  hold(axKinP, 'on');
                plot(axKinP, app.delay, base(iw, :), 'k-', 'LineWidth', 1.2);
                % Vertical line at the current candidate / clicked t,
                % matching the main-window crosshair style (black solid,
                % slightly transparent).
                xline(axKinP, lastClickT, 'k-', 'LineWidth', 1.0, ...
                    'Alpha', 0.7, ...
                    'Label', sprintf('t = %.3g %s', ...
                        lastClickT, app.tUnitAx()), ...
                    'LabelHorizontalAlignment', 'left', ...
                    'LabelVerticalAlignment', 'top');
                if ~isempty(app.ridgeWL)
                    [~, jr] = min(abs(app.ridgeWL - app.wavelength(iw)));
                    xline(axKinP, app.ridgeT(jr), 'g--', 'LineWidth', 1.2, ...
                        'Label', sprintf('t_{ridge} = %.3g %s', ...
                            app.ridgeT(jr), app.tUnitAx()));
                end
                yline(axKinP, 0, 'k:');
                hold(axKinP, 'off');  grid(axKinP, 'on');
                % Match the kinetic preview x-range to the heatmap view
                % window so the user sees the same delay region in both
                % panels.
                tLo = app.chirpViewTMin;  tHi = app.chirpViewTMax;
                if isempty(tLo) || ~isfinite(tLo), tLo = app.delay(1);   end
                if isempty(tHi) || ~isfinite(tHi), tHi = app.delay(end); end
                if tHi <= tLo, tLo = app.delay(1);  tHi = app.delay(end); end
                xlim(axKinP, [tLo, tHi]);
                title(axKinP, sprintf('%s = %.1f %s', ...
                    app.xSym(), app.wavelength(iw), app.xUnitAx()));
            end
            
            function setAndRecompute(field, val)
                switch field
                    case 'method', app.ridgeMethod  = val;
                    case 'tmin',   app.ridgeTmin    = val;
                    case 'tmax',   app.ridgeTmax    = val;
                    case 'smooth', app.ridgeSmoothN = val;
                end
                if cbShowRidge.Value, recomputeRidge(); end
            end
            
            function recomputeRidge()
                app.computeRidge(app.ridgeMethod, ...
                    app.ridgeTmin, app.ridgeTmax, app.ridgeSmoothN);
                delete(hRidgeOnMap(isgraphics(hRidgeOnMap)));
                hRidgeOnMap = gobjects(0);
                if cbShowRidge.Value && ~isempty(app.ridgeT)
                    hRidgeOnMap = line(axMap, app.ridgeWL, app.ridgeT, ...
                        'Color', [0 0.7 0], 'LineWidth', 1.2, ...
                        'HitTest', 'off', 'PickableParts', 'none');
                end
            end
            
            function toggleRidge()
                if cbShowRidge.Value
                    recomputeRidge();
                else
                    delete(hRidgeOnMap(isgraphics(hRidgeOnMap)));
                    hRidgeOnMap = gobjects(0);
                end
            end
            
            function autoPlace(N)
                if isempty(app.ridgeT)
                    app.computeRidge(app.ridgeMethod, ...
                        app.ridgeTmin, app.ridgeTmax, app.ridgeSmoothN);
                end
                if isempty(app.ridgeT)
                    uialert(f, 'Could not compute ridge (check search window).', ...
                        'No ridge');  return;
                end
                wls = linspace(app.wavelength(1), app.wavelength(end), N+2);
                wls = wls(2:end-1);
                newPts = zeros(N, 2);
                for k = 1:N
                    [~, j] = min(abs(app.ridgeWL - wls(k)));
                    newPts(k, :) = [app.ridgeWL(j), app.ridgeT(j)];
                end
                app.chirpPts = newPts;
                redraw();
                if ~cbShowRidge.Value
                    cbShowRidge.Value = true;
                    recomputeRidge();
                end
            end
            
            function applyViewRange()
                tmin = app.chirpViewTMin;
                tmax = app.chirpViewTMax;
                if isempty(tmin) || ~isfinite(tmin), tmin = app.delay(1);   end
                if isempty(tmax) || ~isfinite(tmax), tmax = app.delay(end); end
                if tmax <= tmin
                    tmax = tmin + eps(tmin) * 10;
                    if tmax <= tmin, tmax = tmin + 0.1; end
                end
                % Only update when the limits differ, to avoid unnecessary
                % re-rendering that can cause visible flicker/shift.
                if ~strcmp(axMap.YScale, 'linear')
                    axMap.YScale = 'linear';
                end
                cur = axMap.YLim;
                if abs(cur(1) - tmin) > eps(max(1, abs(tmin))) || ...
                   abs(cur(2) - tmax) > eps(max(1, abs(tmax)))
                    axMap.YLim = [tmin tmax];
                end
                axMap.YLimMode = 'manual';
                axMap.XLimMode = 'manual';
            end
            
            function onViewChange()
                app.chirpViewTMin = edViewTmin.Value;
                app.chirpViewTMax = edViewTmax.Value;
                applyViewRange();
                if ~isempty(lastClickWL), updateKinPreview(); end
            end
            
            function resetView()
                app.chirpViewTMin = app.delay(1);
                app.chirpViewTMax = min(app.defaultChirpTMax(), app.delay(end));
                edViewTmin.Value = app.chirpViewTMin;
                edViewTmax.Value = app.chirpViewTMax;
                applyViewRange();
                if ~isempty(lastClickWL), updateKinPreview(); end
            end
            
            function applyChirp()
                % A second click while the main window is being redrawn
                % is ignored (see doApply of the crop window).
                if applyBusy, return; end
                if size(app.chirpPts, 1) < 4
                    uialert(f, 'At least 4 points are required.', ...
                        'Not enough points');  return;
                end
                try
                    [pp, rmsFit] = app.fitChirpParams(app.chirpPts);
                catch ME
                    uialert(f, sprintf('Fit failed:\n%s', ME.message), ...
                        'Fit error');  return;
                end
                applyBusy = true;
                try
                    app.chirpParams  = pp;
                    app.chirpFitRMS  = rmsFit;
                    app.chirpApplied = true;
                    app.updateAll();
                catch errA
                    applyBusy = false;
                    rethrow(errA);
                end
                applyBusy = false;
                if isvalid(f), close(f); end
            end
        end
        
        % =================================================================
        % Global analysis window
        % =================================================================
        function openGlobalAnalysisWindow(app)
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');  return;
            end
            if ~isempty(app.gaFig) && isvalid(app.gaFig)
                figure(app.gaFig);  return;
            end
            
            f = uifigure('Name', 'Global Analysis', ...
                'Position', app.fitFigurePosition(1420, 1040));
            app.applyLightTheme(f);
            app.gaFig = f;
            f.WindowButtonMotionFcn = @(~,~) onGaMotion();
            f.WindowButtonUpFcn     = @(~,~) onGaUp();
            
            gl = uigridlayout(f, [1, 2]);
            gl.ColumnWidth = {340, '1x'};
            gl.Padding = [10 10 10 10];  gl.ColumnSpacing = 10;
            
            % ---------- Left panel: parameter setup ----------
            left = uipanel(gl, 'Title', 'Fit setup');
            left.Layout.Column = 1;
            lgl = uigridlayout(left, [16, 1]);
            lgl.RowHeight = {30, 22, 150, 30, 30, 30, 30, 30, 30, 32, 32, 30, 30, 22, 22, '1x'};
            lgl.Padding = [8 8 8 8];  lgl.RowSpacing = 6;
            
            row1 = uigridlayout(lgl, [1, 4]);
            row1.ColumnWidth = {'fit', 60, 'fit', 'fit'};
            row1.Padding = [0 0 0 0];  row1.ColumnSpacing = 6;
            uilabel(row1, 'Text', 'Components:');
            ddN = uidropdown(row1, 'Items', {'1','2','3','4','5'}, ...
                'Value', num2str(app.gaNComp));
            cbInf = uicheckbox(row1, 'Text', app.uiStr('Include \tau = \infty'), ...
                'Value', app.gaHasInf);
            uibutton(row1, 'Text', 'SVD?', ...
                'Tooltip', ['Estimate the number of components from the ' ...
                    'singular values of the data (count above the noise ' ...
                    'floor) and set Components accordingly. The ratios ' ...
                    'are printed in the result box; the choice remains ' ...
                    'yours.'], ...
                'ButtonPushedFcn', @(~,~) onSuggestN());
            
            uilabel(lgl, 'Text', sprintf('Initial time constants (%s, positive):', app.tUnitTxt()), ...
                'FontAngle', 'italic', 'FontColor', [0.3 0.3 0.3]);
            
            tauTable = uitable(lgl, ...
                'ColumnName',   app.uiStr({ ...
                                 sprintf('\\tau_{init} (%s)', app.tUnitTxt()), ...
                                 '\tau fixed', '\beta_{init}', 'Stretched'}), ...
                'ColumnFormat', {'numeric', 'logical', 'numeric', 'logical'}, ...
                'ColumnEditable', [true true true true], ...
                'ColumnWidth',  {85, 70, 70, 70}, ...
                'RowName', 'numbered', ...
                'CellEditCallback', @(~,~) refreshSeqRowState());
            
            % IRF t0
            rowT0 = uigridlayout(lgl, [1, 3]);
            rowT0.ColumnWidth = {'fit', 80, 'fit'};
            rowT0.Padding = [0 0 0 0];  rowT0.ColumnSpacing = 6;
            uilabel(rowT0, 'Text', sprintf('IRF t_0 (%s):', app.tUnitTxt()));
            edT0 = uieditfield(rowT0, 'numeric', 'Value', app.gaT0);
            cbT0fix = uicheckbox(rowT0, 'Text', 'Fixed', ...
                'Value', app.gaT0Fixed);
            
            % IRF FWHM
            rowFW = uigridlayout(lgl, [1, 3]);
            rowFW.ColumnWidth = {'fit', 80, 'fit'};
            rowFW.Padding = [0 0 0 0];  rowFW.ColumnSpacing = 6;
            uilabel(rowFW, 'Text', sprintf('IRF FWHM (%s):', app.tUnitTxt()));
            edFW = uieditfield(rowFW, 'numeric', 'Value', app.gaFWHM, ...
                'Limits', [1e-6 Inf]);
            cbFWfix = uicheckbox(rowFW, 'Text', 'Fixed', ...
                'Value', app.gaFWHMFixed);
            
            % IRF model: Gaussian IRF included ("Numerical conv.") or
            % none (pure exponentials, for data cropped to start after
            % the IRF). A value left by an older version ('skip') is
            % mapped to the IRF model first: the dropdown would reject
            % it.
            app.gaIRFTreat = app.irfModeNorm(app.gaIRFTreat);
            rowIRF = uigridlayout(lgl, [1, 3]);
            rowIRF.ColumnWidth = {'fit', 165, '1x'};
            rowIRF.Padding = [0 0 0 0];  rowIRF.ColumnSpacing = 6;
            uilabel(rowIRF, 'Text', 'IRF model:');
            ddIRF = uidropdown(rowIRF, ...
                'Items',     {'None (pure exponentials)', ...
                              'Numerical conv.'}, ...
                'ItemsData', {'none', 'numerical'}, ...
                'Value',     app.gaIRFTreat, ...
                'Tooltip', ['Numerical conv.: the model is convolved ' ...
                    'with the Gaussian IRF (t0, FWHM). Standard exp ' ...
                    'components use the closed-form erfcx expression, ' ...
                    'stretched components a numerical quadrature at ' ...
                    'each delay.' newline ...
                    'None: no IRF convolution -- fit a plain ' ...
                    'sum of exponentials exp(-(t-t0)/tau)*H(t-t0). Use ' ...
                    'after cropping the data to start beyond the IRF; ' ...
                    't0 and FWHM are excluded from the fit (t0 acts as ' ...
                    'a fixed onset).'], ...
                'ValueChangedFcn', @(~,~) onIRFModeChange());
            uilabel(rowIRF, 'Text', '');
            
            % Coherent artifact: extra terms at time zero with free
            % spectra -- a Gaussian and its first and second time
            % derivative. Their width is that of the IRF, or one of
            % their own (fitted or fixed).
            app.gaCA = app.caOptions(app.gaCA);
            rowCA = uigridlayout(lgl, [1, 2]);
            rowCA.ColumnWidth = {'fit', '1x'};
            rowCA.Padding = [0 0 0 0];  rowCA.ColumnSpacing = 6;
            uilabel(rowCA, 'Text', 'Coh. artifact:');
            ddCA = uidropdown(rowCA, ...
                'Items',     {'None', '1 term: Gaussian', ...
                              '2 terms: + 1st derivative', ...
                              '3 terms: + 2nd derivative'}, ...
                'ItemsData', {0, 1, 2, 3}, ...
                'Value',     app.gaCA.order, ...
                'Tooltip', ['Signal that exists only while pump and ' ...
                    'probe overlap (cross-phase modulation, two-photon ' ...
                    'absorption, stimulated Raman of the solvent). It is ' ...
                    'fitted with up to three extra terms at time zero, ' ...
                    'each with a free spectrum: a Gaussian and its first ' ...
                    'and second time derivative.' newline ...
                    'Without these terms an artifact distorts the fitted ' ...
                    'IRF and can appear as a very fast component. With ' ...
                    'them, time constants below the IRF width are less ' ...
                    'well determined: switch them on only when the data ' ...
                    'show an artifact.' newline ...
                    'The artifact spectra are drawn dashed in the DADS ' ...
                    'plot and are not part of the DADS / EADS.'], ...
                'ValueChangedFcn', @(~,~) syncCAEnable());
            caW0 = app.gaCA.fwhm;
            if ~(isfinite(caW0) && caW0 > 0), caW0 = app.gaFWHM; end
            rowCAW = uigridlayout(lgl, [1, 3]);
            rowCAW.ColumnWidth = {'fit', '1x', 70};
            rowCAW.Padding = [0 0 0 0];  rowCAW.ColumnSpacing = 6;
            lblCAW = uilabel(rowCAW, 'Text', 'Artifact width:');
            ddCAW = uidropdown(rowCAW, ...
                'Items',     {'Same as IRF', 'Own, fitted', 'Own, fixed'}, ...
                'ItemsData', {'irf', 'free', 'fixed'}, ...
                'Value',     app.gaCA.mode, ...
                'Tooltip', ['Width (FWHM) of the artifact terms.' newline ...
                    'Same as IRF: the IRF FWHM, so the artifact helps ' ...
                    'to determine the IRF.' newline ...
                    'Own: a separate width, fitted (kept within a ' ...
                    'factor of 5 of the IRF FWHM) or fixed at the value ' ...
                    'on the right. The IRF FWHM itself is then hardly ' ...
                    'determined by the fit: keep it fixed.'], ...
                'ValueChangedFcn', @(~,~) syncCAEnable());
            edCAW = uieditfield(rowCAW, 'numeric', 'Value', caW0, ...
                'Limits', [1e-6 Inf], ...
                'Tooltip', sprintf(['Own width of the artifact terms ' ...
                    '(FWHM, %s): starting value or fixed value.'], ...
                    app.tUnitTxt()));
            
            % Fit time range: delay points outside [t_min, t_max] are
            % excluded from the fit (amplitudes AND loss). Lets you fit
            % e.g. only t >= 300 fs without cropping the dataset.
            rowFitT = uigridlayout(lgl, [1, 5]);
            rowFitT.ColumnWidth = {'fit', 62, 'fit', 62, 'fit'};
            rowFitT.Padding = [0 0 0 0];  rowFitT.ColumnSpacing = 4;
            uilabel(rowFitT, 'Text', sprintf('Fit t range (%s):', ...
                app.tUnitTxt()));
            edFitTmin = uieditfield(rowFitT, 'numeric', ...
                'Value', min(app.delay), ...
                'Tooltip', 'Lower bound of the fitted delay range.');
            uilabel(rowFitT, 'Text', '~', ...
                'HorizontalAlignment', 'center');
            edFitTmax = uieditfield(rowFitT, 'numeric', ...
                'Value', max(app.delay), ...
                'Tooltip', 'Upper bound of the fitted delay range.');
            uibutton(rowFitT, 'Text', 'Full', ...
                'Tooltip', 'Reset to the full delay range.', ...
                'ButtonPushedFcn', @(~,~) setFitTFull());
            
            % Run Fit / Reset
            rowBtn = uigridlayout(lgl, [1, 2]);
            rowBtn.ColumnWidth = {'1x', '1x'};
            rowBtn.Padding = [0 0 0 0];  rowBtn.ColumnSpacing = 5;
            btnRun = uibutton(rowBtn, 'Text', 'Run Fit', ...
                'BackgroundColor', [0.30 0.65 0.35], 'FontColor', 'w');
            btnReset = uibutton(rowBtn, 'Text', 'Reset');
            
            % Export buttons row -- one per target
            rowExp = uigridlayout(lgl, [1, 7]);
            rowExp.ColumnWidth = {'fit', '1x', '1x', '1x', '1x', '1x', '1x'};
            rowExp.Padding = [0 0 0 0];  rowExp.ColumnSpacing = 3;
            uilabel(rowExp, 'Text', 'Export:', ...
                'HorizontalAlignment', 'right');
            btnExpDADS  = uibutton(rowExp, 'Text', 'DADS');
            btnExpEADS = uibutton(rowExp, 'Text', 'EADS');
            btnExpFit  = uibutton(rowExp, 'Text', 'Fit');
            btnExpRes  = uibutton(rowExp, 'Text', 'Resid.');
            btnExpKin  = uibutton(rowExp, 'Text', app.uiStr('Kin \lambda'));
            btnExpPop  = uibutton(rowExp, 'Text', 'Popul.', ...
                'Tooltip', ['Export the sequential-model population ' ...
                    'kinetics c_i(t) shown in the bottom-right panel.']);
            
            % Reading of a stretched step in the sequential model. It
            % concerns the SAS and the populations only (the fit does
            % not depend on it), so a change is applied to the stored
            % result at once.
            if ~any(strcmp(app.gaSeqModel, {'age', 't0'}))
                app.gaSeqModel = 'age';
            end
            rowSeq = uigridlayout(lgl, [1, 2]);
            rowSeq.ColumnWidth = {'fit', '1x'};
            rowSeq.Padding = [0 0 0 0];  rowSeq.ColumnSpacing = 6;
            seqTip = ['How a stretched step is read in the sequential ' ...
                'model S1 -> S2 -> ... (SAS and populations; the fit ' ...
                'itself is not affected, no new fit is needed).' newline ...
                'Clock from entry: the stretched decay of a state ' ...
                'starts when the population enters that state (a ' ...
                'distribution of lifetimes of the state).' newline ...
                'Clock from excitation: all steps share one clock ' ...
                'started by the pump at t0 -- a time-dependent rate ' ...
                'k(t) = (beta/tau)*(t/tau)^(beta-1), as in dispersive ' ...
                'or diffusion-limited kinetics.' newline ...
                'The two differ only when a component other than the ' ...
                'fastest one is stretched.'];
            lblSeq = uilabel(rowSeq, 'Text', 'Stretched step:', ...
                'Tooltip', seqTip);
            ddSeq = uidropdown(rowSeq, ...
                'Items',     {'Clock from entry into state', ...
                              'Clock from excitation (t0)'}, ...
                'ItemsData', {'age', 't0'}, ...
                'Value',     app.gaSeqModel, ...
                'Tooltip',   seqTip, ...
                'ValueChangedFcn', @(~,~) onSeqModelChange());
            
            % Map delay-scale controls (2D maps in right panel)
            rowMap = uigridlayout(lgl, [1, 4]);
            rowMap.ColumnWidth = {'fit', 90, 'fit', 55};
            rowMap.Padding = [0 0 0 0];  rowMap.ColumnSpacing = 6;
            uilabel(rowMap, 'Text', 'Map scale:', ...
                'HorizontalAlignment', 'right');
            ddMapScale = uidropdown(rowMap, ...
                'Items',     {'Linear', 'Log', 'Split'}, ...
                'ItemsData', {'linear', 'log', 'split'}, ...
                'Value', 'log');
            uilabel(rowMap, 'Text', sprintf('Thr. (%s):', ...
                app.tUnitTxt()), 'HorizontalAlignment', 'right');
            edMapThr = uieditfield(rowMap, 'numeric', ...
                'Value', app.splitThreshold, 'Limits', [eps Inf]);
            
            % Map delay-axis range (applies to linear axis: linear mode
            % full, split mode lower row).
            rowMapT = uigridlayout(lgl, [1, 4]);
            rowMapT.ColumnWidth = {'fit', 70, 'fit', 70};
            rowMapT.Padding = [0 0 0 0];  rowMapT.ColumnSpacing = 6;
            uilabel(rowMapT, 'Text', sprintf('t min (%s):', ...
                app.tUnitTxt()), 'HorizontalAlignment', 'right', ...
                'Tooltip', ['Lower bound of the delay axis on the 2D ' ...
                'maps. Applies to linear mode (whole axis) and split ' ...
                'mode (linear lower panel only).']);
            edMapTmin = uieditfield(rowMapT, 'numeric', 'Value', -1);
            uilabel(rowMapT, 'Text', sprintf('t max (%s):', ...
                app.tUnitTxt()), 'HorizontalAlignment', 'right', ...
                'Tooltip', ['Upper bound of the delay axis on the 2D ' ...
                'maps. Default: max delay in the data.']);
            if isempty(app.delay)
                tmaxDef = 1000;
            else
                tmaxDef = max(app.delay);
            end
            edMapTmax = uieditfield(rowMapT, 'numeric', 'Value', tmaxDef);
            
            % Status + results
            lblStatus = uilabel(lgl, 'Text', 'Ready.', ...
                'FontColor', [0.2 0.2 0.2]);
            
            txtResult = uitextarea(lgl, 'Editable', 'off', ...
                'Value', {'Results will appear here after fitting.'}, ...
                'FontName', 'Monospaced');
            
            % ---------- Right panel: all results in one view ----------
            right = uipanel(gl, 'Title', 'Results');
            right.Layout.Column = 2;
            rgl = uigridlayout(right, [5, 1]);
            rgl.RowHeight   = {'1x', '0.85x', '0.85x', 30, '0.95x'};
            rgl.Padding     = [6 6 6 6];
            rgl.RowSpacing  = 6;
            
            % --- Row 1: 2D maps strip (experimental / fit / residual) ---
            mapsRow = uigridlayout(rgl, [1, 3]);
            mapsRow.Layout.Row = 1;
            mapsRow.ColumnWidth = {'1x', '1x', '1x'};
            mapsRow.Padding = [0 0 0 0];
            mapsRow.ColumnSpacing = 8;
            
            % Each map is hosted inside its own uipanel so we can rebuild
            % the inner axes layout when the user switches between
            % Linear/Log (1 axis) and Split (2 stacked axes).
            mapPanelExp = uipanel(mapsRow, 'BorderType', 'none', ...
                'Title', app.uiStr('Experimental  (click to pick \lambda)'));
            mapPanelFit = uipanel(mapsRow, 'BorderType', 'none', ...
                'Title', 'Global fit');
            mapPanelRes = uipanel(mapsRow, 'BorderType', 'none', ...
                'Title', 'Residual');
            
            % Axes handles, filled in by buildMapAxes(). Each field is
            % either a 1x1 (Linear/Log) or 1x2 (Split: [top log, bot lin])
            % array of uiaxes.  The cbar* companions hold dedicated
            % colorbar hosts in the same grid cells so all three data
            % columns keep the same horizontal extent.
            mapAxExp  = gobjects(0);   cbarAxExp = gobjects(0);
            mapAxFit  = gobjects(0);   cbarAxFit = gobjects(0); %#ok<NASGU>
            mapAxRes  = gobjects(0);   cbarAxRes = gobjects(0);
            
            % Vertical crosshair handles, one per axes. Rebuilt alongside
            % axes in buildMapAxes().
            hMapVLine = gobjects(0);
            
            % --- Row 2: DADS (raw) and DADS (peak-normalized) ---
            dadsRow = uigridlayout(rgl, [1, 2]);
            dadsRow.Layout.Row = 2;
            dadsRow.ColumnWidth = {'1x', '1x'};
            dadsRow.Padding = [0 0 0 0];
            dadsRow.ColumnSpacing = 8;
            
            axDADS = uiaxes(dadsRow);
            xlabel(axDADS, app.xAxLabel());
            ylabel(axDADS, 'Amplitude (\DeltaA or Intensity)');
            title(axDADS, 'Decay-associated difference spectra');
            grid(axDADS, 'on');
            
            axDADSnorm = uiaxes(dadsRow);
            xlabel(axDADSnorm, app.xAxLabel());
            ylabel(axDADSnorm, 'Normalized amplitude');
            title(axDADSnorm, 'DADS (peak-normalized)');
            grid(axDADSnorm, 'on');
            
            % --- Row 3: EADS (raw) and EADS (peak-normalized) ---
            eadsRow = uigridlayout(rgl, [1, 2]);
            eadsRow.Layout.Row = 3;
            eadsRow.ColumnWidth = {'1x', '1x'};
            eadsRow.Padding = [0 0 0 0];
            eadsRow.ColumnSpacing = 8;
            
            axEADS = uiaxes(eadsRow);
            xlabel(axEADS, app.xAxLabel());
            ylabel(axEADS, 'Amplitude (\DeltaA or Intensity)');
            title(axEADS, ['Evolution-associated difference spectra  ' ...
                '(sequential model)']);
            grid(axEADS, 'on');
            
            axEADSnorm = uiaxes(eadsRow);
            xlabel(axEADSnorm, app.xAxLabel());
            ylabel(axEADSnorm, 'Normalized amplitude');
            title(axEADSnorm, 'EADS (peak-normalized)');
            grid(axEADSnorm, 'on');
            
            % --- Row 4: component visibility checkboxes ---
            cbHost = uipanel(rgl, 'BorderType', 'none');
            cbHost.Layout.Row = 4;
            cbComps = gobjects(0);
            cbCAshow = [];      % "artifact" checkbox of the same row
            
            % --- Row 5: kinetics controls + trace ---
            kinBlock = uigridlayout(rgl, [2, 1]);
            kinBlock.Layout.Row = 5;
            kinBlock.RowHeight = {32, '1x'};
            kinBlock.Padding = [0 0 0 0];  kinBlock.RowSpacing = 3;
            
            kinCtrl = uigridlayout(kinBlock, [1, 7]);
            kinCtrl.Layout.Row = 1;
            kinCtrl.ColumnWidth = {'fit', 90, 'fit', 'fit', 80, 'fit', '1x'};
            kinCtrl.Padding = [0 0 0 0];  kinCtrl.ColumnSpacing = 8;
            uilabel(kinCtrl, 'Text', sprintf('%s (%s):', app.xSymTxt(), app.xUnitTxt()), ...
                'HorizontalAlignment', 'right');
            edKinWL = uieditfield(kinCtrl, 'numeric', 'Value', app.selWL, ...
                'ValueDisplayFormat', '%.2f');
            uibutton(kinCtrl, 'Text', 'Use main selection', ...
                'ButtonPushedFcn', @(~,~) syncKinWLFromMain());
            uilabel(kinCtrl, 'Text', 'Scale:', ...
                'HorizontalAlignment', 'right');
            ddKinScale = uidropdown(kinCtrl, ...
                'Items', {'Log', 'Linear'}, 'Value', 'Log');
            lblPickHint = uilabel(kinCtrl, ...
                'Text', '   (or click on the experimental map above)', ...
                'FontColor', [0.35 0.35 0.35], 'FontAngle', 'italic');  %#ok<NASGU>
            uilabel(kinCtrl, 'Text', '');
            
            % Side-by-side layout: left ~60% data+fit, right ~40% residual.
            kinAxRow = uigridlayout(kinBlock, [1, 3]);
            kinAxRow.Layout.Row = 2;
            kinAxRow.ColumnWidth   = {'3x', '2x', '2x'};
            kinAxRow.Padding       = [0 0 0 0];
            kinAxRow.ColumnSpacing = 8;
            
            axKinFit = uiaxes(kinAxRow);
            axKinFit.Layout.Column = 1;
            xlabel(axKinFit, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axKinFit, '\DeltaA or Intensity');
            title(axKinFit, 'Kinetics fit at selected \lambda');
            grid(axKinFit, 'on');
            
            axKinRes = uiaxes(kinAxRow);
            axKinRes.Layout.Column = 2;
            xlabel(axKinRes, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axKinRes, 'Residual (\DeltaA or Intensity)');
            title(axKinRes, 'Fit residual');
            grid(axKinRes, 'on');
            
            axPop = uiaxes(kinAxRow);
            axPop.Layout.Column = 3;
            xlabel(axPop, sprintf('Delay time (%s)', app.tUnitAx()));
            ylabel(axPop, 'Population');
            title(axPop, 'EADS populations c_i(t)');
            grid(axPop, 'on');
            
            % Local state for nested callbacks
            curKinWL = app.selWL;
            draggingG = false;     % GA map drag state
            dragGaAx  = [];        % axes the drag started on
            
            edKinWL.ValueChangedFcn    = @(src,~) onKinWLChange(src.Value);
            ddKinScale.ValueChangedFcn = @(~,~) onKinScaleChange();
            
            % Initialize table
            fillTable();
            refreshSeqRowState();
            updateEadsPanelTitle();
            ddN.ValueChangedFcn    = @(~,~) onNCompChange();
            btnRun.ButtonPushedFcn    = @(~,~) doRun();
            btnReset.ButtonPushedFcn  = @(~,~) doReset();
            btnExpDADS.ButtonPushedFcn  = @(~,~) exportDADS();
            btnExpEADS.ButtonPushedFcn = @(~,~) exportEADS();
            btnExpFit.ButtonPushedFcn  = @(~,~) exportMatrix2D( ...
                'fit', 'Global fit', 'fit');
            btnExpRes.ButtonPushedFcn  = @(~,~) exportMatrix2D( ...
                'residual', 'Residual', 'residual');
            btnExpKin.ButtonPushedFcn  = @(~,~) exportKinetics();
            btnExpPop.ButtonPushedFcn  = @(~,~) exportPopulation();
            % Apply the current IRF-model mode to the t0 / FWHM controls
            onIRFModeChange();
            ddMapScale.ValueChangedFcn = @(~,~) onMapScaleChange();
            edMapThr.ValueChangedFcn   = @(~,~) onMapScaleChange();
            edMapTmin.ValueChangedFcn  = @(~,~) onMapScaleChange();
            edMapTmax.ValueChangedFcn  = @(~,~) onMapScaleChange();
            
            % Restore a previously used fit time range (session or
            % earlier run in this MATLAB session).
            [edFitTmin.Value, edFitTmax.Value] = app.clampFitRange( ...
                app.gaFitTMin, app.gaFitTMax);
            
            % ---- Restore previous results (if any) ----
            % All fit results live in app.gaResult* properties, which
            % survive the window being closed. If a fit was done in a
            % previous session of this window (or restored from a
            % saved session file), re-render it now so closing /
            % reopening the window does not lose the analysis.
            gaOK = ~isempty(app.gaResultDADS) && ...
                size(app.gaResultDADS, 1) == numel(app.wavelength) && ...
                (isempty(app.gaResultFit) || ...
                 isequal(size(app.gaResultFit), [numel(app.wavelength), numel(app.delay)]));
            if ~isempty(app.gaResultDADS) && ~gaOK
                % Stale result from a different axis (crop / reload
                % without passing through invalidateChildResults).
                app.clearGaResults();
            end
            if gaOK
                lblStatus.Text = sprintf( ...
                    'Restored previous fit. RMS = %.3g', app.gaResultRMS);
                txtResult.Value = {sprintf( ...
                    ['Previous fit restored (RMS = %.3g).  ' ...
                     'Run again to refit.'], app.gaResultRMS)};
                plotResults();
                updateEadsPanelTitle();
            end
            refreshSeqRowState();
            
            % ---------- Nested callbacks ----------
            function fillTable()
                N = app.gaNComp;
                if numel(app.gaTauInit)  ~= N
                    if numel(app.gaTauInit) < N
                        extra = N - numel(app.gaTauInit);
                        if ~isempty(app.gaTauInit)
                            lastTau = app.gaTauInit(end);
                        else
                            lastTau = 1;
                        end
                        app.gaTauInit  = [app.gaTauInit(:).',  ...
                            lastTau * 10.^(1:extra)];
                        app.gaTauFixed = [app.gaTauFixed(:).', false(1, extra)];
                    else
                        app.gaTauInit  = app.gaTauInit(1:N);
                        app.gaTauFixed = app.gaTauFixed(1:N);
                    end
                end
                if numel(app.gaBetaInit) ~= N
                    if numel(app.gaBetaInit) < N
                        app.gaBetaInit = [app.gaBetaInit(:).', ...
                            ones(1, N - numel(app.gaBetaInit))];
                    else
                        app.gaBetaInit = app.gaBetaInit(1:N);
                    end
                end
                if numel(app.gaStretchOn) ~= N
                    if numel(app.gaStretchOn) < N
                        app.gaStretchOn = [app.gaStretchOn(:).', ...
                            false(1, N - numel(app.gaStretchOn))];
                    else
                        app.gaStretchOn = app.gaStretchOn(1:N);
                    end
                end
                tauTable.Data = [num2cell(app.gaTauInit(:)), ...
                                 num2cell(logical(app.gaTauFixed(:))), ...
                                 num2cell(app.gaBetaInit(:)), ...
                                 num2cell(logical(app.gaStretchOn(:)))];
            end
            
            function onSuggestN()
                [nSug, ~, txt] = app.suggestNComponents();
                if isempty(txt), return; end
                nSug = max(1, min(nSug, 5));     % dropdown range
                ddN.Value = num2str(nSug);
                onNCompChange();
                lblStatus.Text = sprintf('SVD suggests N = %d', nSug);
                txtResult.Value = {sprintf('SVD rank estimate: N = %d', nSug), ...
                    txt, '', ['Components set to this value; adjust ' ...
                    'and press Run Fit.']};
            end

            function onNCompChange()
                app.gaNComp = str2double(ddN.Value);
                fillTable();
                refreshSeqRowState();
            end
            
            function refreshSeqRowState()
                % The "Stretched step" selector only matters with
                % stretched components -- in the table (the next fit) or
                % in the stored result (what is shown now).
                anyStretch = false;
                for ii = 1:size(tauTable.Data, 1)
                    if logical(tauTable.Data{ii, 4})
                        anyStretch = true;
                        break;
                    end
                end
                if ~anyStretch && ~isempty(app.gaResultTau)
                    anyStretch = any(resStretch());
                end
                if anyStretch
                    ddSeq.Enable     = 'on';
                    lblSeq.FontColor = [0.10 0.10 0.10];
                else
                    ddSeq.Enable     = 'off';
                    lblSeq.FontColor = [0.55 0.55 0.55];
                end
            end
            
            function readTable()
                if isempty(tauTable.Data), return; end
                D = tauTable.Data;
                N = size(D, 1);
                tau = zeros(1, N);  fix = false(1, N);
                bet = zeros(1, N);  son = false(1, N);
                for ii = 1:N
                    tau(ii) = double(D{ii, 1});
                    fix(ii) = logical(D{ii, 2});
                    bet(ii) = double(D{ii, 3});
                    son(ii) = logical(D{ii, 4});
                end
                app.gaTauInit   = tau;
                app.gaTauFixed  = fix;
                app.gaBetaInit  = bet;
                app.gaStretchOn = son;
            end
            
            function onMapScaleChange()
                % Rebuild the three 2D maps with the current scale choice.
                % plotMaps() internally calls buildMapAxes() which destroys
                % and re-creates the inner axes when the mode changes.
                % The bottom kinetics / population panels follow the map
                % time-scale choice (Split maps use a log kinetics axis).
                if strcmp(ddMapScale.Value, 'linear')
                    ddKinScale.Value = 'Linear';
                else
                    ddKinScale.Value = 'Log';
                end
                if isempty(app.gaResultFit), return; end
                plotMaps();
                onKinScaleChange();
            end
            
            function doReset()
                app.gaNComp    = 3;
                app.gaTauInit  = [1 10 100];
                app.gaTauFixed = [false false false];
                app.gaBetaInit  = [1 1 1];
                app.gaStretchOn = [false false false];
                app.gaIRFTreat  = 'numerical';
                app.gaHasInf   = false;
                app.gaT0       = 0;       app.gaT0Fixed   = true;
                app.gaFWHM     = 0.15;    app.gaFWHMFixed = true;
                ddN.Value = num2str(app.gaNComp);
                cbInf.Value    = app.gaHasInf;
                edT0.Value     = app.gaT0;    cbT0fix.Value = app.gaT0Fixed;
                edFW.Value     = app.gaFWHM;  cbFWfix.Value = app.gaFWHMFixed;
                ddIRF.Value    = app.gaIRFTreat;
                app.gaCA       = app.caOptions([]);
                ddCA.Value     = 0;
                ddCAW.Value    = 'irf';
                edCAW.Value    = app.gaFWHM;
                if ~strcmp(ddSeq.Value, 'age')
                    ddSeq.Value = 'age';
                    onSeqModelChange();
                end
                app.gaFitTMin  = [];
                app.gaFitTMax  = [];
                setFitTFull();
                fillTable();
                refreshSeqRowState();
                onIRFModeChange();
                updateEadsPanelTitle();
                cla(axPop);
                title(axPop, 'EADS populations c_i(t)');
                lblStatus.Text = 'Reset to defaults.';
            end
            
            function doRun()
                % Collect params
                readTable();
                app.gaHasInf    = cbInf.Value;
                app.gaT0        = edT0.Value;
                app.gaT0Fixed   = cbT0fix.Value;
                app.gaFWHM      = edFW.Value;
                app.gaFWHMFixed = cbFWfix.Value;
                app.gaIRFTreat  = ddIRF.Value;
                app.gaSeqModel  = ddSeq.Value;
                app.gaCA = struct('order', ddCA.Value, ...
                    'mode', ddCAW.Value, 'fwhm', edCAW.Value);
                if strcmp(app.gaIRFTreat, 'none')
                    % t0 degenerate with amplitudes; FWHM absent.
                    app.gaT0Fixed   = true;
                    app.gaFWHMFixed = true;
                end
                [app.gaFitTMin, app.gaFitTMax] = app.fitRangeToStore( ...
                    edFitTmin.Value, edFitTmax.Value);
                tSelChk = app.delay >= min([edFitTmin.Value, edFitTmax.Value]) & ...
                          app.delay <= max([edFitTmin.Value, edFitTmax.Value]);
                if nnz(tSelChk) < (app.gaNComp + 2)
                    uialert(f, ['Fit t range contains too few delay ' ...
                        'points. Widen the range.'], 'Invalid range');
                    return;
                end
                
                if any(app.gaTauInit <= 0)
                    uialert(f, 'All tau values must be positive.', ...
                        'Invalid input');  return;
                end
                if any(app.gaBetaInit <= 0) || any(app.gaBetaInit > 2)
                    uialert(f, ['beta_init must be in (0, 2]. ' ...
                        'Use 1 for ordinary exponential.'], ...
                        'Invalid input');  return;
                end
                if app.gaFWHM <= 0
                    uialert(f, 'IRF FWHM must be positive.', ...
                        'Invalid input');  return;
                end
                
                lblStatus.Text = 'Fitting...';  drawnow;
                try
                    [tauOut, betaOut, t0Out, fwhmOut, A, fitM, info] = ...
                        app.fitGlobalAnalysis( ...
                        app.deltaA, app.delay, ...
                        app.gaTauInit, app.gaT0, app.gaFWHM, ...
                        app.gaTauFixed, app.gaT0Fixed, app.gaFWHMFixed, ...
                        app.gaHasInf, ...
                        app.gaBetaInit, app.gaStretchOn, app.gaIRFTreat, ...
                        [app.gaFitTMin, app.gaFitTMax], app.gaCA);
                catch ME
                    uialert(f, sprintf('Fit failed:\n%s', ME.message), ...
                        'Fit error');
                    lblStatus.Text = 'Fit failed.';
                    app.clearGaResults();
                    refreshSeqRowState();
                    return;
                end
                
                app.gaResultTau    = tauOut;
                app.gaResultBeta   = betaOut;
                app.gaResultT0     = t0Out;
                app.gaResultFWHM   = fwhmOut;
                app.gaResultHasInf = app.gaHasInf;
                app.gaResultDADS   = A;
                app.gaResultFit    = fitM;
                app.gaResultRMS    = info.rms;
                app.gaResultTauStd = info.tauStd;
                % Settings this result was obtained with: the table /
                % dropdown may be edited afterwards without a new fit.
                app.gaResultStretchOn = app.gaStretchOn;
                app.gaResultIRF       = app.gaIRFTreat;
                app.gaResultTRange    = [app.gaFitTMin, app.gaFitTMax];
                if info.caOrder > 0
                    app.gaResultCA = struct('order', info.caOrder, ...
                        'spec', info.caSpec, 'fwhm', info.caFwhm, ...
                        'fwhmStd', info.caFwhmStd, 'mode', info.caMode);
                else
                    app.gaResultCA = [];
                end
                
                % EADS (1st-order Bateman) vs SAS (sequential, possibly
                % stretched) -- pick the right one based on whether any
                % component is stretched.
                useStretched = any(app.gaStretchOn);
                sasErr = '';
                try
                    if useStretched
                        app.gaResultEADS = app.computeSASsequential( ...
                            app.gaDataNoCA(), app.delay, ...
                            tauOut, betaOut, app.gaStretchOn, ...
                            app.gaHasInf, t0Out, fwhmOut, ...
                            app.gaIRFTreat, app.gaResultTRange, ...
                            app.gaSeqModel);
                    else
                        app.gaResultEADS = app.computeEADSFromDADS( ...
                            A, tauOut, app.gaHasInf);
                    end
                catch ME
                    app.gaResultEADS = [];
                    sasErr = ME.message;
                end
                if ~isempty(app.gaResultEADS) && ...
                        all(~isfinite(app.gaResultEADS(:)))
                    sasErr = 'all-NaN result';
                    app.gaResultEADS = [];
                end
                
                % Update result text
                lines = {};
                lines{end+1} = sprintf('Fit converged (%d iters)', info.iters);
                lines{end+1} = sprintf('Initial RMS : %.4g', info.initialRMS);
                lines{end+1} = sprintf('Final   RMS : %.4g', info.rms);
                if info.rms >= 0.99 * info.initialRMS && info.iters > 2
                    lines{end+1} = '';
                    lines{end+1} = '** NOT CONVERGED - try different';
                    lines{end+1} = '   initial tau values, or fix IRF/t0.';
                end
                lines{end+1} = '';
                lines{end+1} = 'Time constants:';
                for ii = 1:numel(tauOut)
                    if app.gaStretchOn(ii)
                        tag = sprintf(' [stretched, beta=%.3f]', betaOut(ii));
                        if isfield(info, 'betaStd') && ...
                                isfinite(info.betaStd(ii))
                            tag = sprintf( ...
                                ' [stretched, beta=%.3f+/-%.3f]', ...
                                betaOut(ii), info.betaStd(ii));
                        end
                    else
                        tag = '';
                    end
                    if app.gaTauFixed(ii)
                        tag = [tag, ' (tau fixed)']; %#ok<AGROW>
                    end
                    tag = [tag, app.undetNote(info, 'tauStd', ii)]; %#ok<AGROW>
                    if isfield(info, 'tauStd') && isfinite(info.tauStd(ii))
                        lines{end+1} = sprintf( ...
                            '  tau_%d = %10.4g +/- %.2g %s%s', ...
                            ii, tauOut(ii), info.tauStd(ii), ...
                            app.tUnitTxt(), tag);  %#ok<AGROW>
                    else
                        lines{end+1} = sprintf('  tau_%d = %10.4g %s%s', ...
                            ii, tauOut(ii), app.tUnitTxt(), tag);  %#ok<AGROW>
                    end
                end
                if app.gaHasInf
                    lines{end+1} = '  tau_inf  (offset component)';
                end
                lines{end+1} = '';
                if isfield(info, 't0Std') && isfinite(info.t0Std)
                    lines{end+1} = sprintf('IRF t0   = %.4g +/- %.2g %s', ...
                        t0Out, info.t0Std, app.tUnitTxt());
                else
                    lines{end+1} = sprintf('IRF t0   = %.4g %s%s', ...
                        t0Out, app.tUnitTxt(), ...
                        [ternary(app.gaT0Fixed, ' (fixed)', ''), ...
                         app.undetNote(info, 't0Std', 1)]);
                end
                if isfield(info, 'fwhmAtLimit') && info.fwhmAtLimit
                    lines{end+1} = sprintf('IRF FWHM = %.4g %s  [not resolved]', ...
                        fwhmOut, app.tUnitTxt());
                elseif isfield(info, 'fwhmStd') && isfinite(info.fwhmStd)
                    lines{end+1} = sprintf('IRF FWHM = %.4g +/- %.2g %s', ...
                        fwhmOut, info.fwhmStd, app.tUnitTxt());
                else
                    lines{end+1} = sprintf('IRF FWHM = %.4g %s%s', ...
                        fwhmOut, app.tUnitTxt(), ...
                        [ternary(app.gaFWHMFixed, ' (fixed)', ''), ...
                         app.undetNote(info, 'fwhmStd', 1)]);
                end
                lines = [lines, app.irfLimitNote(info)];
                lines{end+1} = ['IRF model: ' ternary( ...
                    strcmp(app.gaIRFTreat, 'none'), ...
                    'none (pure exponentials)', 'numerical conv.')];
                if info.caOrder > 0
                    caNm = {'Gaussian', '1st derivative', '2nd derivative'};
                    lines{end+1} = '';
                    lines{end+1} = sprintf( ...
                        'Coherent artifact: %d term(s)', info.caOrder);
                    switch info.caMode
                        case 'irf'
                            lines{end+1} = '  width = IRF FWHM';
                        case 'fixed'
                            lines{end+1} = sprintf( ...
                                '  width = %.4g %s (fixed)', ...
                                info.caFwhm, app.tUnitTxt());
                        otherwise
                            if isfinite(info.caFwhmStd)
                                lines{end+1} = sprintf( ...
                                    '  width = %.4g +/- %.2g %s', ...
                                    info.caFwhm, info.caFwhmStd, ...
                                    app.tUnitTxt());
                            else
                                lines{end+1} = sprintf( ...
                                    '  width = %.4g %s%s', info.caFwhm, ...
                                    app.tUnitTxt(), ...
                                    app.undetNote(info, 'caFwhmStd', 1));
                            end
                    end
                    for ii = 1:info.caOrder
                        [~, iMx] = max(abs(info.caSpec(:, ii)));
                        lines{end+1} = sprintf( ...
                            '  %-14s largest %.3g at %.4g %s', ...
                            caNm{ii}, info.caSpec(iMx, ii), ...
                            app.wavelength(iMx), app.xUnitTxt()); %#ok<AGROW>
                    end
                    if any(tauOut < fwhmOut)
                        lines{end+1} = '  A tau below the IRF FWHM is';
                        lines{end+1} = '  hard to tell from the artifact.';
                    end
                    if ~strcmp(info.caMode, 'irf') && ~app.gaFWHMFixed ...
                            && info.caOrder >= 2
                        lines{end+1} = '  With a width of their own the';
                        lines{end+1} = '  derivative terms absorb a change of';
                        lines{end+1} = '  the IRF FWHM: it is then hardly';
                        lines{end+1} = '  determined. Fix it, or use';
                        lines{end+1} = '  "Same as IRF".';
                    end
                    if (~app.gaFWHMFixed && (fwhmOut > 3 * app.gaFWHM || ...
                            fwhmOut < app.gaFWHM / 3)) || ...
                            (~app.gaT0Fixed && abs(t0Out - app.gaT0) > ...
                            3 * max(app.gaFWHM, fwhmOut))
                        lines{end+1} = '  t0 or FWHM moved far from the start';
                        lines{end+1} = '  value. With too few components the';
                        lines{end+1} = '  artifact terms can leave time zero';
                        lines{end+1} = '  and stand in for kinetics: check the';
                        lines{end+1} = '  spectra, or fix the IRF.';
                    end
                end
                if useStretched
                    lines{end+1} = '';
                    lines{end+1} = 'Right-panel "EADS" -> "SAS"';
                    lines{end+1} = '(sequential w/ stretched; clock';
                    lines{end+1} = 'set by "Stretched step").';
                    if ~isempty(sasErr)
                        lines{end+1} = '';
                        lines{end+1} = '*** SAS computation FAILED ***';
                        lines{end+1} = sasErr;
                    end
                end
                txtResult.Value = lines;
                
                lblStatus.Text = sprintf('Done. RMS = %.3g', info.rms);
                plotResults();
                updateEadsPanelTitle();
                refreshSeqRowState();
            end
            
            function onKinScaleChange()
                plotKinetics();
                plotPopulation();
            end
            
            function setFitTFull()
                edFitTmin.Value = min(app.delay);
                edFitTmax.Value = max(app.delay);
            end
            
            function onIRFModeChange()
                % 'none' removes the IRF from the model entirely: FWHM
                % is meaningless and t0 is degenerate with amplitudes,
                % so both are locked (t0 stays editable as a fixed
                % onset value).
                isNone = strcmp(ddIRF.Value, 'none');
                if isNone
                    edFW.Enable    = 'off';
                    cbFWfix.Enable = 'off';
                    cbT0fix.Value  = true;
                    cbT0fix.Enable = 'off';
                else
                    edFW.Enable    = 'on';
                    cbFWfix.Enable = 'on';
                    cbT0fix.Enable = 'on';
                end
                syncCAEnable();
            end
            
            function syncCAEnable()
                % The artifact terms need the IRF (time zero and a
                % width); their own width is used only when asked for.
                hasIRF = ~strcmp(ddIRF.Value, 'none');
                useCA  = hasIRF && ddCA.Value > 0;
                ddCA.Enable   = app.pick(hasIRF, 'on', 'off');
                ddCAW.Enable  = app.pick(useCA, 'on', 'off');
                lblCAW.Enable = app.pick(useCA, 'on', 'off');
                edCAW.Enable  = app.pick(useCA && ...
                    ~strcmp(ddCAW.Value, 'irf'), 'on', 'off');
            end
            
            function onSeqModelChange()
                % Other reading of the stretched steps: the fit stays as
                % it is; the SAS of the stored result are recomputed with
                % the settings of that result, and the SAS and population
                % panels are redrawn.
                if strcmp(app.gaSeqModel, ddSeq.Value), return; end
                app.gaSeqModel = ddSeq.Value;
                if isempty(app.gaResultTau) || ~any(resStretch())
                    return;
                end
                lblStatus.Text = 'Recomputing SAS...';  drawnow;
                try
                    E = app.computeSASsequential( ...
                        app.gaDataNoCA(), app.delay, ...
                        app.gaResultTau, app.gaResultBeta, ...
                        resStretch(), app.gaResultHasInf, ...
                        app.gaResultT0, app.gaResultFWHM, ...
                        resIRF(), app.gaResultTRange, app.gaSeqModel);
                    if all(~isfinite(E(:)))
                        error('TRSpec:sasAllNaN', 'all-NaN result');
                    end
                    app.gaResultEADS = E;
                    lblStatus.Text = sprintf('SAS recomputed (%s).', ...
                        seqModelTxt());
                catch ME
                    app.gaResultEADS = [];
                    lblStatus.Text = sprintf( ...
                        'SAS computation failed: %s', ME.message);
                end
                plotEADS();
                plotPopulation();
                updateEadsPanelTitle();
            end
            
            function s = seqModelTxt()
                % Short name of the selected reading of stretched steps.
                if strcmp(app.gaSeqModel, 't0')
                    s = 'clock from t0';
                else
                    s = 'clock from entry';
                end
            end
            
            function [popC, popT, nSp] = computePopulation()
                % Sequential-model populations c_i(t) matching the EADS
                % decomposition. For the all-standard-exp case use the
                % Bateman matrix on the SAME basis columns as the fit
                % (linearity: convolution with the IRF commutes with
                % the Bateman superposition). With stretched components
                % use the sequential chain of computeSASsequential.
                popC = [];  popT = [];  nSp = 0;
                if isempty(app.gaResultTau), return; end
                popT = app.delay(:);
                if any(resStretch())
                    try
                        [~, popC] = app.computeSASsequential( ...
                            app.gaDataNoCA(), app.delay, ...
                            app.gaResultTau, app.gaResultBeta, ...
                            resStretch(), app.gaResultHasInf, ...
                            app.gaResultT0, app.gaResultFWHM, ...
                            resIRF(), app.gaResultTRange, ...
                            app.gaSeqModel);
                    catch
                        popC = [];
                    end
                else
                    try
                        [~, tauSorted, Bm] = app.computeEADSFromDADS( ...
                            app.gaResultDADS, app.gaResultTau, ...
                            app.gaResultHasInf);
                        if app.gaResultHasInf
                            tauCols = [tauSorted(:); Inf];
                        else
                            tauCols = tauSorted(:);
                        end
                        basis = app.buildGaBasisExt(popT, ...
                            tauCols.', ones(1, numel(tauCols)), ...
                            false(1, numel(tauCols)), ...
                            app.gaResultT0, app.gaResultFWHM, ...
                            false, resIRF());
                        popC = (Bm * basis.').';   % N x nSpecies
                    catch
                        popC = [];
                    end
                end
                nSp = size(popC, 2);
            end
            
            function plotPopulation()
                cla(axPop);
                [popC, popT, nSp] = computePopulation();
                if isempty(popC)
                    if isempty(app.gaResultTau)
                        title(axPop, 'EADS populations c_i(t)');
                    else
                        title(axPop, 'Populations c_i(t): not available');
                    end
                    return;
                end
                cmapP = lines(max(nSp, 3));
                hold(axPop, 'on');
                for jj = 1:nSp
                    plot(axPop, popT, popC(:, jj), '-', ...
                        'LineWidth', 1.4, 'Color', cmapP(jj, :), ...
                        'DisplayName', sprintf('S%d', jj));
                end
                hold(axPop, 'off');
                if strcmp(ddKinScale.Value, 'Log')
                    tPos = popT(popT > 0);
                    if ~isempty(tPos)
                        axPop.XScale = 'log';
                        xlP = [min(tPos), max(popT)];
                    else
                        xlP = [min(popT), max(popT)];
                    end
                else
                    axPop.XScale = 'linear';
                    xlP = [min(popT), max(popT)];
                end
                frP = getFitTRangeIfRestricted();
                if ~isempty(frP)
                    if strcmp(axPop.XScale, 'log')
                        frP(1) = max(frP(1), xlP(1));
                    end
                    if frP(1) < frP(2), xlP = frP; end
                end
                xlim(axPop, xlP);
                ylim(axPop, 'auto');
                legend(axPop, 'Location', 'best', 'Box', 'off');
                if any(resStretch())
                    title(axPop, 'SAS populations c_i(t)');
                else
                    title(axPop, 'EADS populations c_i(t)');
                end
            end
            
            function exportPopulation()
                [popC, popT, nSp] = computePopulation();
                if isempty(popC)
                    uialert(f, 'Run a fit first.', 'No result');
                    return;
                end
                [fn, fp] = app.putFile('*.csv', ...
                    'Export population kinetics', 'ga_populations.csv');
                if isequal(fn, 0), return; end
                hdr = ['t_' app.tUnitTxt()];
                for jj = 1:nSp
                    hdr = [hdr, sprintf(',S%d', jj)]; %#ok<AGROW>
                end
                try
                    fid = fopen(fullfile(fp, fn), 'w');
                    if fid < 0
                        error('Cannot write %s', fullfile(fp, fn));
                    end
                    fprintf(fid, '%s\n', hdr);
                    fclose(fid);
                    writematrix([popT, popC], fullfile(fp, fn), ...
                        'WriteMode', 'append');
                    lblStatus.Text = sprintf('Populations saved: %s', fn);
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Export error');
                end
            end
            
            function so = resStretch()
                % Stretched flags that belong to the STORED result. The
                % live table settings (app.gaStretchOn) can differ once
                % the table or N has been edited without a new fit.
                if isempty(app.gaResultTau)
                    so = logical(app.gaStretchOn(:).');   % nothing stored
                    return;
                end
                so = app.gaResultStretchOn;
                if numel(so) ~= numel(app.gaResultTau)
                    so = app.gaStretchOn;       % older session files
                end
                if numel(so) ~= numel(app.gaResultTau)
                    so = false(1, numel(app.gaResultTau));
                end
                so = logical(so(:).');
            end
            
            function m = resIRF()
                % IRF model of the stored result (see resStretch).
                m = app.gaResultIRF;
                if isempty(m), m = app.gaIRFTreat; end
                m = app.irfModeNorm(m);
            end
            
            function updateEadsPanelTitle()
                % EADS axes title becomes "SAS" when any component is
                % stretched (EADS Bateman doesn't apply to stretched).
                if ~isvalid(axEADS) || ~isvalid(axEADSnorm), return; end
                if any(resStretch())
                    title(axEADS, ['SAS (sequential; stretched steps: ' ...
                        seqModelTxt() ')']);
                    title(axEADSnorm, 'SAS (peak-normalized)');
                    if isvalid(btnExpEADS), btnExpEADS.Text = 'SAS'; end
                else
                    title(axEADS, ...
                        'Evolution-associated difference spectra (sequential model)');
                    title(axEADSnorm, 'EADS (peak-normalized)');
                    if isvalid(btnExpEADS), btnExpEADS.Text = 'EADS'; end
                end
            end
            
            function plotResults()
                if isempty(app.gaResultDADS), return; end
                rebuildCompCheckboxes();
                plotDADS();
                plotEADS();
                plotMaps();
                plotKinetics();
                plotPopulation();
            end
            
            function onCompVis()
                plotDADS();
                plotEADS();
            end
            
            function syncKinWLFromMain()
                curKinWL = app.selWL;
                edKinWL.Value = curKinWL;
                plotKinetics();
            end
            
            function onKinWLChange(val)
                curKinWL = val;
                plotKinetics();
            end
            
            function rebuildCompCheckboxes()
                delete(cbHost.Children);
                if isempty(app.gaResultDADS), return; end
                nc    = size(app.gaResultDADS, 2);
                nMain = numel(app.gaResultTau);
                hasCA = isstruct(app.gaResultCA) && ...
                    isfield(app.gaResultCA, 'spec') && ...
                    ~isempty(app.gaResultCA.spec);
                cbCAshow = [];
                cg = uigridlayout(cbHost, [1, nc + 2 + double(hasCA)]);
                colW = [{'fit'}, repmat({'fit'}, 1, nc + double(hasCA)), {'1x'}];
                cg.ColumnWidth = colW;
                cg.Padding = [0 4 0 4];  cg.ColumnSpacing = 10;
                uilabel(cg, 'Text', 'Show components:', 'FontAngle', 'italic');
                cbComps = gobjects(1, nc);
                for ii = 1:nMain
                    txt = app.uiStr(sprintf('\\tau_%d = %.3g %s', ii, ...
                        app.gaResultTau(ii), app.tUnitAx()));
                    cbComps(ii) = uicheckbox(cg, 'Text', txt, ...
                        'Value', true, ...
                        'ValueChangedFcn', @(~,~) onCompVis());
                end
                if app.gaResultHasInf
                    cbComps(nc) = uicheckbox(cg, 'Text', app.uiStr('\tau = \infty'), ...
                        'Value', true, 'ValueChangedFcn', @(~,~) onCompVis());
                end
                if hasCA
                    cbCAshow = uicheckbox(cg, 'Text', 'artifact', ...
                        'Value', true, ...
                        'Tooltip', ['Show the spectra of the coherent-' ...
                            'artifact terms (dashed) in the DADS plot.'], ...
                        'ValueChangedFcn', @(~,~) onCompVis());
                end
                uilabel(cg, 'Text', '');  % right spacer
            end
            
            function vis = getCompVisible()
                if isempty(cbComps) || ~all(arrayfun(@isvalid, cbComps))
                    vis = true(1, size(app.gaResultDADS, 2));
                else
                    vis = arrayfun(@(cb) logical(cb.Value), cbComps);
                end
            end
            
            function plotDADS()
                if isempty(app.gaResultDADS), return; end
                A     = app.gaResultDADS;
                nc    = size(A, 2);
                nMain = numel(app.gaResultTau);
                cmap  = lines(max(8, nc));
                vis   = getCompVisible();
                
                % ---- raw DADS ----
                cla(axDADS);  hold(axDADS, 'on');
                legH = gobjects(0);  legS = {};
                for ii = 1:nc
                    if ~vis(ii), continue; end
                    h = plot(axDADS, app.wavelength, A(:, ii), '-', ...
                        'Color', cmap(ii, :), 'LineWidth', 1.4);
                    legH(end+1) = h;  %#ok<AGROW>
                    if ii <= nMain
                        legS{end+1} = sprintf('\\tau_%d = %.3g %s', ...
                            ii, app.gaResultTau(ii), app.tUnitAx());  %#ok<AGROW>
                    else
                        legS{end+1} = '\tau = \infty';  %#ok<AGROW>
                    end
                end
                % Spectra of the coherent-artifact terms (dashed): part
                % of the fit, not of the DADS.
                caR = app.gaResultCA;
                showCA = isempty(cbCAshow) || ~isvalid(cbCAshow) || ...
                    cbCAshow.Value;
                if showCA && isstruct(caR) && isfield(caR, 'spec') && ...
                        size(caR.spec, 1) == numel(app.wavelength)
                    caNm  = {'artifact: Gaussian', 'artifact: 1st deriv.', ...
                             'artifact: 2nd deriv.'};
                    caCol = [0.30 0.30 0.30; 0.50 0.50 0.50; 0.68 0.68 0.68];
                    for ii = 1:min(3, size(caR.spec, 2))
                        h = plot(axDADS, app.wavelength, caR.spec(:, ii), ...
                            '--', 'Color', caCol(ii, :), 'LineWidth', 1.1);
                        legH(end+1) = h;  %#ok<AGROW>
                        legS{end+1} = caNm{ii};  %#ok<AGROW>
                    end
                end
                yline(axDADS, 0, 'k:');
                hold(axDADS, 'off');  grid(axDADS, 'on');
                xlabel(axDADS, app.xAxLabel());
                ylabel(axDADS, 'Amplitude (\DeltaA or Intensity)');
                title(axDADS, sprintf('Decay-associated difference spectra   (RMS = %.3g)', ...
                    app.gaResultRMS));
                if ~isempty(legS)
                    legend(axDADS, legH, legS, 'Location', 'best');
                else
                    legend(axDADS, 'off');
                end
                xlim(axDADS, [min(app.wavelength) max(app.wavelength)]);
                
                % ---- normalized DADS (peak-normalized per component) ----
                cla(axDADSnorm);  hold(axDADSnorm, 'on');
                legHN = gobjects(0);  legSN = {};
                for ii = 1:nc
                    if ~vis(ii), continue; end
                    col = A(:, ii);
                    pk  = max(abs(col), [], 'omitnan');
                    if isempty(pk) || pk == 0 || isnan(pk), pk = 1; end
                    h = plot(axDADSnorm, app.wavelength, col / pk, '-', ...
                        'Color', cmap(ii, :), 'LineWidth', 1.4);
                    legHN(end+1) = h;  %#ok<AGROW>
                    if ii <= nMain
                        legSN{end+1} = sprintf('\\tau_%d', ii);  %#ok<AGROW>
                    else
                        legSN{end+1} = '\tau = \infty';  %#ok<AGROW>
                    end
                end
                yline(axDADSnorm, 0, 'k:');
                hold(axDADSnorm, 'off');  grid(axDADSnorm, 'on');
                xlabel(axDADSnorm, app.xAxLabel());
                ylabel(axDADSnorm, 'Normalized amplitude');
                title(axDADSnorm, 'DADS (peak-normalized)');
                if ~isempty(legSN)
                    legend(axDADSnorm, legHN, legSN, 'Location', 'best');
                else
                    legend(axDADSnorm, 'off');
                end
                xlim(axDADSnorm, [min(app.wavelength) max(app.wavelength)]);
                ylim(axDADSnorm, [-1.1 1.1]);
            end
            
            function plotEADS()
                % Sequential-model spectra: axEADS (raw) + axEADSnorm (normalized).
                % Same component color as DADS is used when the species and
                % the fit-order DADS column share the same tau value (so the
                % visibility checkboxes still control both plots coherently).
                
                cla(axEADS);  cla(axEADSnorm);
                if isempty(app.gaResultEADS)
                    if any(resStretch())
                        title(axEADS, ['Species-associated spectra ' ...
                            '(sequential, with stretched) - run fit to compute']);
                        title(axEADSnorm, 'SAS (peak-normalized)');
                    else
                        title(axEADS, ['Evolution-associated difference spectra  ' ...
                            '(sequential) - run fit to compute']);
                        title(axEADSnorm, 'EADS (peak-normalized)');
                    end
                    return;
                end
                
                E     = app.gaResultEADS;
                nc    = size(E, 2);
                nMain = numel(app.gaResultTau);
                cmap  = lines(max(8, nc));
                
                % Sorted tau order (ascending) matches EADS column order.
                [tauSorted, sortIdx] = sort(app.gaResultTau(:), 'ascend');
                % Map sorted species index -> fit-order index (for checkbox
                % lookup and color consistency with DADS).
                if app.gaResultHasInf
                    permFit = [sortIdx(:); nMain + 1];
                else
                    permFit = sortIdx(:);
                end
                vis = getCompVisible();   % in fit-order
                
                % ---- raw EADS ----
                hold(axEADS, 'on');
                legH = gobjects(0);  legS = {};
                for iS = 1:nc
                    fitIdx = permFit(iS);
                    if fitIdx <= numel(vis) && ~vis(fitIdx), continue; end
                    col = E(:, iS);
                    h = plot(axEADS, app.wavelength, col, '-', ...
                        'Color', cmap(fitIdx, :), 'LineWidth', 1.4);
                    legH(end+1) = h;  %#ok<AGROW>
                    if iS <= nMain
                        legS{end+1} = sprintf('S_%d  \\tau = %.3g %s', ...
                            iS, tauSorted(iS), app.tUnitAx());  %#ok<AGROW>
                    else
                        legS{end+1} = sprintf('S_%d  \\tau = \\infty', iS);  %#ok<AGROW>
                    end
                end
                yline(axEADS, 0, 'k:');
                hold(axEADS, 'off');  grid(axEADS, 'on');
                xlabel(axEADS, app.xAxLabel());
                ylabel(axEADS, 'Amplitude (\DeltaA or Intensity)');
                if ~isempty(legS)
                    legend(axEADS, legH, legS, 'Location', 'best');
                else
                    legend(axEADS, 'off');
                end
                xlim(axEADS, [min(app.wavelength) max(app.wavelength)]);
                % Compute ylim explicitly from the visible columns so
                % stale autolimits / legend offsets can't squish the
                % traces against an empty band.
                visIdx = false(1, nc);
                for iS = 1:nc
                    fitIdx = permFit(iS);
                    if fitIdx > numel(vis) || vis(fitIdx)
                        visIdx(iS) = true;
                    end
                end
                if any(visIdx)
                    Eshown = E(:, visIdx);
                    yMin = min(Eshown(:), [], 'omitnan');
                    yMax = max(Eshown(:), [], 'omitnan');
                    if isfinite(yMin) && isfinite(yMax) && yMin < yMax
                        pad = 0.08 * (yMax - yMin);
                        ylim(axEADS, [yMin - pad, yMax + pad]);
                    else
                        ylim(axEADS, 'auto');
                    end
                else
                    ylim(axEADS, 'auto');
                end
                
                % ---- normalized EADS ----
                hold(axEADSnorm, 'on');
                legHN = gobjects(0);  legSN = {};
                for iS = 1:nc
                    fitIdx = permFit(iS);
                    if fitIdx <= numel(vis) && ~vis(fitIdx), continue; end
                    col = E(:, iS);
                    pk  = max(abs(col), [], 'omitnan');
                    if isempty(pk) || pk == 0 || isnan(pk), pk = 1; end
                    h = plot(axEADSnorm, app.wavelength, col / pk, '-', ...
                        'Color', cmap(fitIdx, :), 'LineWidth', 1.4);
                    legHN(end+1) = h;  %#ok<AGROW>
                    legSN{end+1} = sprintf('S_%d', iS);  %#ok<AGROW>
                end
                yline(axEADSnorm, 0, 'k:');
                hold(axEADSnorm, 'off');  grid(axEADSnorm, 'on');
                xlabel(axEADSnorm, app.xAxLabel());
                ylabel(axEADSnorm, 'Normalized amplitude');
                if ~isempty(legSN)
                    legend(axEADSnorm, legHN, legSN, 'Location', 'best');
                else
                    legend(axEADSnorm, 'off');
                end
                xlim(axEADSnorm, [min(app.wavelength) max(app.wavelength)]);
                ylim(axEADSnorm, [-1.1 1.1]);
                % Titles (EADS / SAS, and the reading of stretched steps)
                % come from one place, whatever triggered the redraw.
                updateEadsPanelTitle();
            end
            
            function plotMaps()
                if isempty(app.gaResultFit), return; end
                buildMapAxes();   % (re)creates axes inside the panels
                
                wl = app.wavelength(:).';
                t  = app.delay(:).';
                
                cmap = app.getColormap();
                cl = max(abs(app.deltaA(:)), [], 'omitnan');
                if isempty(cl) || cl == 0 || isnan(cl), cl = 1; end
                
                R = app.deltaA - app.gaResultFit;
                clr = max(abs(R(:)), [], 'omitnan');
                if isempty(clr) || clr == 0 || isnan(clr), clr = cl/10; end
                resRMS = sqrt(mean(R(isfinite(R)).^2));
                
                drawMapStack(mapAxExp, cbarAxExp, t, wl, app.deltaA.',    [-cl cl],  cmap, true);
                drawMapStack(mapAxFit, cbarAxFit, t, wl, app.gaResultFit.', [-cl cl],  cmap, false);
                drawMapStack(mapAxRes, cbarAxRes, t, wl, R.',             [-clr clr], cmap, true);
                
                % Panels carry the human-readable titles (axes-level titles
                % are hidden in split mode to save vertical space).
                mapPanelExp.Title = app.uiStr('Experimental  (click to pick \lambda)');
                mapPanelFit.Title = 'Global fit';
                mapPanelRes.Title = sprintf('Residual  (RMS = %.3g)', resRMS);
                
                % Install click pickers and crosshair on EVERY axes.
                allAx = [mapAxExp(:); mapAxFit(:); mapAxRes(:)];
                hMapVLine = gobjects(1, numel(allAx));
                for kk = 1:numel(allAx)
                    ax = allAx(kk);
                    kids = ax.Children;
                    for jj = 1:numel(kids)
                        if isa(kids(jj), 'matlab.graphics.primitive.Surface')
                            kids(jj).ButtonDownFcn = @(~,ev) onMapClick(ev);
                        end
                    end
                    ax.ButtonDownFcn = @(~,ev) onMapClick(ev);
                    hMapVLine(kk) = xline(ax, curKinWL, 'k-', ...
                        'LineWidth', 1.0, 'Alpha', 0.7, ...
                        'HitTest', 'off', 'PickableParts', 'none');
                end
            end
            
            function buildMapAxes()
                % (Re)build the axes inside mapPanelExp/Fit/Res according
                % to the currently selected map scale. Each panel hosts a
                % [nRows x 2] uigridlayout: column 1 holds the data axes
                % (one per delay-band), column 2 is a fixed-width host
                % for the colorbar. This guarantees all three data axes
                % have the same horizontal extent, regardless of whether
                % a particular map actually draws its own colorbar.
                mode = ddMapScale.Value;    % 'linear' | 'log' | 'split'
                if strcmp(mode, 'split')
                    nRows = 2;
                else
                    nRows = 1;
                end
                
                hosts   = [mapPanelExp, mapPanelFit, mapPanelRes];
                dataOut = cell(1, 3);
                cbarOut = cell(1, 3);
                for h = 1:3
                    hp = hosts(h);
                    delete(hp.Children);
                    g = uigridlayout(hp, [nRows, 2]);
                    g.ColumnWidth   = {'1x', 50};
                    g.RowHeight     = repmat({'1x'}, 1, nRows);
                    g.Padding       = [0 0 0 0];
                    g.ColumnSpacing = 2;
                    g.RowSpacing    = 2;
                    
                    dAx = gobjects(1, nRows);
                    cAx = gobjects(1, nRows);
                    for r = 1:nRows
                        dAx(r) = uiaxes(g);
                        dAx(r).Layout.Row    = r;
                        dAx(r).Layout.Column = 1;
                        dAx(r).Color = [0.93 0.93 0.93];
                        xlabel(dAx(r), app.xAxLabel());
                        ylabel(dAx(r), sprintf('Delay (%s)', ...
                            app.tUnitAx()));
                        
                        cAx(r) = uiaxes(g);
                        cAx(r).Layout.Row    = r;
                        cAx(r).Layout.Column = 2;
                        cAx(r).Visible = 'off';
                    end
                    dataOut{h} = dAx;
                    cbarOut{h} = cAx;
                end
                mapAxExp = dataOut{1};  cbarAxExp = cbarOut{1};
                mapAxFit = dataOut{2};  cbarAxFit = cbarOut{2};  %#ok<NASGU>
                mapAxRes = dataOut{3};  cbarAxRes = cbarOut{3};
            end
            
            function drawMapStack(axArr, cbarArr, t, wl, C, clim, cmap, showCbar)
                % Draw the (N_delay x N_wl) color map C on either a single
                % axes or a stacked pair of axes according to the current
                % mode. If showCbar is true, a colorbar is hosted on the
                % first cbar cell (the top cell in split mode).
                mode = ddMapScale.Value;
                thr  = max(edMapThr.Value, eps);
                tmin = edMapTmin.Value;
                tmax = edMapTmax.Value;
                if ~isfinite(tmin), tmin = min(t); end
                if ~isfinite(tmax), tmax = max(t); end
                if tmax <= tmin, tmax = max(t); end
                % When the fit was restricted to a delay sub-range, the
                % maps follow it: clamp the displayed window to the fit
                % range (the user's Map t range can still narrow it
                % further, never widen it).
                frM = getFitTRangeIfRestricted();
                if ~isempty(frM)
                    tmin = max(tmin, frM(1));
                    tmax = min(tmax, frM(2));
                    if tmax <= tmin
                        tmin = frM(1);  tmax = frM(2);
                    end
                end
                
                switch mode
                    case 'linear'
                        masks  = {(t >= tmin) & (t <= tmax)};
                        scales = {'linear'};
                        ylims  = {[tmin, tmax]};
                    case 'log'
                        % In log mode tmin only matters when the fit
                        % range restricts it (log otherwise starts at
                        % the first positive delay); tmax is always an
                        % upper bound.
                        pos = (t > 0) & (t >= tmin) & (t <= tmax);
                        masks  = {pos};
                        scales = {'log'};
                        if any(pos)
                            ylims = {[min(t(pos)), max(t(pos))]};
                        else
                            ylims = {[1, 10]};
                        end
                    case 'split'
                        % Bottom: linear panel from tmin to thr.
                        % Top:    log panel from thr to tmax.
                        bottomMask = (t <= thr) & (t >= tmin);
                        topMask    = (t >= thr) & (t <= tmax);
                        masks  = {topMask,    bottomMask};
                        scales = {'log',      'linear'};
                        ylims  = {[thr, tmax], [tmin, thr]};
                end
                
                for k = 1:numel(axArr)
                    ax = axArr(k);
                    cla(ax);
                    % Reset tick modes so leftover manual ytick labels
                    % from a previous draw cannot bleed into the new
                    % range. We re-apply the leading-space padding at
                    % the bottom of this function from the auto ticks.
                    ax.YTickMode      = 'auto';
                    ax.YTickLabelMode = 'auto';
                    ax.XTickMode      = 'auto';
                    ax.XTickLabelMode = 'auto';
                    
                    if k > numel(masks), continue; end
                    m = masks{k};
                    set(ax, 'YScale', scales{k});
                    
                    if any(m)
                        tSub = t(m);
                        cSub = C(m, :);
                        [W, T] = meshgrid(wl, tSub);
                        surface(ax, W, T, zeros(size(cSub)), cSub, ...
                            'EdgeColor', 'none', 'FaceColor', 'flat');
                        view(ax, 2);
                        ax.YDir = 'normal';
                        ax.CLim = clim;
                        colormap(ax, cmap);
                    end
                    xlim(ax, [min(wl), max(wl)]);
                    
                    % Always set ylim from the requested range, regardless
                    % of whether any data points fell inside the mask.
                    yl = ylims{k};
                    if all(isfinite(yl)) && yl(1) < yl(2)
                        if strcmp(scales{k}, 'log') && yl(1) <= 0
                            yl(1) = max(eps, yl(1));
                        end
                        ylim(ax, yl);
                    end
                    
                    if k == numel(axArr)
                        xlabel(ax, app.xAxLabel());
                    else
                        xlabel(ax, '');
                        ax.XTickLabel = [];
                    end
                    ylabel(ax, sprintf('Delay (%s)', app.tUnitAx()));
                end
                
                % Set the final ylim FIRST (before label padding) so
                % MATLAB recomputes ytick positions appropriate to the
                % requested [tmin, tmax] window, not the surface data
                % bounding box.
                for k = 1:numel(axArr)
                    if k > numel(ylims), continue; end
                    yl = ylims{k};
                    if all(isfinite(yl)) && yl(1) < yl(2)
                        if strcmp(scales{k}, 'log') && yl(1) <= 0
                            yl(1) = max(eps, yl(1));
                        end
                        ylim(axArr(k), yl);
                        % Cap the number of ticks on linear axes so
                        % they don't pile on top of each other when
                        % the cell is short (e.g. the bottom row of
                        % split view, which is much shorter than the
                        % top log-scaled row).  We pick "nice"
                        % rounded steps automatically.
                        if strcmp(scales{k}, 'linear')
                            niceTicks = chooseNiceTicks(yl(1), yl(2), 4);
                            if ~isempty(niceTicks)
                                axArr(k).YTick = niceTicks;
                            end
                        end
                    end
                end
                drawnow limitrate;
                
                % Force identical left-margin across all rows by
                % padding the ytick labels with leading spaces so
                % they all have the same visual width. This stops
                % the wavelength axis from drifting when one row's
                % ticks are short ("0,1,2") and another's are wide
                % ("10^{2}" or "-0.5,1.5").  Leading spaces DO take
                % up horizontal space in MATLAB tick rendering.
                %
                % Critical: capture YTick (positions) AND YTickLabel
                % AFTER the final ylim has been applied, then set them
                % BOTH so the manual-mode ticks line up with the visible
                % range -- if we only set YTickLabel, MATLAB sticks the
                % labels onto whatever YTick positions existed before
                % the most recent ylim() call (those came from the
                % surface auto-fit, not from [tmin,tmax]).
                allLabels = {};
                for k = 1:numel(axArr)
                    if k > numel(masks), continue; end
                    if isempty(axArr(k).YTickLabel), continue; end
                    allLabels = [allLabels; axArr(k).YTickLabel(:)]; %#ok<AGROW>
                end
                if ~isempty(allLabels)
                    maxLen = max(cellfun(@(s) numel(strtrim(s)), allLabels));
                    for k = 1:numel(axArr)
                        if k > numel(masks), continue; end
                        ytPos = axArr(k).YTick;
                        lbls  = axArr(k).YTickLabel;
                        if isempty(lbls), continue; end
                        for jj = 1:numel(lbls)
                            s = strtrim(lbls{jj});
                            lbls{jj} = [repmat(' ', 1, maxLen - numel(s)), s];
                        end
                        % Set BOTH together so position and label are
                        % paired correctly under manual mode.
                        axArr(k).YTick      = ytPos;
                        axArr(k).YTickLabel = lbls;
                    end
                end
                
                % Always render a colorbar in EVERY cbar cell so the
                % column reservations stay aligned across the three
                % map panels.  When showCbar is false, or for the
                % lower split rows of a "true colorbar" panel, draw a
                % fully-transparent placeholder colorbar so the cell
                % keeps its 50 px width and the data axes stay flush.
                for k = 1:numel(cbarArr)
                    hc = cbarArr(k);
                    cla(hc);
                    hc.Visible = 'off';
                    if showCbar && k == 1
                        hc.CLim = clim;
                        colormap(hc, cmap);
                        cb = colorbar(hc);
                        cb.Label.String = '\DeltaA or Intensity';
                    else
                        % Invisible placeholder colorbar to hold the
                        % column width identically to the visible case.
                        hc.CLim = clim;
                        colormap(hc, cmap);
                        cb = colorbar(hc);
                        cb.Visible = 'off';
                    end
                end
            end
            
            function ticks = chooseNiceTicks(lo, hi, nMax)
                %CHOOSENICETICKS  Pick at most nMax "nice" tick values
                %   in [lo, hi] using a 1/2/5 * 10^k step ladder.
                if hi <= lo || ~isfinite(lo) || ~isfinite(hi)
                    ticks = [];
                    return;
                end
                rng = hi - lo;
                rawStep = rng / max(nMax, 1);
                p = 10^floor(log10(rawStep));
                m = rawStep / p;
                if m <= 1
                    nice = 1;
                elseif m <= 2
                    nice = 2;
                elseif m <= 5
                    nice = 5;
                else
                    nice = 10;
                end
                step = nice * p;
                tLo  = ceil(lo / step) * step;
                tHi  = floor(hi / step) * step;
                if tHi < tLo
                    ticks = [];
                else
                    ticks = tLo:step:tHi;
                end
            end
            
            function onMapClick(ev)
                pt = ev.IntersectionPoint;
                if isempty(pt), return; end
                if isprop(ev, 'Source') && isa(ev.Source, ...
                        'matlab.ui.control.UIAxes')
                    dragGaAx = ev.Source;
                else
                    dragGaAx = [];
                end
                draggingG = true;
                applyGaPoint(pt(1));
            end
            function onGaMotion()
                if ~draggingG, return; end
                if ~app.mouseInsideFigure(f)
                    onGaUp();  return;
                end
                ax = dragGaAx;
                if isempty(ax) || ~isvalid(ax), return; end
                cp = ax.CurrentPoint;
                if isempty(cp) || size(cp,1) < 1, return; end
                applyGaPoint(cp(1, 1));
            end
            function onGaUp()
                draggingG = false;
            end
            function applyGaPoint(wlClick)
                wlR = [min(app.wavelength), max(app.wavelength)];
                curKinWL = min(max(wlClick, wlR(1)), wlR(2));
                edKinWL.Value = curKinWL;
                for kk = 1:numel(hMapVLine)
                    if isgraphics(hMapVLine(kk))
                        hMapVLine(kk).Value = curKinWL;
                    end
                end
                plotKinetics();
                drawnow limitrate;
            end
            
            function plotKinetics()
                if isempty(app.gaResultFit), return; end
                [~, iw] = min(abs(app.wavelength - curKinWL));
                wlActual = app.wavelength(iw);
                
                dataRow = app.deltaA(iw, :);
                fitRow  = app.gaResultFit(iw, :);
                resRow  = dataRow - fitRow;
                
                % --- Left panel: Data + Fit ---
                app.claAll(axKinFit);  hold(axKinFit, 'on');
                plot(axKinFit, app.delay, dataRow, 'k.', ...
                    'MarkerSize', 7, 'DisplayName', 'Data');
                plot(axKinFit, app.delay, fitRow, 'r-', ...
                    'LineWidth', 1.6, 'DisplayName', 'Fit');
                yline(axKinFit, 0, 'k:', 'HandleVisibility', 'off');
                hold(axKinFit, 'off');  grid(axKinFit, 'on');
                xlabel(axKinFit, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axKinFit, '\DeltaA or Intensity');
                legend(axKinFit, 'Location', 'best');
                title(axKinFit, sprintf(['Kinetics at  %s = %.2f %s   ' ...
                    '(requested %.2f %s)'], app.xSym(), wlActual, app.xUnitAx(), ...
                    curKinWL, app.xUnitAx()));
                
                % --- Right panel: Residual only ---
                cla(axKinRes);  hold(axKinRes, 'on');
                plot(axKinRes, app.delay, resRow, '-', ...
                    'Color', [0.45 0.45 0.45], 'LineWidth', 0.9);
                yline(axKinRes, 0, 'k:');
                hold(axKinRes, 'off');  grid(axKinRes, 'on');
                xlabel(axKinRes, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axKinRes, 'Residual (\DeltaA or Intensity)');
                title(axKinRes, sprintf('Fit residual at %s = %.2f %s', ...
                    app.xSym(), wlActual, app.xUnitAx()));
                
                % Apply x-scale to BOTH panels (linked behavior).
                useLog = strcmpi(ddKinScale.Value, 'Log') && ...
                    any(app.delay > 0);
                if useLog
                    set(axKinFit, 'XScale', 'log');
                    set(axKinRes, 'XScale', 'log');
                    xl = [min(app.delay(app.delay > 0)), max(app.delay)];
                else
                    set(axKinFit, 'XScale', 'linear');
                    set(axKinRes, 'XScale', 'linear');
                    xl = [min(app.delay), max(app.delay)];
                end
                % When the fit was restricted to a sub-range of the
                % delay axis, zoom the display onto that range so only
                % the fitted portion of the data is shown.
                fr = getFitTRangeIfRestricted();
                if ~isempty(fr)
                    if useLog, fr(1) = max(fr(1), xl(1)); end
                    if fr(1) < fr(2), xl = fr; end
                end
                xlim(axKinFit, xl);
                xlim(axKinRes, xl);
            end
            
            function fr = getFitTRangeIfRestricted()
                % Returns [tmin tmax] when a fit-time-range narrower
                % than the full delay axis is stored; [] otherwise.
                fr = [];
                if isempty(app.gaFitTMin) || isempty(app.gaFitTMax) || ...
                        ~isfinite(app.gaFitTMin) || ~isfinite(app.gaFitTMax)
                    return;
                end
                lo = min(app.gaFitTMin, app.gaFitTMax);
                hi = max(app.gaFitTMin, app.gaFitTMax);
                dLo = min(app.delay);  dHi = max(app.delay);
                span = dHi - dLo;
                if span <= 0, return; end
                % Treat as restricted only if meaningfully narrower.
                if (lo - dLo) > 1e-9 * span || (dHi - hi) > 1e-9 * span
                    fr = [lo, hi];
                end
            end
            
            function exportDADS()
                if isempty(app.gaResultDADS)
                    uialert(f, 'Run the fit first.', 'No results');  return;
                end
                nc    = size(app.gaResultDADS, 2);
                nMain = numel(app.gaResultTau);
                hdr = cell(1, 1 + nc);
                hdr{1} = app.xHdr();
                for ii = 1:nMain
                    hdr{1 + ii} = sprintf('DADS_tau=%.4g_%s', ...
                        app.gaResultTau(ii), app.tUnitHdr());
                end
                if app.gaResultHasInf, hdr{end} = 'DADS_tau=inf'; end
                
                mat = [app.wavelength(:), app.gaResultDADS];
                % Spectra of the coherent-artifact terms, if any, follow
                % the DADS.
                caR = app.gaResultCA;
                if isstruct(caR) && isfield(caR, 'spec') && ...
                        size(caR.spec, 1) == numel(app.wavelength)
                    caHdr = {'artifact_Gaussian', 'artifact_1st_derivative', ...
                             'artifact_2nd_derivative'};
                    kCA = min(3, size(caR.spec, 2));
                    mat = [mat, caR.spec(:, 1:kCA)];
                    hdr = [hdr, caHdr(1:kCA)];
                end
                
                [file, path] = app.putFile( ...
                    {'*.csv','CSV (*.csv)'; '*.dat','DAT (tab, *.dat)'}, ...
                    'Export DADS', 'DADS.csv');
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    T = array2table(mat, 'VariableNames', hdr);
                    writetable(T, full, 'Delimiter', dlm);
                    uialert(f, sprintf('Saved: %s', full), ...
                        'DADS exported', 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Export error');
                end
            end
            
            function exportEADS()
                if isempty(app.gaResultEADS)
                    uialert(f, ['Spectra not available. Either the fit ' ...
                        'has not been run yet, or the computation ' ...
                        'failed (near-degenerate lifetimes).'], ...
                        'No results');
                    return;
                end
                
                useStretched = any(resStretch());
                if useStretched
                    label = 'SAS';   sciTag = 'Species-associated spectra';
                else
                    label = 'EADS';  sciTag = 'Evolution-associated DS';
                end
                
                nDecay = numel(app.gaResultTau);
                nTot   = size(app.gaResultEADS, 2);
                tauSorted = sort(app.gaResultTau(:), 'ascend');
                
                hdr = cell(1, 1 + nTot);
                hdr{1} = app.xHdr();
                for ii = 1:nDecay
                    hdr{1 + ii} = sprintf('%s_S%d_tau=%.4g_%s', ...
                        label, ii, tauSorted(ii), app.tUnitHdr());
                end
                if app.gaResultHasInf
                    hdr{end} = sprintf('%s_S%d_tau=inf', label, nTot);
                end
                
                mat = [app.wavelength(:), app.gaResultEADS];
                
                [file, path] = app.putFile( ...
                    {'*.csv','CSV (*.csv)'; '*.dat','DAT (tab, *.dat)'}, ...
                    sprintf('Export %s', label), [label '.csv']);
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    T = array2table(mat, 'VariableNames', hdr);
                    writetable(T, full, 'Delimiter', dlm);
                    uialert(f, sprintf( ...
                        ['Saved %s (%s, %d species).\n\nFile: %s'], ...
                        label, sciTag, nTot, full), ...
                        sprintf('%s exported', label), 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Export error');
                end
            end
            
            function exportMatrix2D(which, label, tag)
                if isempty(app.gaResultFit)
                    uialert(f, 'Run the fit first.', 'No results');  return;
                end
                switch which
                    case 'fit'
                        Mdata = app.gaResultFit;
                    case 'residual'
                        Mdata = app.deltaA - app.gaResultFit;
                    otherwise
                        return;
                end
                % Export a (wavelength x delay) matrix in the same layout as
                % input data files: first row = [0, t1..tN], first col =
                % [0, lambda1..lambdaM], interior = Mdata. The resulting file
                % is compatible with Load Data.
                suggested = sprintf('GA_%s.csv', tag);
                [file, path] = app.putFile( ...
                    {'*.csv','CSV comma-delimited (*.csv)'; ...
                     '*.dat','DAT tab-delimited (*.dat)'; ...
                     '*.txt','TXT tab-delimited (*.txt)'}, ...
                    sprintf('Export %s (2D matrix)', label), suggested);
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    app.writeDataFile(full, ...
                        app.wavelength, app.delay, Mdata, dlm);
                    uialert(f, sprintf( ...
                        ['Saved %s matrix  (%d ' app.uiStr('\lambda') ...
                         ' x %d t).\n\n' ...
                         'File: %s'], ...
                        label, numel(app.wavelength), numel(app.delay), full), ...
                        [label ' exported'], 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Export error');
                end
            end
            
            function exportKinetics()
                if isempty(app.gaResultFit)
                    uialert(f, 'Run the fit first.', 'No results');  return;
                end
                % Export the kinetics trace currently displayed in the bottom
                % panel: data / fit / residual vs delay at the selected lambda.
                [~, iw] = min(abs(app.wavelength - curKinWL));
                wlActual = app.wavelength(iw);
                
                dataCol = app.deltaA(iw, :).';
                fitCol  = app.gaResultFit(iw, :).';
                resCol  = dataCol - fitCol;
                tCol    = app.delay(:);
                mat     = [tCol, dataCol, fitCol, resCol];
                hdr     = {sprintf('delay_%s', app.tUnitHdr()), ...
                           'data_dA', 'fit_dA', 'residual_dA'};
                
                suggested = sprintf('kinetics_%.1f%s.csv', wlActual, app.xUnitTxt());
                [file, path] = app.putFile( ...
                    {'*.csv','CSV (*.csv)'; '*.dat','DAT (tab, *.dat)'}, ...
                    sprintf('Export kinetics @ %.2f %s', wlActual, app.xUnitTxt()), suggested);
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    T = array2table(mat, 'VariableNames', hdr);
                    writetable(T, full, 'Delimiter', dlm);
                    uialert(f, sprintf( ...
                        ['Saved kinetics at %s = %.2f %s.\n\n' ...
                         'Columns: delay, data, fit, residual.\n' ...
                         'File: %s'], app.pick(strcmp(app.specUnit,'cm-1'),'nu','lambda'), ...
                         wlActual, app.xUnitTxt(), full), ...
                        'Kinetics exported', 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Export error');
                end
            end
            
            function s = ternary(cond, a, b)
                if cond, s = a; else, s = b; end
            end
        end
        
        function openTargetWindow(app)
            %OPENTARGETWINDOW  Target analysis: fit a compartment model
            %   (kinetic scheme) to the whole data set and show the
            %   species-associated spectra and populations it implies.
            %   The scheme is edited as text or in a K-matrix table and
            %   is shown as a box-and-arrow diagram.
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');
                return;
            end
            if ~isempty(app.targetFig) && isvalid(app.targetFig)
                figure(app.targetFig);  return;
            end

            f = uifigure('Name', 'Target Analysis', ...
                'Position', app.fitFigurePosition(1480, 980));
            app.applyLightTheme(f);
            app.targetFig = f;
            f.CloseRequestFcn = @(~,~) onClose();
            f.WindowButtonMotionFcn = @(~,~) onPickMotion();
            f.WindowButtonUpFcn     = @(~,~) onPickUp();

            gl = uigridlayout(f, [1, 2]);
            gl.ColumnWidth = {470, '1x'};
            gl.Padding = [10 10 10 10];  gl.ColumnSpacing = 10;

            % ---- Presets: generic models first, then examples ----
            presetKeys  = {'seq3', 'seq3gs', 'seq4', 'par3', 'branch', ...
                           'branchLoss', 'rev', 'exISC', 'exET', 'exSF', ...
                           'custom'};
            presetNames = { ...
                'Sequential:  A -> B -> C', ...
                'Sequential, all decay:  A -> B -> C -> GS', ...
                'Sequential, 4 compartments:  A -> B -> C -> D', ...
                'Parallel:  A, B, C decay independently', ...
                'Branching:  A -> B,  A -> C', ...
                'Branching with loss:  A -> B -> GS,  A -> GS', ...
                'Reversible:  A <-> B -> GS', ...
                'Example - intersystem crossing:  S1 -> T1, S1 -> S0', ...
                'Example - energy transfer:  D -> A, D -> GS', ...
                'Example - singlet fission:  S1 -> TT -> 2 T', ...
                'Custom'};
            presetTxt = { ...
                {'A -> B'; 'B -> C'}, ...
                {'A -> B'; 'B -> C'; 'C -> GS'}, ...
                {'A -> B'; 'B -> C'; 'C -> D'}, ...
                {'A -> GS'; 'B -> GS'; 'C -> GS'; 'init: A = 1, B = 1, C = 1'}, ...
                {'A -> B'; 'A -> C'; 'B -> GS'; 'C -> GS'}, ...
                {'A -> B'; 'B -> GS'; 'A -> GS'}, ...
                {'A <-> B : kf, kb'; 'B -> GS'}, ...
                {'S1 -> T1 : kisc'; 'S1 -> S0 : kS'; 'T1 -> S0 : kT'}, ...
                {'D -> A : ket'; 'D -> GS : kD'; 'A -> GS : kA'}, ...
                {'S1 -> TT'; 'TT -> 2 T'}, ...
                {}};
            % Syntax summary: one list for the tooltip of the text box and
            % for the "Syntax..." window (a table there, so that the two
            % columns stay aligned whatever the font and window width).
            syntaxIntro = ['One statement per line, or several on a line ' ...
                'separated by ";".' newline ...
                'Names of compartments and rate parameters are free. ' ...
                'Text after % is a comment.'];
            syntaxRows = { ...
                'A -> B',                   'first-order step with its own rate parameter'; ...
                'A -> B : k1',              '... with a named rate parameter'; ...
                'B -> C : k1',              'the same name = the same parameter (shared)'; ...
                'B -> C : 2*k1',            'rate tied to a parameter (also k1/2)'; ...
                'A -> B -> C',              'a chain is split into single steps'; ...
                'A <-> B : kf, kb',         'reversible step: two rates'; ...
                'B -> 2 C',                 'two C per B (stoichiometric factor)'; ...
                'A -> GS',                  'loss to the ground state; GS, S0 and 0 are the dark sink (no spectrum)'; ...
                'species: A, B, C',         'declare compartments (sets the order; may name one without steps)'; ...
                'init: A = 0.7, B = 0.3',   'initial populations (default: first compartment = 1)'; ...
                'zero: B 700 900',          'spectrum of B = 0 from 700 to 900'; ...
                'equal: C = B',             'spectrum of C equals that of B'; ...
                'equal: C = 0.5*B 600 700', '... proportional, in a range'; ...
                'nonneg',                   'all spectra >= 0'; ...
                'nonneg: A, B',             '... only those of A and B'};
            syntaxNote = ['Ranges are in the units of the spectral axis. ' ...
                'Only first-order steps are supported.'];
            syntaxHelp = [syntaxIntro newline];
            for qh = 1:size(syntaxRows, 1)
                syntaxHelp = [syntaxHelp newline '  ' syntaxRows{qh, 1} ...
                    '    ' syntaxRows{qh, 2}]; %#ok<AGROW>
            end
            syntaxHelp = [syntaxHelp newline newline syntaxNote];

            % ============== Left: scheme and fit setup ==============
            left = uipanel(gl, 'Title', 'Kinetic scheme and fit setup');
            left.Layout.Column = 1;
            lgl = uigridlayout(left, [15, 1]);
            % rows: preset, scheme tabs, apply row, parameter table, t0,
            %       FWHM, IRF model, coherent artifact, artifact width,
            %       fit t range, run row, export label, export row,
            %       status, result text.  15 rows.
            lgl.RowHeight = {30, 200, 30, 140, 28, 28, 28, 28, 28, 28, 34, ...
                             20, 30, 22, '1x'};
            lgl.Padding = [8 8 8 8];  lgl.RowSpacing = 5;

            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Model:', 'FontWeight', 'bold');
            ddPreset = uidropdown(row, 'Items', presetNames, ...
                'ItemsData', presetKeys, 'Value', 'seq3');

            % Three views of the same scheme: text, K matrix, diagram.
            tabs = uitabgroup(lgl);
            tabTxt = uitab(tabs, 'Title', 'Scheme (text)');
            gT = uigridlayout(tabTxt, [1, 1]);
            gT.Padding = [4 4 4 4];
            txtScheme = uitextarea(gT, 'Value', presetTxt{1}, ...
                'FontName', 'Courier New', 'FontSize', 12, ...
                'Tooltip', syntaxHelp);
            tabK = uitab(tabs, 'Title', 'K matrix');
            gK = uigridlayout(tabK, [2, 1]);
            gK.RowHeight = {'1x', 26};
            gK.Padding = [4 4 4 4];  gK.RowSpacing = 4;
            tblK = uitable(gK, 'ColumnEditable', true, ...
                'Tooltip', ['Column = compartment a step starts from, ' ...
                    'row = where it goes. A cell holds the rate of that ' ...
                    'step: a parameter name (k1), a tied rate (2*k1), ' ...
                    'with "[x2]" appended for two products per step. ' ...
                    'Empty = no step. The same name in several cells ' ...
                    'shares the parameter. Last row: initial population.']);
            row = uigridlayout(gK, [1, 4]);
            row.ColumnWidth = {'1x', 'fit', '1x', 'fit'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 5;
            edNewSp = uieditfield(row, 'text', ...
                'Tooltip', 'Name of a new compartment.');
            btnAddSp = uibutton(row, 'Text', 'Add compartment');
            ddRemSp = uidropdown(row, 'Items', {'-'});
            btnRemSp = uibutton(row, 'Text', 'Remove');
            tabDia = uitab(tabs, 'Title', 'Diagram');
            gD = uigridlayout(tabDia, [1, 1]);
            gD.Padding = [4 4 4 4];
            try, gD.BackgroundColor = 'w'; catch, end
            axDia = uiaxes(gD);
            styleDiaAxes(axDia);
            try, axDia.Toolbar.Visible = 'off'; catch, end

            row = uigridlayout(lgl, [1, 3]);
            row.ColumnWidth = {'1x', '1x', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            btnApply = uibutton(row, 'Text', 'Apply scheme', ...
                'Tooltip', ['Read the scheme text and list its rate ' ...
                    'parameters in the table below. Parameters that ' ...
                    'were already listed keep their values. "Run fit" ' ...
                    'does this as well.']);
            btnFromGA = uibutton(row, 'Text', app.uiStr('\tau from GA'), ...
                'Tooltip', ['Copy the time constants of the last Global ' ...
                    'Analysis fit into the table (fastest first, in ' ...
                    'table order), together with its t0 and FWHM.']);
            btnHelp = uibutton(row, 'Text', 'Syntax...');

            tblPar = uitable(lgl, ...
                'ColumnName', app.uiStr({'Rate', ...
                    sprintf('\\tau_{init} (%s)', app.tUnitTxt()), 'Fixed', ...
                    '\tau min', '\tau max', 'Used in'}), ...
                'ColumnFormat', {'char', 'numeric', 'logical', ...
                                 'char', 'char', 'char'}, ...
                'ColumnEditable', [false true true true true false], ...
                'ColumnWidth', {50, 85, 48, 58, 58, 'auto'}, ...
                'RowName', {}, ...
                'Tooltip', app.uiStr(['One row per rate parameter k; the ' ...
                    'table shows its time constant \tau = 1/k. "\tau min" ' ...
                    'and "\tau max" are optional limits (leave empty for ' ...
                    'none).']));

            % IRF defaults: the last GA result if there is one.
            if ~isempty(app.gaResultT0), t0Def = app.gaResultT0; else, t0Def = app.gaT0; end
            if ~isempty(app.gaResultFWHM), fwDef = app.gaResultFWHM; else, fwDef = app.gaFWHM; end
            if isempty(t0Def) || ~isfinite(t0Def), t0Def = 0; end
            if isempty(fwDef) || ~isfinite(fwDef) || fwDef <= 0, fwDef = 0.1; end

            row = uigridlayout(lgl, [1, 3]);
            row.ColumnWidth = {110, '1x', 'fit'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', app.uiStr(sprintf('IRF t_0 (%s):', app.tUnitTxt())));
            edT0 = uieditfield(row, 'numeric', 'Value', t0Def);
            cbT0fix = uicheckbox(row, 'Text', 'Fixed', 'Value', true);

            row = uigridlayout(lgl, [1, 3]);
            row.ColumnWidth = {110, '1x', 'fit'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', sprintf('IRF FWHM (%s):', app.tUnitTxt()));
            edFW = uieditfield(row, 'numeric', 'Limits', [1e-6 Inf], ...
                'Value', max(fwDef, 1e-6));
            cbFWfix = uicheckbox(row, 'Text', 'Fixed', 'Value', true);

            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {110, '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'IRF model:');
            ddIRF = uidropdown(row, ...
                'Items', {'Gaussian IRF', 'None (IRF-free)'}, ...
                'ItemsData', {'irf', 'none'}, 'Value', 'irf', ...
                'Tooltip', ['Gaussian IRF: every decay mode is the ' ...
                    'closed-form exponential convolved with the IRF.' newline ...
                    'None: pure kinetics starting at t0, no convolution; ' ...
                    't0 and FWHM are not fitted. Intended for a Fit t ' ...
                    'range that starts after the IRF.']);

            % Coherent artifact: extra terms at time zero with free
            % spectra (see caBasis). Default: as in the stored global
            % fit, with a width of its own kept at the value found there.
            caDef = app.caOptions([]);
            if isstruct(app.gaResultCA) && isfield(app.gaResultCA, 'order') ...
                    && app.gaResultCA.order > 0
                caDef.order = app.gaResultCA.order;
                if ~strcmp(app.gaResultCA.mode, 'irf')
                    caDef.mode = 'fixed';
                    caDef.fwhm = app.gaResultCA.fwhm;
                end
            end
            if ~(isfinite(caDef.fwhm) && caDef.fwhm > 0), caDef.fwhm = fwDef; end
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {110, '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Coh. artifact:');
            ddCA = uidropdown(row, ...
                'Items',     {'None', '1 term: Gaussian', ...
                              '2 terms: + 1st derivative', ...
                              '3 terms: + 2nd derivative'}, ...
                'ItemsData', {0, 1, 2, 3}, 'Value', caDef.order, ...
                'Tooltip', ['Signal that exists only while pump and ' ...
                    'probe overlap. It is fitted with up to three extra ' ...
                    'terms at time zero, each with a free spectrum: a ' ...
                    'Gaussian and its first and second time derivative. ' ...
                    'The spectral constraints of the scheme do not apply ' ...
                    'to them.' newline ...
                    'With these terms, steps faster than the IRF width ' ...
                    'are less well determined: switch them on only when ' ...
                    'the data show an artifact.'], ...
                'ValueChangedFcn', @(~,~) syncCAEnable());
            row = uigridlayout(lgl, [1, 3]);
            row.ColumnWidth = {110, '1x', 70};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            lblCAW = uilabel(row, 'Text', 'Artifact width:');
            ddCAW = uidropdown(row, ...
                'Items',     {'Same as IRF', 'Own, fitted', 'Own, fixed'}, ...
                'ItemsData', {'irf', 'free', 'fixed'}, 'Value', caDef.mode, ...
                'Tooltip', ['Width (FWHM) of the artifact terms: the IRF ' ...
                    'FWHM, or a width of their own, fitted (within a ' ...
                    'factor of 5 of the IRF FWHM) or fixed at the value ' ...
                    'on the right. With a width of their own the IRF ' ...
                    'FWHM is hardly determined by the fit: keep it fixed.'], ...
                'ValueChangedFcn', @(~,~) syncCAEnable());
            edCAW = uieditfield(row, 'numeric', 'Limits', [1e-6 Inf], ...
                'Value', max(caDef.fwhm, 1e-6), ...
                'Tooltip', sprintf(['Own width of the artifact terms ' ...
                    '(FWHM, %s): starting value or fixed value.'], ...
                    app.tUnitTxt()));

            row = uigridlayout(lgl, [1, 5]);
            row.ColumnWidth = {'fit', 70, 'fit', 70, 'fit'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 4;
            uilabel(row, 'Text', sprintf('Fit t range (%s):', app.tUnitTxt()));
            edFitTmin = uieditfield(row, 'numeric', 'Value', min(app.delay), ...
                'Tooltip', ['Delay points outside this range are left ' ...
                    'out of the fit (the model is still drawn on the ' ...
                    'full axis).']);
            uilabel(row, 'Text', '~', 'HorizontalAlignment', 'center');
            edFitTmax = uieditfield(row, 'numeric', 'Value', max(app.delay), ...
                'Tooltip', 'Upper bound of the fitted delay range.');
            btnFull = uibutton(row, 'Text', 'Full', ...
                'Tooltip', 'Reset to the full delay range.');

            row = uigridlayout(lgl, [1, 3]);
            row.ColumnWidth = {'1.2x', '1x', '0.8x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            btnRun = uibutton(row, 'Text', 'Run fit', ...
                'BackgroundColor', [0.60 0.25 0.45], ...
                'FontColor', 'w', 'FontWeight', 'bold');
            btnUse = uibutton(row, 'Text', app.uiStr('Fitted \rightarrow start'), ...
                'Enable', 'off', ...
                'Tooltip', ['Copy the fitted time constants, t0 and FWHM ' ...
                    'into the starting values, e.g. to refit with ' ...
                    'another choice of fixed parameters.']);
            btnReset = uibutton(row, 'Text', 'Reset');

            uilabel(lgl, 'Text', 'Export:', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 5]);
            row.ColumnWidth = {'1x', '1x', '1.3x', '1.2x', '1.3x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 5;
            btnExpSAS = uibutton(row, 'Text', 'SAS', 'Enable', 'off');
            btnExpDAS = uibutton(row, 'Text', 'DAS', 'Enable', 'off');
            btnExpPop = uibutton(row, 'Text', 'Populations', 'Enable', 'off');
            btnExpFit = uibutton(row, 'Text', 'Fit matrix', 'Enable', 'off');
            btnExpPar = uibutton(row, 'Text', 'Parameters', 'Enable', 'off');

            lblStatus = uilabel(lgl, 'Text', 'Ready.', ...
                'FontColor', [0.2 0.2 0.2]);
            idleText = {'Choose or write a scheme, set the starting time'; ...
                        'constants, then click "Run fit".'};
            txtResult = uitextarea(lgl, 'Editable', 'off', ...
                'Value', idleText, 'FontName', 'Monospaced');

            % ============== Right: results ==============
            right = uipanel(gl, 'Title', 'Results');
            right.Layout.Column = 2;
            rgl = uigridlayout(right, [3, 2]);
            rgl.ColumnWidth = {'1x', '1x'};
            rgl.RowHeight   = {30, '1x', '1x'};
            rgl.Padding     = [8 8 8 8];  rgl.RowSpacing = 8;
            rgl.ColumnSpacing = 10;

            ctl = uigridlayout(rgl, [1, 10]);
            ctl.Layout.Row = 1;  ctl.Layout.Column = [1 2];
            ctl.ColumnWidth = {'fit', 70, 'fit', 75, 'fit', 70, 'fit', ...
                               'fit', 'fit', '1x'};
            ctl.Padding = [0 0 0 0];  ctl.ColumnSpacing = 8;
            wlLo = min(app.wavelength);  wlHi = max(app.wavelength);
            if ~(wlHi > wlLo), wlHi = wlLo + 1; end
            uilabel(ctl, 'Text', 'Spectra:', 'HorizontalAlignment', 'right');
            ddSpec = uidropdown(ctl, 'Items', {'SAS', 'DAS'}, 'Value', 'SAS', ...
                'Tooltip', ['SAS: spectrum of every compartment of the ' ...
                    'scheme. DAS: amplitude spectrum of every decay ' ...
                    'time of the same fit (what a global fit with ' ...
                    'these lifetimes gives).']);
            uilabel(ctl, 'Text', 'Time axis:', 'HorizontalAlignment', 'right');
            ddScale = uidropdown(ctl, 'Items', {'Log', 'Linear'}, 'Value', 'Log');
            uilabel(ctl, 'Text', sprintf('Kinetics at %s (%s):', ...
                app.xSymTxt(), app.xUnitTxt()), 'HorizontalAlignment', 'right');
            edWL = uieditfield(ctl, 'numeric', 'Limits', [wlLo, wlHi], ...
                'Value', min(max(app.selWL, wlLo), wlHi), ...
                'Tooltip', ['Probe position of the kinetics panel. It ' ...
                    'is marked by a line on the spectra panel and on ' ...
                    'the residual map: click or drag on either of ' ...
                    'them to move it.']);
            cbNorm = uicheckbox(ctl, 'Text', 'Normalize', 'Value', false, ...
                'Tooltip', 'Scale every spectrum to a maximum of 1.');
            cbContrib = uicheckbox(ctl, 'Text', 'Contributions', ...
                'Value', true, ...
                'Tooltip', ['Show spectrum x population of every ' ...
                    'compartment under the fitted trace.']);
            btnSVD = uibutton(ctl, 'Text', 'Residual SVD...', ...
                'Enable', 'off', ...
                'Tooltip', ['Singular value decomposition of the ' ...
                    'residual matrix. Structure in the first singular ' ...
                    'vectors is signal the model does not describe; ' ...
                    'pure noise has none.']);
            uilabel(ctl, 'Text', '');

            axSAS = uiaxes(rgl);
            axSAS.Layout.Row = 2;  axSAS.Layout.Column = 1;
            axPop = uiaxes(rgl);
            axPop.Layout.Row = 2;  axPop.Layout.Column = 2;
            axKin = uiaxes(rgl);
            axKin.Layout.Row = 3;  axKin.Layout.Column = 1;
            axRes = uiaxes(rgl);
            axRes.Layout.Row = 3;  axRes.Layout.Column = 2;

            % ---- State ----
            sch = [];      % parsed scheme shown in the tables
            res = [];      % last fit result ([] = none)
            svdFig = [];   % residual-SVD window of this result
            synFig = [];   % syntax summary window
            diaFig = [];   % larger view of the scheme diagram ...
            axDiaWin = []; % ... and its axes
            % Probe-position marker of the kinetics panel: a line on the
            % spectra panel and on the residual map, moved by click / drag.
            hWLspec   = [];      % marker on the spectra panel (vertical)
            hWLres    = [];      % marker on the residual map (horizontal)
            draggingW = false;   % a drag is in progress
            dragAx    = [];      % axes the drag started on
            dragIsY   = false;   % ... with the spectral axis on y
            pickTimer = [];      % ends a drag whose release was missed

            % ---- Wire up ----
            ddPreset.ValueChangedFcn  = @(~,~) safeCb(@onPreset, 'Model');
            txtScheme.ValueChangedFcn = @(~,~) onSchemeEdited();
            tabs.SelectionChangedFcn  = @(~,~) safeCb(@onTabChange, 'Scheme');
            axDia.ButtonDownFcn       = @(~,~) safeCb(@openDiagramWindow, 'Diagram');
            tblK.CellEditCallback     = @(~,ev) safeCb(@() onKEdit(ev), 'K matrix');
            btnAddSp.ButtonPushedFcn  = @(~,~) safeCb(@addCompartment, 'Add compartment');
            btnRemSp.ButtonPushedFcn  = @(~,~) safeCb(@removeCompartment, 'Remove compartment');
            btnApply.ButtonPushedFcn  = @(~,~) safeCb(@applyScheme, 'Apply scheme');
            btnFromGA.ButtonPushedFcn = @(~,~) safeCb(@fromGA, 'Copy from GA');
            btnHelp.ButtonPushedFcn   = @(~,~) safeCb(@showSyntaxHelp, 'Syntax');
            ddIRF.ValueChangedFcn     = @(~,~) onIRFChange();
            btnFull.ButtonPushedFcn   = @(~,~) setFitTFull();
            btnRun.ButtonPushedFcn    = @(~,~) safeCb(@doRun, 'Run fit');
            btnUse.ButtonPushedFcn    = @(~,~) safeCb(@useFitted, 'Fitted to start');
            btnReset.ButtonPushedFcn  = @(~,~) safeCb(@doReset, 'Reset');
            btnExpSAS.ButtonPushedFcn = @(~,~) safeCb(@exportSAS, 'Export SAS');
            btnExpDAS.ButtonPushedFcn = @(~,~) safeCb(@exportDAS, 'Export DAS');
            btnExpPop.ButtonPushedFcn = @(~,~) safeCb(@exportPop, 'Export populations');
            btnExpFit.ButtonPushedFcn = @(~,~) safeCb(@exportFit, 'Export fit');
            btnExpPar.ButtonPushedFcn = @(~,~) safeCb(@exportPar, 'Export parameters');
            btnSVD.ButtonPushedFcn    = @(~,~) safeCb(@showResidualSVD, 'Residual SVD');
            edWL.ValueChangedFcn      = @(~,~) safeCb(@onWLChange, 'Kinetics');
            ddScale.ValueChangedFcn   = @(~,~) safeCb(@plotTimePanels, 'Time axis');
            ddSpec.ValueChangedFcn    = @(~,~) safeCb(@plotSpec, 'Spectra');
            cbNorm.ValueChangedFcn    = @(~,~) safeCb(@plotSpec, 'Spectra');
            cbContrib.ValueChangedFcn = @(~,~) safeCb(@plotKin, 'Kinetics');

            % ---- Restore the previous settings ----
            if isstruct(app.tgSet) && isfield(app.tgSet, 'schemeTxt')
                try
                    st = app.tgSet;
                    if any(strcmp(presetKeys, st.preset))
                        ddPreset.Value = st.preset;
                    else
                        ddPreset.Value = 'custom';
                    end
                    txtScheme.Value = st.schemeTxt;
                    if iscell(st.table) && size(st.table, 2) == 6
                        tblPar.Data = st.table;
                    end
                    edT0.Value      = st.t0;
                    edFW.Value      = max(st.fwhm, 1e-6);
                    cbT0fix.Value   = st.t0Fixed;
                    cbFWfix.Value   = st.fwhmFixed;
                    ddIRF.Value     = st.irf;
                    if isfield(st, 'ca')
                        caSt = app.caOptions(st.ca);
                        ddCA.Value  = caSt.order;
                        ddCAW.Value = caSt.mode;
                        if isfinite(caSt.fwhm) && caSt.fwhm > 0
                            edCAW.Value = caSt.fwhm;
                        end
                    end
                    [edFitTmin.Value, edFitTmax.Value] = ...
                        app.clampFitRange(st.fitTMin, st.fitTMax);
                    ddScale.Value   = st.scale;
                    cbNorm.Value    = st.norm;
                    cbContrib.Value = st.contrib;
                    if isfield(st, 'spec') && any(strcmp({'SAS', 'DAS'}, st.spec))
                        ddSpec.Value = st.spec;
                    end
                catch
                    lblStatus.Text = 'Stored settings could not be restored.';
                end
            end
            onIRFChange();
            try
                applyScheme();
            catch ME0
                lblStatus.Text = ['Scheme: ' ME0.message];
            end

            % ---- Restore the previous result (if it fits the data) ----
            if isstruct(app.tgSaved) && isfield(app.tgSaved, 'res') && ...
                    isstruct(app.tgSaved.res) && ...
                    isfield(app.tgSaved.res, 'DAS') && ...
                    isfield(app.tgSaved.res.sch, 'params') && ...
                    size(app.tgSaved.res.SAS, 1) == numel(app.wavelength) && ...
                    size(app.tgSaved.res.conc, 1) == numel(app.delay)
                try
                    res = app.tgSaved.res;
                    edWL.Value = min(max(app.tgSaved.curWL, wlLo), wlHi);
                    txtResult.Value = app.tgSaved.lines;
                    lblStatus.Text  = sprintf('%s  [restored]', ...
                        app.tgSaved.status);
                    setResultButtons('on');
                catch
                    res = [];
                    lblStatus.Text = 'Stored result could not be restored -- run again.';
                end
            end
            try
                plotAll();
            catch
                lblStatus.Text = 'Stored result could not be drawn -- run again.';
            end

            % ============== Nested callbacks ==============
            function safeCb(fcn, label)
                try
                    fcn();
                catch ME
                    uialert(f, sprintf('"%s" failed:\n%s\n\n(In %s line %d)', ...
                        label, ME.message, ...
                        ME.stack(1).name, ME.stack(1).line), ...
                        'Callback error');
                    lblStatus.Text = sprintf('Error in %s.', label);
                end
            end

            function onPreset()
                k = find(strcmp(presetKeys, ddPreset.Value), 1);
                if isempty(k) || isempty(presetTxt{k}), return; end
                txtScheme.Value = presetTxt{k};
                applyScheme();
            end

            function onSchemeEdited()
                % Hand-edited text no longer is the model shown.
                ddPreset.Value = 'custom';
            end

            function onTabChange()
                % The K matrix and the diagram always show the scheme
                % of the text.
                if tabs.SelectedTab == tabK || tabs.SelectedTab == tabDia
                    applyScheme();
                end
            end

            function applyScheme()
                % Parse the text and (re)build the parameter table and
                % the K matrix. Parameters that were already listed keep
                % their time constant, Fixed flag and limits.
                newSch = app.parseTargetScheme(txtScheme.Value);
                old = tblPar.Data;
                nPar = newSch.nP;
                spanA = max(app.delay) - min(app.delay);
                if ~(spanA > 0), spanA = 1; end
                tLo = max(2 * edFW.Value, 1e-3 * spanA);
                tHi = max(spanA / 3, 2 * tLo);
                defTau = logspace(log10(tLo), log10(tHi), max(nPar, 2));
                dat = cell(nPar, 6);
                for q = 1:nPar
                    pn = newSch.params{q};
                    use = {};
                    for w = 1:numel(newSch.steps)
                        if newSch.steps(w).par ~= q, continue; end
                        if newSch.steps(w).mult == 1
                            use{end+1} = newSch.steps(w).label; %#ok<AGROW>
                        else
                            use{end+1} = sprintf('%s (%s)', ...
                                newSch.steps(w).label, newSch.steps(w).expr); %#ok<AGROW>
                        end
                    end
                    dat{q, 1} = pn;
                    dat{q, 2} = round(defTau(q), 3, 'significant');
                    dat{q, 3} = false;
                    dat{q, 4} = '';
                    dat{q, 5} = '';
                    dat{q, 6} = strjoin(use, ';  ');
                    if iscell(old) && size(old, 2) == 6 && ~isempty(old)
                        hit = find(strcmp(old(:, 1), pn), 1);
                        if ~isempty(hit)
                            dat{q, 2} = old{hit, 2};
                            dat{q, 3} = logical(old{hit, 3});
                            dat{q, 4} = old{hit, 4};
                            dat{q, 5} = old{hit, 5};
                        end
                    end
                end
                tblPar.Data = dat;
                sch = newSch;
                refreshK();
                lblStatus.Text = sprintf(['Scheme: %d compartment(s) (%s), ' ...
                    '%d step(s), %d rate parameter(s).'], sch.nS, ...
                    strjoin(sch.species, ', '), numel(sch.steps), sch.nP);
                drawDiagrams();
            end

            function refreshK()
                % K-matrix view of the current scheme: column = from,
                % row = to; last rows: loss to the ground state and the
                % initial populations.
                nC = sch.nS;
                dat = repmat({''}, nC + 2, nC);
                for w = 1:numel(sch.steps)
                    sp = sch.steps(w);
                    if sp.to == 0, r = nC + 1; else, r = sp.to; end
                    if sp.coef == 1
                        dat{r, sp.from} = sp.expr;
                    else
                        dat{r, sp.from} = sprintf('%s [x%g]', sp.expr, sp.coef);
                    end
                end
                for q = 1:nC
                    dat{nC + 2, q} = sprintf('%g', sch.c0(q));
                end
                tblK.Data = dat;
                tblK.ColumnName = cellfun(@(n) ['from ' n], sch.species, ...
                    'UniformOutput', false);
                tblK.RowName = [cellfun(@(n) ['to ' n], sch.species(:), ...
                    'UniformOutput', false); {['to ' sch.sinkName]}; {'initial'}];
                tblK.ColumnFormat = repmat({'char'}, 1, nC);
                tblK.ColumnEditable = true(1, nC);
                ddRemSp.Items = sch.species;
                ddRemSp.Value = sch.species{end};
            end

            function useSchemeStruct(s2)
                % Adopt an edited scheme: write it as text and re-read it,
                % so that the text stays the single description. If the
                % text is not valid, nothing changes.
                newTxt = app.targetSchemeText(s2);
                app.parseTargetScheme(newTxt);       % throws if invalid
                txtScheme.Value = newTxt;
                ddPreset.Value  = 'custom';
                applyScheme();
            end

            function onKEdit(ev)
                % One cell of the K matrix was edited.
                try
                    r = ev.Indices(1);  c = ev.Indices(2);
                    val = strtrim(char(ev.NewData));
                    nC = sch.nS;
                    s2 = sch;
                    if r == nC + 2
                        % initial population
                        if isempty(val), v0 = 0; else, v0 = str2double(val); end
                        if ~isfinite(v0) || v0 < 0
                            error('The initial population must be a number >= 0.');
                        end
                        s2.c0(c) = v0;
                        if ~any(s2.c0 > 0)
                            error('At least one compartment needs a population > 0.');
                        end
                    elseif r == c
                        error('A compartment has no step to itself.');
                    else
                        if r == nC + 1, iTo = 0; else, iTo = r; end
                        hit = find([s2.steps.from] == c & [s2.steps.to] == iTo, 1);
                        if isempty(val)
                            if ~isempty(hit), s2.steps(hit) = []; end
                        else
                            cf = 1;  ex = val;
                            tk = regexp(val, ...
                                '^(.*?)\s*\[\s*x\s*([0-9]*\.?[0-9]+)\s*\]$', ...
                                'tokens', 'once');
                            if ~isempty(tk)
                                ex = strtrim(tk{1});
                                cf = str2double(tk{2});
                            end
                            if iTo == 0 && cf ~= 1
                                error('A loss to the ground state takes no "[x..]" factor.');
                            end
                            if ~isfinite(cf) || cf <= 0
                                error('The factor in "[x..]" must be positive.');
                            end
                            stp = struct('from', c, 'to', iTo, 'coef', cf, ...
                                'par', 0, 'mult', 1, 'label', '', 'expr', ex);
                            if isempty(hit)
                                s2.steps(end+1) = stp;
                            else
                                s2.steps(hit) = stp;
                            end
                        end
                    end
                    useSchemeStruct(s2);
                catch ME
                    refreshK();          % put the table back
                    rethrow(ME);
                end
            end

            function addCompartment()
                nm = strtrim(edNewSp.Value);
                if isempty(nm)
                    error('Type a name for the new compartment first.');
                end
                if any(strcmp(sch.species, nm))
                    error('A compartment called "%s" already exists.', nm);
                end
                s2 = sch;
                s2.species{end+1} = nm;
                s2.c0(end+1, 1)   = 0;
                s2.nonneg(end+1)  = all(sch.nonneg);
                useSchemeStruct(s2);
                edNewSp.Value = '';
                lblStatus.Text = sprintf(['Compartment %s added: give it ' ...
                    'steps in the K matrix.'], nm);
            end

            function removeCompartment()
                k = find(strcmp(sch.species, ddRemSp.Value), 1);
                if isempty(k), return; end
                if sch.nS < 2
                    error('The last compartment cannot be removed.');
                end
                nm = sch.species{k};
                s2 = sch;
                keep = [s2.steps.from] ~= k & [s2.steps.to] ~= k;
                s2.steps = s2.steps(keep);
                for w = 1:numel(s2.steps)
                    if s2.steps(w).from > k, s2.steps(w).from = s2.steps(w).from - 1; end
                    if s2.steps(w).to   > k, s2.steps(w).to   = s2.steps(w).to   - 1; end
                end
                s2.species(k) = [];
                s2.c0(k)      = [];
                s2.nonneg(k)  = [];
                if ~any(s2.c0 > 0), s2.c0(1) = 1; end
                if ~isempty(s2.zero)
                    s2.zero = s2.zero([s2.zero.sp] ~= k);
                    for w = 1:numel(s2.zero)
                        if s2.zero(w).sp > k, s2.zero(w).sp = s2.zero(w).sp - 1; end
                    end
                end
                if ~isempty(s2.rel)
                    s2.rel = s2.rel([s2.rel.dep] ~= k & [s2.rel.master] ~= k);
                    for w = 1:numel(s2.rel)
                        if s2.rel(w).dep    > k, s2.rel(w).dep    = s2.rel(w).dep    - 1; end
                        if s2.rel(w).master > k, s2.rel(w).master = s2.rel(w).master - 1; end
                    end
                end
                useSchemeStruct(s2);
                lblStatus.Text = sprintf(['Compartment %s removed, with its ' ...
                    'steps and constraints.'], nm);
            end

            function [tauV, fixV, loV, hiV] = readTable()
                dat = tblPar.Data;
                nRow = size(dat, 1);
                tauV = zeros(1, nRow);
                fixV = false(1, nRow);
                loV  = nan(1, nRow);
                hiV  = nan(1, nRow);
                for q = 1:nRow
                    v = dat{q, 2};
                    if ischar(v) || isstring(v), v = str2double(v); end
                    if isempty(v) || ~isscalar(v) || ~isfinite(v) || v <= 0
                        error(['Parameter %s: the time constant must be ' ...
                            'a positive number.'], dat{q, 1});
                    end
                    tauV(q) = double(v);
                    fixV(q) = logical(dat{q, 3});
                    loV(q)  = readLimit(dat{q, 4}, dat{q, 1});
                    hiV(q)  = readLimit(dat{q, 5}, dat{q, 1});
                end
            end

            function v = readLimit(c, pname)
                % '' -> NaN (no limit); otherwise a positive number.
                v = NaN;
                if isnumeric(c)
                    if isempty(c) || isnan(c), return; end
                    v = double(c);
                else
                    c = strtrim(char(c));
                    if isempty(c), return; end
                    v = str2double(c);
                end
                if ~isfinite(v) || v <= 0
                    error(['Parameter %s: a limit must be a positive ' ...
                        'number (or empty for none).'], pname);
                end
            end

            function fromGA()
                if isempty(app.gaResultTau)
                    uialert(f, 'Run a Global Analysis fit first.', ...
                        'No GA result');
                    return;
                end
                dat = tblPar.Data;
                tg = sort(app.gaResultTau(:).', 'ascend');
                nCopy = min(size(dat, 1), numel(tg));
                for q = 1:nCopy
                    dat{q, 2} = round(tg(q), 4, 'significant');
                end
                tblPar.Data = dat;
                if ~isempty(app.gaResultT0) && isfinite(app.gaResultT0)
                    edT0.Value = app.gaResultT0;
                end
                if ~isempty(app.gaResultFWHM) && app.gaResultFWHM > 0
                    edFW.Value = max(app.gaResultFWHM, 1e-6);
                end
                lblStatus.Text = sprintf(['%d time constant(s), t0 and FWHM ' ...
                    'copied from the GA result.'], nCopy);
            end

            function useFitted()
                if isempty(res), return; end
                dat = tblPar.Data;
                for q = 1:size(dat, 1)
                    hit = find(strcmp(res.sch.params, dat{q, 1}), 1);
                    if ~isempty(hit)
                        dat{q, 2} = round(res.tau(hit), 5, 'significant');
                    end
                end
                tblPar.Data = dat;
                edT0.Value = res.t0;
                if res.fwhm > 0, edFW.Value = max(res.fwhm, 1e-6); end
                if isfield(res, 'ca') && isstruct(res.ca) && ...
                        ~strcmp(res.ca.mode, 'irf') && res.ca.fwhm > 0
                    edCAW.Value = res.ca.fwhm;
                end
                lblStatus.Text = 'Fitted values copied to the starting values.';
            end

            function onIRFChange()
                % Without IRF, t0 is only the onset of the kinetics and
                % FWHM has no meaning: neither is fitted.
                if strcmp(ddIRF.Value, 'none')
                    edFW.Enable    = 'off';
                    cbFWfix.Enable = 'off';
                    cbT0fix.Value  = true;
                    cbT0fix.Enable = 'off';
                else
                    edFW.Enable    = 'on';
                    cbFWfix.Enable = 'on';
                    cbT0fix.Enable = 'on';
                end
                syncCAEnable();
            end

            function syncCAEnable()
                % The artifact terms need the IRF; their own width is
                % used only when asked for.
                hasIRF = ~strcmp(ddIRF.Value, 'none');
                useCA  = hasIRF && ddCA.Value > 0;
                ddCA.Enable   = app.pick(hasIRF, 'on', 'off');
                ddCAW.Enable  = app.pick(useCA, 'on', 'off');
                lblCAW.Enable = app.pick(useCA, 'on', 'off');
                edCAW.Enable  = app.pick(useCA && ...
                    ~strcmp(ddCAW.Value, 'irf'), 'on', 'off');
            end

            function setFitTFull()
                edFitTmin.Value = min(app.delay);
                edFitTmax.Value = max(app.delay);
            end

            function setResultButtons(state)
                btnUse.Enable    = state;
                btnExpSAS.Enable = state;
                btnExpDAS.Enable = state;
                btnExpPop.Enable = state;
                btnExpFit.Enable = state;
                btnExpPar.Enable = state;
                btnSVD.Enable    = state;
            end

            function doRun()
                applyScheme();
                [tauV, fixV, loV, hiV] = readTable();
                irfM = ddIRF.Value;
                [rLo, rHi] = app.fitRangeToStore(edFitTmin.Value, edFitTmax.Value);
                tR = [rLo, rHi];
                caS = struct('order', ddCA.Value, 'mode', ddCAW.Value, ...
                    'fwhm', edCAW.Value);
                lblStatus.Text = 'Fitting...';  drawnow;
                try
                    [tauO, t0O, fwO, SAS, fitM, info] = app.fitTargetAnalysis( ...
                        app.deltaA, app.wavelength, app.delay, sch, ...
                        tauV, fixV, edT0.Value, edFW.Value, ...
                        cbT0fix.Value, cbFWfix.Value, irfM, tR, loV, hiV, caS);
                catch ME
                    uialert(f, sprintf('Target fit failed:\n%s', ME.message), ...
                        'Fit error');
                    lblStatus.Text = 'Fit failed.';
                    return;
                end
                r = struct();
                r.sch        = sch;
                r.tau        = tauO;
                r.tauFixed   = fixV;
                r.tauStd     = info.tauStd;
                r.atLimit    = info.atLimit;
                r.t0         = t0O;
                r.fwhm       = fwO;
                r.t0Std      = info.t0Std;
                r.fwhmStd    = info.fwhmStd;
                r.fwhmLo      = info.fwhmLo;
                r.fwhmAtLimit = info.fwhmAtLimit;
                r.irfRetry    = info.irfRetry;
                r.t0Fixed    = cbT0fix.Value || strcmp(irfM, 'none');
                r.fwhmFixed  = cbFWfix.Value || strcmp(irfM, 'none');
                r.irf        = irfM;
                r.tRange     = tR;
                r.SAS        = SAS;
                r.DAS        = info.DAS;
                r.conc       = info.conc;
                r.fit        = fitM;
                r.rms        = info.rms;
                r.initialRMS = info.initialRMS;
                r.iters      = info.iters;
                r.converged  = info.converged;
                r.obsTau     = info.obsTau;
                r.nObs       = info.nObs;
                r.nFree      = info.nFree;
                r.reordered  = info.reordered;
                % coherent-artifact terms of this fit ([] = none)
                r.ca = [];
                if info.caOrder > 0
                    r.ca = struct('order', info.caOrder, ...
                        'spec', info.caSpec, 'fwhm', info.caFwhm, ...
                        'fwhmStd', info.caFwhmStd, 'mode', info.caMode);
                end
                res = r;
                closeSVD();

                txtResult.Value = resultLines();
                if res.converged
                    lblStatus.Text = sprintf('Done. RMS = %.4g', res.rms);
                else
                    lblStatus.Text = sprintf( ...
                        'Stopped at the iteration limit. RMS = %.4g', res.rms);
                end
                setResultButtons('on');
                plotAll();

                % Persist so the analysis survives closing this window.
                saveSettings();
                app.tgSaved = struct('res', res, 'curWL', edWL.Value, ...
                    'status', lblStatus.Text, 'lines', {txtResult.Value});
            end

            function L = resultLines()
                u = app.tUnitTxt();
                st1 = res.sch.steps;
                nSt = numel(st1);
                % time constant of every step: tau_par / mult
                tauStep = zeros(1, nSt);
                for q = 1:nSt
                    tauStep(q) = res.tau(st1(q).par) / st1(q).mult;
                end
                L = {'Scheme:'};
                wLab = 1;
                if nSt > 0, wLab = max(cellfun(@numel, {st1.label})); end
                for q = 1:nSt
                    lab = st1(q).label;
                    lab = [lab, repmat(' ', 1, wLab - numel(lab))]; %#ok<AGROW>
                    L{end+1} = sprintf('  %s  : %s', lab, st1(q).expr); %#ok<AGROW>
                end
                if nSt == 0
                    L{end+1} = '  (no steps)';
                end
                c0s = {};
                for q = 1:res.sch.nS
                    if res.sch.c0(q) ~= 0
                        c0s{end+1} = sprintf('%s = %g', ...
                            res.sch.species{q}, res.sch.c0(q)); %#ok<AGROW>
                    end
                end
                L{end+1} = ['Initial population: ' strjoin(c0s, ', ')];
                for q = 1:numel(res.sch.zero)
                    z = res.sch.zero(q);
                    nIn = nnz(app.wavelength >= z.lo & app.wavelength <= z.hi);
                    L{end+1} = sprintf(['Spectrum of %s = 0 for %g .. %g %s ' ...
                        '(%d points)'], res.sch.species{z.sp}, z.lo, z.hi, ...
                        app.xUnitTxt(), nIn); %#ok<AGROW>
                end
                for q = 1:numel(res.sch.rel)
                    rl = res.sch.rel(q);
                    if rl.ratio == 1
                        rhs = res.sch.species{rl.master};
                    else
                        rhs = sprintf('%g x %s', rl.ratio, res.sch.species{rl.master});
                    end
                    if isfinite(rl.lo) && isfinite(rl.hi)
                        L{end+1} = sprintf('Spectrum of %s = %s for %g .. %g %s', ...
                            res.sch.species{rl.dep}, rhs, rl.lo, rl.hi, ...
                            app.xUnitTxt()); %#ok<AGROW>
                    else
                        L{end+1} = sprintf('Spectrum of %s = %s', ...
                            res.sch.species{rl.dep}, rhs); %#ok<AGROW>
                    end
                end
                if any(res.sch.nonneg)
                    L{end+1} = ['Non-negative spectra: ' ...
                        strjoin(res.sch.species(res.sch.nonneg), ', ')];
                end
                L{end+1} = '';
                L{end+1} = sprintf('RMS = %.4g  (start: %.4g), %d iterations', ...
                    res.rms, res.initialRMS, res.iters);
                if ~res.converged
                    L{end+1} = '** NOT CONVERGED: iteration limit reached.';
                    L{end+1} = '   Try other starting values or fix more.';
                end
                if ~isempty(res.tRange)
                    L{end+1} = sprintf('Fit t range: %.4g .. %.4g %s', ...
                        res.tRange(1), res.tRange(2), u);
                end
                L{end+1} = '';
                L{end+1} = 'Rate parameters (tau = 1/k):';
                wPar = 1;
                if res.sch.nP > 0, wPar = max(cellfun(@numel, res.sch.params)); end
                anyUndet = false;
                for q = 1:res.sch.nP
                    pn = res.sch.params{q};
                    pn = [pn, repmat(' ', 1, wPar - numel(pn))]; %#ok<AGROW>
                    sd = res.tauStd(q);
                    kTxt = sprintf('k = %.4g /%s', 1 / res.tau(q), u);
                    if res.tauFixed(q)
                        L{end+1} = sprintf('  %s  tau = %.4g %s (fixed)   %s', ...
                            pn, res.tau(q), u, kTxt); %#ok<AGROW>
                    elseif res.atLimit(q)
                        L{end+1} = sprintf('  %s  tau = %.4g %s  [at its limit]', ...
                            pn, res.tau(q), u); %#ok<AGROW>
                    elseif isinf(sd)
                        anyUndet = true;
                        L{end+1} = sprintf('  %s  tau = %.4g %s  [not determined]', ...
                            pn, res.tau(q), u); %#ok<AGROW>
                    elseif isfinite(sd)
                        L{end+1} = sprintf('  %s  tau = %.4g +/- %.2g %s   %s', ...
                            pn, res.tau(q), sd, u, kTxt); %#ok<AGROW>
                    else
                        L{end+1} = sprintf('  %s  tau = %.4g %s   %s', ...
                            pn, res.tau(q), u, kTxt); %#ok<AGROW>
                    end
                end
                % Steps whose rate is a multiple of a parameter.
                tied = find([st1.mult] ~= 1);
                if ~isempty(tied)
                    L{end+1} = 'Steps with a tied rate:';
                    for q = tied
                        L{end+1} = sprintf('  %s  tau = %.4g %s  (%s)', ...
                            st1(q).label, tauStep(q), u, st1(q).expr); %#ok<AGROW>
                    end
                end
                L{end+1} = '';
                if strcmp(res.irf, 'none')
                    L{end+1} = sprintf('No IRF; kinetics start at t0 = %.4g %s', ...
                        res.t0, u);
                else
                    L{end+1} = irfLine('IRF t0  ', res.t0, res.t0Std, res.t0Fixed, u);
                    if isfield(res, 'fwhmAtLimit') && isequal(res.fwhmAtLimit, true)
                        L{end+1} = sprintf('IRF FWHM = %.4g %s  [not resolved]', ...
                            res.fwhm, u);
                    else
                        L{end+1} = irfLine('IRF FWHM', res.fwhm, res.fwhmStd, res.fwhmFixed, u);
                    end
                    L = [L, app.irfLimitNote(res)];
                end
                if isfield(res, 'ca') && isstruct(res.ca) && res.ca.order > 0
                    caNm = {'Gaussian', '1st derivative', '2nd derivative'};
                    L{end+1} = '';
                    L{end+1} = sprintf('Coherent artifact: %d term(s)', ...
                        res.ca.order);
                    if strcmp(res.ca.mode, 'irf')
                        L{end+1} = '  width = IRF FWHM';
                    else
                        L{end+1} = ['  ' irfLine('width', res.ca.fwhm, ...
                            res.ca.fwhmStd, strcmp(res.ca.mode, 'fixed'), u)];
                    end
                    for q = 1:res.ca.order
                        [~, iMx] = max(abs(res.ca.spec(:, q)));
                        L{end+1} = sprintf('  %-14s largest %.3g at %.4g %s', ...
                            caNm{q}, res.ca.spec(iMx, q), ...
                            app.wavelength(iMx), app.xUnitTxt()); %#ok<AGROW>
                    end
                    if ~strcmp(res.ca.mode, 'irf') && ~res.fwhmFixed && ...
                            res.ca.order >= 2
                        L{end+1} = '  With a width of their own the derivative';
                        L{end+1} = '  terms absorb a change of the IRF FWHM: it';
                        L{end+1} = '  is then hardly determined. Fix it, or use';
                        L{end+1} = '  "Same as IRF".';
                    end
                end

                % Lifetimes the data can actually show (eigenvalues of
                % the rate matrix) -- the time constants of the DAS.
                ot = res.obsTau(isfinite(res.obsTau));
                if ~isempty(ot)
                    L{end+1} = '';
                    L{end+1} = sprintf('Observable lifetimes (%s):', u);
                    L{end+1} = ['  ' strjoin(arrayfun(@(x) sprintf('%.4g', x), ...
                        ot, 'UniformOutput', false), ', ')];
                end

                % Branching fractions of compartments with several ways out.
                first = true;
                for q = 1:res.sch.nS
                    outI = find([st1.from] == q);
                    if numel(outI) < 2, continue; end
                    kOut = 1 ./ tauStep(outI);
                    fr = 100 * kOut / sum(kOut);
                    parts = cell(1, numel(outI));
                    for w = 1:numel(outI)
                        if st1(outI(w)).to == 0
                            dst = 'ground state';
                        else
                            dst = res.sch.species{st1(outI(w)).to};
                        end
                        parts{w} = sprintf('-> %s %.1f %%', dst, fr(w));
                    end
                    if first
                        L{end+1} = ''; %#ok<AGROW>
                        L{end+1} = 'Branching:'; %#ok<AGROW>
                        first = false;
                    end
                    pIdx = [st1(outI).par];
                    note = '';
                    if any(isinf(res.tauStd(pIdx)))
                        note = '   [not determined]';
                    elseif all(res.tauFixed(pIdx)) || numel(unique(pIdx)) == 1
                        note = '   (fixed by the model)';
                    end
                    L{end+1} = sprintf('  %s: %s%s', res.sch.species{q}, ...
                        strjoin(parts, ', '), note); %#ok<AGROW>
                end

                L{end+1} = '';
                L{end+1} = 'Notes:';
                if anyUndet
                    L{end+1} = sprintf(['- %d free rate parameter(s), %d ' ...
                        'observable lifetime(s):'], res.nFree, res.nObs);
                    L{end+1} = '  the marked ones only act in combination and';
                    L{end+1} = '  cannot be separated by these data. Fix one,';
                    L{end+1} = '  tie it to another (k2 = 2*k1), or add a';
                    L{end+1} = '  spectral constraint (zero: / equal:).';
                end
                if any(res.atLimit)
                    L{end+1} = '- A parameter ended on one of its limits: the';
                    L{end+1} = '  data would move it further. No error is';
                    L{end+1} = '  given for it.';
                end
                if res.reordered
                    L{end+1} = '- The same fit is obtained with time constants';
                    L{end+1} = '  interchanged; the order of your starting';
                    L{end+1} = '  values was kept.';
                end
                L{end+1} = '- The fit tests the decay rates, not the scheme:';
                L{end+1} = '  other schemes with the same rates fit equally';
                L{end+1} = '  well. The spectra are those implied by the';
                L{end+1} = '  scheme chosen here.';
                L = L(:);
            end

            function s = irfLine(name, val, sd, isFixed, u)
                if isFixed
                    s = sprintf('%s = %.4g %s (fixed)', name, val, u);
                elseif isinf(sd)
                    s = sprintf('%s = %.4g %s  [not determined]', name, val, u);
                elseif isfinite(sd)
                    s = sprintf('%s = %.4g +/- %.2g %s', name, val, sd, u);
                else
                    s = sprintf('%s = %.4g %s', name, val, u);
                end
            end

            function plotAll()
                plotSpec();
                plotPop();
                plotKin();
                plotRes();
                drawDiagrams();
            end

            function plotTimePanels()
                plotPop();
                plotKin();
                plotRes();
            end

            function applyTimeAxis(ax)
                tA = app.delay(:);
                isLog = strcmp(ddScale.Value, 'Log') && any(tA > 0);
                if isLog
                    ax.XScale = 'log';
                    xl = [min(tA(tA > 0)), max(tA)];
                else
                    ax.XScale = 'linear';
                    xl = [min(tA), max(tA)];
                end
                % A restricted fit range also is the displayed range.
                if ~isempty(res) && ~isempty(res.tRange)
                    fr = res.tRange;
                    if isLog, fr(1) = max(fr(1), xl(1)); end
                    if fr(1) < fr(2), xl = fr; end
                end
                if xl(2) > xl(1), xlim(ax, xl); end
            end

            function plotSpec()
                % Top-left panel: SAS (one per compartment) or DAS (one
                % per decay time) of the same fit.
                app.claAll(axSAS);
                xlabel(axSAS, app.xAxLabel());
                grid(axSAS, 'on');
                isDAS = strcmp(ddSpec.Value, 'DAS');
                if isempty(res)
                    ylabel(axSAS, 'Amplitude (\DeltaA or Intensity)');
                    title(axSAS, 'Spectra (run fit)');
                    legend(axSAS, 'off');
                    xlim(axSAS, [wlLo, wlHi]);
                    hWLspec = wlMarker(axSAS, false);
                    return;
                end
                if isDAS
                    Y = res.DAS;
                    nm = cell(1, size(Y, 2));
                    for q = 1:size(Y, 2)
                        if isfinite(res.obsTau(q))
                            nm{q} = sprintf('\\tau = %.3g %s', ...
                                res.obsTau(q), app.tUnitAx());
                        else
                            nm{q} = '\tau = \infty';
                        end
                    end
                else
                    Y = res.SAS;
                    nm = res.sch.species;
                end
                cols = lines(max(size(Y, 2), 3));
                hold(axSAS, 'on');
                for q = 1:size(Y, 2)
                    y = Y(:, q);
                    if cbNorm.Value
                        m = max(abs(y), [], 'omitnan');
                        if ~isempty(m) && isfinite(m) && m > 0, y = y / m; end
                    end
                    % (not pickable: a click on a curve has to reach the
                    % axes, which moves the probe marker)
                    plot(axSAS, app.wavelength, y, '-', 'LineWidth', 1.4, ...
                        'Color', cols(q, :), 'DisplayName', nm{q}, ...
                        'HitTest', 'off', 'PickableParts', 'none');
                end
                % Spectra of the coherent-artifact terms (dashed): part
                % of the fit, not of the SAS / DAS.
                if isfield(res, 'ca') && isstruct(res.ca) && ...
                        size(res.ca.spec, 1) == numel(app.wavelength)
                    caNm  = {'artifact: Gaussian', 'artifact: 1st deriv.', ...
                             'artifact: 2nd deriv.'};
                    caCol = [0.30 0.30 0.30; 0.50 0.50 0.50; 0.68 0.68 0.68];
                    for q = 1:min(3, size(res.ca.spec, 2))
                        y = res.ca.spec(:, q);
                        if cbNorm.Value
                            m = max(abs(y), [], 'omitnan');
                            if ~isempty(m) && isfinite(m) && m > 0, y = y / m; end
                        end
                        plot(axSAS, app.wavelength, y, '--', 'LineWidth', 1.1, ...
                            'Color', caCol(q, :), 'DisplayName', caNm{q}, ...
                            'HitTest', 'off', 'PickableParts', 'none');
                    end
                end
                yline(axSAS, 0, 'k:', 'HandleVisibility', 'off', ...
                    'HitTest', 'off', 'PickableParts', 'none');
                hWLspec = wlMarker(axSAS, false);
                hold(axSAS, 'off');
                if cbNorm.Value
                    ylabel(axSAS, 'Normalized amplitude');
                else
                    ylabel(axSAS, 'Amplitude (\DeltaA or Intensity)');
                end
                if isDAS
                    title(axSAS, 'Decay-associated spectra (DAS) of this fit');
                else
                    title(axSAS, 'Species-associated spectra (SAS)');
                end
                xlim(axSAS, [min(app.wavelength), max(app.wavelength)]);
                legend(axSAS, 'Location', 'best', 'Box', 'off');
            end

            function plotPop()
                app.claAll(axPop);
                xlabel(axPop, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axPop, 'Population');
                grid(axPop, 'on');
                if isempty(res)
                    title(axPop, 'Populations (run fit)');
                    legend(axPop, 'off');
                    return;
                end
                cols = lines(max(res.sch.nS, 3));
                hold(axPop, 'on');
                for q = 1:res.sch.nS
                    plot(axPop, app.delay, res.conc(:, q), '-', ...
                        'LineWidth', 1.4, 'Color', cols(q, :), ...
                        'DisplayName', res.sch.species{q});
                end
                hold(axPop, 'off');
                title(axPop, 'Populations c(t)');
                applyTimeAxis(axPop);
                ylim(axPop, 'auto');
                legend(axPop, 'Location', 'best', 'Box', 'off');
            end

            function plotKin()
                app.claAll(axKin);
                [~, iw] = min(abs(app.wavelength - edWL.Value));
                hold(axKin, 'on');
                plot(axKin, app.delay, app.deltaA(iw, :), 'k.', ...
                    'MarkerSize', 8, 'DisplayName', 'Data');
                if ~isempty(res)
                    if cbContrib.Value
                        cols = lines(max(res.sch.nS, 3));
                        for q = 1:res.sch.nS
                            a = res.SAS(iw, q);
                            if ~isfinite(a) || a == 0, continue; end
                            plot(axKin, app.delay, a * res.conc(:, q), '-', ...
                                'LineWidth', 0.9, 'Color', cols(q, :), ...
                                'DisplayName', res.sch.species{q});
                        end
                        if isfield(res, 'ca') && isstruct(res.ca) && ...
                                size(res.ca.spec, 1) == numel(app.wavelength) ...
                                && all(isfinite(res.ca.spec(iw, :)))
                            caB = app.caBasis(app.delay, res.t0, ...
                                res.ca.fwhm, res.ca.order);
                            plot(axKin, app.delay, caB * res.ca.spec(iw, :).', ...
                                '--', 'LineWidth', 0.9, ...
                                'Color', [0.4 0.4 0.4], ...
                                'DisplayName', 'artifact');
                        end
                    end
                    plot(axKin, app.delay, res.fit(iw, :), 'r-', ...
                        'LineWidth', 1.5, 'DisplayName', 'Fit');
                end
                hold(axKin, 'off');
                grid(axKin, 'on');
                xlabel(axKin, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axKin, '\DeltaA or Intensity');
                title(axKin, sprintf('Kinetics at %s = %.4g %s', ...
                    app.xSym(), app.wavelength(iw), app.xUnitAx()));
                applyTimeAxis(axKin);
                ylim(axKin, 'auto');
                legend(axKin, 'Location', 'best', 'Box', 'off');
            end

            function [Rk, tK, isLog] = residualInView()
                % Residual matrix on the delays that are shown: the
                % fitted range, and t > 0 only on a log axis.
                tA = app.delay(:).';
                isLog = strcmp(ddScale.Value, 'Log') && any(tA > 0);
                keep = true(size(tA));
                if isLog, keep = tA > 0; end
                if ~isempty(res.tRange)
                    keep = keep & tA >= res.tRange(1) & tA <= res.tRange(2);
                end
                if nnz(keep) < 2
                    keep = true(size(tA));  isLog = false;
                end
                Rk = app.deltaA(:, keep) - res.fit(:, keep);
                tK = tA(keep);
            end

            function plotRes()
                app.claAll(axRes);
                colorbar(axRes, 'off');
                xlabel(axRes, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(axRes, app.xAxLabel());
                if isempty(res)
                    title(axRes, 'Residual map (run fit)');
                    axRes.ButtonDownFcn = '';     % nothing to pick on
                    return;
                end
                [Rk, tK, isLog] = residualInView();
                [Tg, Wg] = meshgrid(tK, app.wavelength(:));
                surface(axRes, Tg, Wg, zeros(size(Rk)), Rk, ...
                    'EdgeColor', 'none', 'FaceColor', 'flat', ...
                    'HitTest', 'off', 'PickableParts', 'none');
                view(axRes, 2);
                cl = max(abs(Rk(:)), [], 'omitnan');
                if isempty(cl) || ~isfinite(cl) || cl == 0, cl = 1; end
                axRes.CLim = [-cl, cl];
                colormap(axRes, app.getColormap());
                if isLog
                    axRes.XScale = 'log';
                else
                    axRes.XScale = 'linear';
                end
                xlim(axRes, [min(tK), max(tK)]);
                ylim(axRes, [min(app.wavelength), max(app.wavelength)]);
                hold(axRes, 'on');
                hWLres = wlMarker(axRes, true);
                hold(axRes, 'off');
                title(axRes, sprintf('Residual (data - fit),  RMS = %.3g', ...
                    res.rms));
                cb = colorbar(axRes);
                cb.Label.String = '\DeltaA or Intensity';
            end

            function h = wlMarker(ax, isY)
                % Line at the probe position of the kinetics panel
                % (vertical on the spectra panel, horizontal on the
                % residual map, where the spectral axis is y), and the
                % click / drag wiring of that axes. The line is kept
                % out of the legend and is transparent to the mouse.
                % The built-in drag-to-pan / scroll-to-zoom / data-tip
                % interactions are taken off this axes: they would
                % swallow the mouse sequence of the drag (zooming stays
                % available from the axes toolbar). Redone at every
                % redraw, since plotting can reset both.
                if isY
                    h = yline(ax, edWL.Value, 'k-', 'LineWidth', 1.0, ...
                        'Alpha', 0.7, 'HandleVisibility', 'off', ...
                        'HitTest', 'off', 'PickableParts', 'none');
                else
                    h = xline(ax, edWL.Value, 'k-', 'LineWidth', 1.0, ...
                        'Alpha', 0.7, 'HandleVisibility', 'off', ...
                        'HitTest', 'off', 'PickableParts', 'none');
                end
                try, disableDefaultInteractivity(ax); catch, end
                try, ax.Interactions = []; catch, end
                ax.ButtonDownFcn = @(~, ev) onPickDown(ax, isY, ev);
            end

            function onWLChange()
                % Probe position typed, clicked or dragged: move the
                % markers, redraw the kinetics, and keep the position
                % with the stored result.
                wl = edWL.Value;
                if ~isempty(hWLspec) && isvalid(hWLspec)
                    hWLspec.Value = wl;
                end
                if ~isempty(hWLres) && isvalid(hWLres)
                    hWLres.Value = wl;
                end
                plotKin();
                if isstruct(app.tgSaved) && isfield(app.tgSaved, 'curWL')
                    app.tgSaved.curWL = wl;
                end
            end

            function onPickDown(ax, isY, ev)
                % Mouse pressed on the spectra panel or the residual map.
                pt = ev.IntersectionPoint;
                if isempty(pt), return; end
                % Flag and timer BEFORE the apply (which yields via
                % drawnow): a fast release must be allowed to win.
                dragAx    = ax;
                dragIsY   = isY;
                draggingW = true;
                startPickTimer();
                applyPick(pt);
            end

            function onPickMotion()
                if ~draggingW, return; end
                if isempty(dragAx) || ~isvalid(dragAx), return; end
                cp = dragAx.CurrentPoint;
                if isempty(cp) || size(cp, 1) < 1, return; end
                applyPick(cp(1, :));
            end

            function onPickUp()
                draggingW = false;
                stopPickTimer();
            end

            function applyPick(pt)
                if numel(pt) < 2, return; end
                if dragIsY
                    wl = pt(2);
                else
                    wl = pt(1);
                end
                if ~isfinite(wl), return; end
                lim = edWL.Limits;
                edWL.Value = min(max(wl, lim(1)), lim(2));
                onWLChange();
                drawnow limitrate;
            end

            function startPickTimer()
                % Safety net of the drag: the mouse-up event is lost
                % when the button is released outside the window.
                if isempty(pickTimer) || ~isvalid(pickTimer)
                    try
                        pickTimer = timer('ExecutionMode', 'fixedRate', ...
                            'Period', 0.1, 'BusyMode', 'drop', ...
                            'TimerFcn', @(~,~) checkPick());
                    catch
                        pickTimer = [];  return;
                    end
                end
                try
                    if strcmp(pickTimer.Running, 'off'), start(pickTimer); end
                catch
                end
            end

            function stopPickTimer()
                if ~isempty(pickTimer) && isvalid(pickTimer)
                    try
                        if strcmp(pickTimer.Running, 'on'), stop(pickTimer); end
                    catch
                    end
                end
            end

            function checkPick()
                if ~draggingW, stopPickTimer();  return; end
                if ~app.mouseInsideFigure(f), onPickUp(); end
            end

            function showSyntaxHelp()
                % Syntax summary in a window of its own. It can stay
                % open next to the scheme while typing.
                if ~isempty(synFig) && isvalid(synFig)
                    figure(synFig);
                    return;
                end
                synFig = uifigure('Name', 'Scheme syntax', ...
                    'Position', app.fitFigurePosition(800, 540));
                app.applyLightTheme(synFig);
                synFig.CloseRequestFcn = @(~,~) closeSyntaxHelp();
                g3 = uigridlayout(synFig, [3, 1]);
                g3.RowHeight = {36, '1x', 22};
                g3.Padding = [12 12 12 12];  g3.RowSpacing = 8;
                uilabel(g3, 'Text', syntaxIntro, 'WordWrap', 'on', ...
                    'VerticalAlignment', 'top');
                tb3 = uitable(g3, 'Data', syntaxRows, ...
                    'ColumnName', {'Statement', 'Meaning'}, ...
                    'RowName', {}, 'ColumnEditable', [false false], ...
                    'ColumnWidth', {230, 'auto'});
                try
                    addStyle(tb3, uistyle('FontName', 'Courier New'), ...
                        'column', 1);
                catch
                end
                uilabel(g3, 'Text', syntaxNote, 'FontColor', [0.3 0.3 0.3]);
            end

            function closeSyntaxHelp()
                if ~isempty(synFig) && isvalid(synFig)
                    delete(synFig);
                end
                synFig = [];
            end

            function styleDiaAxes(ax)
                % Plain canvas for the scheme diagram: no ticks, no
                % frame, no zoom or pan.
                ax.XTick = [];  ax.YTick = [];
                ax.XColor = 'none';  ax.YColor = 'none';
                ax.Box = 'off';  ax.Color = 'w';
                try, disableDefaultInteractivity(ax); catch, end
                try, ax.Interactions = []; catch, end
            end

            function [tauS, und] = diagramValues()
                % Time constant of every step of the scheme shown, taken
                % from the last result -- only if that result was fitted
                % with exactly this scheme (steps, initial populations
                % and spectral constraints).
                tauS = [];  und = [];
                if isempty(res) || isempty(sch), return; end
                try
                    rs = res.sch;
                    same = isequal(rs.species, sch.species) && ...
                        numel(rs.steps) == numel(sch.steps) && ...
                        isequal(rs.c0(:), sch.c0(:)) && ...
                        isequal(rs.zero, sch.zero) && ...
                        isequal(rs.rel, sch.rel) && ...
                        isequal(rs.nonneg, sch.nonneg);
                    if same && ~isempty(rs.steps)
                        same = isequal([rs.steps.from], [sch.steps.from]) && ...
                            isequal([rs.steps.to], [sch.steps.to]) && ...
                            isequal([rs.steps.coef], [sch.steps.coef]) && ...
                            isequal({rs.steps.expr}, {sch.steps.expr});
                    end
                    if ~same, return; end
                    nSd = numel(rs.steps);
                    tv = zeros(1, nSd);  uv = false(1, nSd);
                    for q = 1:nSd
                        ip = rs.steps(q).par;
                        tv(q) = res.tau(ip) / rs.steps(q).mult;
                        uv(q) = ~res.tauFixed(ip) && isinf(res.tauStd(ip));
                    end
                    tauS = tv;  und = uv;
                catch
                    tauS = [];  und = [];
                end
            end

            function drawDiagram(ax, Wd, Hd, fsMax, hint)
                % Draw the applied scheme on ax, a canvas of about
                % Wd x Hd pixels. The geometry is that of
                % targetDiagramLayout; the boxes take the colours of the
                % spectra and population curves.
                if isempty(ax) || ~isvalid(ax), return; end
                app.claAll(ax);
                if isempty(sch), return; end
                [tauS, und] = diagramValues();
                Dg = app.targetDiagramLayout(sch, Wd, Hd, fsMax, tauS, und, ...
                    app.tUnitTxt());
                cols  = lines(max(sch.nS, 3));
                inkC  = [0.25 0.25 0.25];
                noHit = {'HitTest', 'off', 'PickableParts', 'none'};
                for q = 1:numel(Dg.line)
                    line('Parent', ax, 'XData', Dg.line(q).x, ...
                        'YData', Dg.line(q).y, 'Color', inkC, ...
                        'LineWidth', 1, noHit{:});
                end
                for q = 1:numel(Dg.head)
                    patch('Parent', ax, 'XData', Dg.head(q).x, ...
                        'YData', Dg.head(q).y, 'FaceColor', inkC, ...
                        'EdgeColor', 'none', noHit{:});
                end
                for q = 1:numel(Dg.rect)
                    rc = Dg.rect(q);
                    cq = cols(rc.sp, :);
                    rr = 0.4 * rc.pos(4);          % corner diameter
                    rectangle('Parent', ax, 'Position', rc.pos, ...
                        'Curvature', [min(1, rr / rc.pos(3)), min(1, rr / rc.pos(4))], ...
                        'FaceColor', 1 - 0.18 * (1 - cq), 'EdgeColor', cq, ...
                        'LineWidth', 1 + 1.2 * double(rc.init), noHit{:});
                end
                for q = 1:numel(Dg.text)
                    tx = Dg.text(q);
                    tc = [0.40 0.40 0.40];  tw = 'normal';
                    switch tx.kind
                        case 'name',  tc = [0.07 0.07 0.07];  tw = 'bold';
                        case 'value', tc = [0.60 0.25 0.45];
                        case 'rate',  tc = [0.20 0.20 0.20];
                    end
                    ht = text('Parent', ax, 'Position', [tx.x, tx.y, 0], ...
                        'String', tx.str, 'Interpreter', 'none', ...
                        'FontUnits', 'pixels', 'FontSize', tx.size, ...
                        'FontWeight', tw, 'Color', tc, ...
                        'HorizontalAlignment', tx.ha, ...
                        'VerticalAlignment', tx.va, noHit{:});
                    if tx.bg
                        ht.BackgroundColor = 'w';
                        ht.Margin = 1;
                    end
                end
                if ~isempty(hint)
                    text('Parent', ax, ...
                        'Position', [Dg.W - 0.6 * Dg.fs, Dg.noteY, 0], ...
                        'String', hint, 'Interpreter', 'none', ...
                        'FontUnits', 'pixels', 'FontSize', Dg.noteSize, ...
                        'Color', [0.45 0.45 0.45], ...
                        'HorizontalAlignment', 'right', ...
                        'VerticalAlignment', 'bottom', noHit{:});
                end
                xlim(ax, [0, Dg.W]);
                ylim(ax, [0, Dg.H]);
            end

            function drawDiagrams()
                % Redraw the scheme diagram: the tab and, if it is open,
                % the larger window. A picture that cannot be drawn must
                % not keep the scheme from being applied or fitted.
                try
                    drawDiagram(axDia, 430, 152, 11, 'click to enlarge');
                    if ~isempty(diaFig) && isvalid(diaFig)
                        ps = diaFig.Position;
                        drawDiagram(axDiaWin, max(ps(3) - 40, 200), ...
                            max(ps(4) - 40, 120), 20, '');
                    end
                catch MEd
                    lblStatus.Text = ['Diagram not drawn: ' MEd.message];
                end
            end

            function openDiagramWindow()
                % Larger view of the scheme diagram in a window of its
                % own (a click on the diagram tab opens it).
                if ~isempty(diaFig) && isvalid(diaFig)
                    figure(diaFig);
                else
                    diaFig = uifigure('Name', 'Kinetic scheme', ...
                        'Position', app.fitFigurePosition(900, 520));
                    app.applyLightTheme(diaFig);
                    diaFig.CloseRequestFcn = @(~,~) closeDiagramWindow();
                    g2 = uigridlayout(diaFig, [1, 1]);
                    g2.Padding = [10 10 10 10];
                    try, g2.BackgroundColor = 'w'; catch, end
                    axDiaWin = uiaxes(g2);
                    styleDiaAxes(axDiaWin);
                end
                drawDiagrams();
            end

            function closeDiagramWindow()
                if ~isempty(diaFig) && isvalid(diaFig)
                    delete(diaFig);
                end
                diaFig = [];
                axDiaWin = [];
            end

            function closeSVD()
                if ~isempty(svdFig) && isvalid(svdFig)
                    delete(svdFig);
                end
                svdFig = [];
            end

            function showResidualSVD()
                % Singular value decomposition of the residual matrix on
                % the fitted delays. What the model does not describe
                % shows up as structure in the first singular vectors.
                if isempty(res), return; end
                closeSVD();
                tA = app.delay(:).';
                keep = true(size(tA));
                if ~isempty(res.tRange)
                    keep = tA >= res.tRange(1) & tA <= res.tRange(2);
                end
                Rm = app.deltaA(:, keep) - res.fit(:, keep);
                rowOK = any(isfinite(Rm), 2) & all(isfinite(res.SAS), 2);
                Rm(~isfinite(Rm)) = 0;
                Rm(~rowOK, :) = 0;
                [Us, Ss, Vs] = svd(Rm, 'econ');
                sv = diag(Ss);
                nShow = min(2, numel(sv));
                if nShow < 1, return; end
                tK = tA(keep);

                svdFig = uifigure('Name', 'Residual SVD', ...
                    'Position', app.fitFigurePosition(1150, 400));
                app.applyLightTheme(svdFig);
                gs = uigridlayout(svdFig, [1, 3]);
                gs.ColumnWidth = {'0.8x', '1x', '1x'};
                gs.Padding = [10 10 10 10];  gs.ColumnSpacing = 10;
                a1 = uiaxes(gs);  a2 = uiaxes(gs);  a3 = uiaxes(gs);
                cols = lines(3);

                nSv = min(20, numel(sv));
                semilogy(a1, 1:nSv, max(sv(1:nSv), realmin), 'o-', ...
                    'LineWidth', 1.2, 'Color', [0.2 0.2 0.2], ...
                    'MarkerFaceColor', [0.2 0.2 0.2]);
                grid(a1, 'on');
                xlabel(a1, 'Index');  ylabel(a1, 'Singular value');
                title(a1, 'Singular values of the residual');

                hold(a2, 'on');  hold(a3, 'on');
                for q = 1:nShow
                    plot(a2, app.wavelength, Us(:, q), '-', 'LineWidth', 1.3, ...
                        'Color', cols(q, :), ...
                        'DisplayName', sprintf('#%d  (s = %.3g)', q, sv(q)));
                    plot(a3, tK, Vs(:, q), '-', 'LineWidth', 1.3, ...
                        'Color', cols(q, :), ...
                        'DisplayName', sprintf('#%d', q));
                end
                yline(a2, 0, 'k:', 'HandleVisibility', 'off');
                yline(a3, 0, 'k:', 'HandleVisibility', 'off');
                hold(a2, 'off');  hold(a3, 'off');
                grid(a2, 'on');  grid(a3, 'on');
                xlabel(a2, app.xAxLabel());
                ylabel(a2, 'Left singular vector');
                title(a2, 'Spectral structure of the residual');
                xlim(a2, [min(app.wavelength), max(app.wavelength)]);
                legend(a2, 'Location', 'best', 'Box', 'off');
                xlabel(a3, sprintf('Delay time (%s)', app.tUnitAx()));
                ylabel(a3, 'Right singular vector');
                title(a3, 'Temporal structure of the residual');
                if strcmp(ddScale.Value, 'Log') && any(tK > 0)
                    a3.XScale = 'log';
                    xlim(a3, [min(tK(tK > 0)), max(tK)]);
                elseif max(tK) > min(tK)
                    xlim(a3, [min(tK), max(tK)]);
                end
                legend(a3, 'Location', 'best', 'Box', 'off');
            end

            function doReset()
                res = [];
                app.tgSaved = [];
                closeSVD();
                setResultButtons('off');
                txtResult.Value = idleText;
                plotAll();
                lblStatus.Text = 'Reset.';
            end

            function saveSettings()
                [sLo, sHi] = app.fitRangeToStore(edFitTmin.Value, edFitTmax.Value);
                sv2 = struct();
                sv2.preset    = ddPreset.Value;
                sv2.schemeTxt = txtScheme.Value;
                sv2.table     = tblPar.Data;
                sv2.t0        = edT0.Value;
                sv2.fwhm      = edFW.Value;
                sv2.t0Fixed   = cbT0fix.Value;
                sv2.fwhmFixed = cbFWfix.Value;
                sv2.irf       = ddIRF.Value;
                sv2.ca        = struct('order', ddCA.Value, ...
                    'mode', ddCAW.Value, 'fwhm', edCAW.Value);
                sv2.fitTMin   = sLo;
                sv2.fitTMax   = sHi;
                sv2.scale     = ddScale.Value;
                sv2.norm      = cbNorm.Value;
                sv2.contrib   = cbContrib.Value;
                sv2.spec      = ddSpec.Value;
                app.tgSet = sv2;
            end

            function hdr = speciesHeaders()
                % Column names for exports: letters, digits, underscore.
                hdr = cell(1, res.sch.nS);
                for q = 1:res.sch.nS
                    nm = regexprep(res.sch.species{q}, '[^A-Za-z0-9_]', '_');
                    if isempty(regexp(nm, '^[A-Za-z]', 'once'))
                        nm = ['S_' nm]; %#ok<AGROW>
                    end
                    hdr{q} = nm;
                end
                % Keep the names unique after the clean-up.
                for q = 2:numel(hdr)
                    while any(strcmp(hdr(1:q-1), hdr{q}))
                        hdr{q} = [hdr{q} '_']; %#ok<AGROW>
                    end
                end
            end

            function exportSAS()
                if isempty(res), return; end
                mat = [app.wavelength(:), res.SAS];
                hdr = [{app.xHdr()}, speciesHeaders()];
                % spectra of the coherent-artifact terms, if any, follow
                if isfield(res, 'ca') && isstruct(res.ca) && ...
                        size(res.ca.spec, 1) == numel(app.wavelength)
                    caHdr = {'artifact_Gaussian', 'artifact_1st_derivative', ...
                             'artifact_2nd_derivative'};
                    kCA = min(3, size(res.ca.spec, 2));
                    mat = [mat, res.ca.spec(:, 1:kCA)];
                    hdr = [hdr, caHdr(1:kCA)];
                end
                doExportTable(mat, hdr, ...
                    'target_SAS', 'Export species-associated spectra');
            end

            function exportDAS()
                if isempty(res), return; end
                hdr = cell(1, numel(res.obsTau));
                for q = 1:numel(res.obsTau)
                    if isfinite(res.obsTau(q))
                        hdr{q} = sprintf('tau%d=%.4g_%s', q, res.obsTau(q), ...
                            app.tUnitHdr());
                    else
                        hdr{q} = sprintf('tau%d=inf', q);
                    end
                end
                doExportTable([app.wavelength(:), res.DAS], ...
                    [{app.xHdr()}, hdr], ...
                    'target_DAS', 'Export decay-associated spectra');
            end

            function exportPop()
                if isempty(res), return; end
                doExportTable([app.delay(:), res.conc], ...
                    [{sprintf('t_%s', app.tUnitHdr())}, speciesHeaders()], ...
                    'target_populations', 'Export populations');
            end

            function exportFit()
                if isempty(res), return; end
                tA = app.delay(:).';
                hdr = [{app.xHdr()}, ...
                    arrayfun(@(x) sprintf('t=%.4g_%s', x, app.tUnitHdr()), ...
                        tA, 'UniformOutput', false)];
                doExportTable([app.wavelength(:), res.fit], hdr, ...
                    'target_fit', 'Export fitted data matrix');
            end

            function exportPar()
                if isempty(res), return; end
                [file, path] = app.putFile({'*.txt', 'Text (*.txt)'}, ...
                    'Export fit parameters', 'target_parameters.txt');
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                fid = fopen(full, 'w');
                if fid < 0
                    uialert(f, sprintf('Cannot write:\n%s', full), 'Error');
                    return;
                end
                txtL = txtResult.Value;
                for q = 1:numel(txtL)
                    fprintf(fid, '%s\n', txtL{q});
                end
                fclose(fid);
                lblStatus.Text = sprintf('Parameters saved: %s', file);
            end

            function doExportTable(mat, hdr, defaultName, dlgTitle)
                [file, path] = app.putFile( ...
                    {'*.csv','CSV (*.csv)'; '*.dat','DAT (tab, *.dat)'}, ...
                    dlgTitle, [defaultName '.csv']);
                if isequal(file, 0), return; end
                full = fullfile(path, file);
                [~, ~, ext] = fileparts(full);
                if strcmpi(ext, '.csv'), dlm = ','; else, dlm = sprintf('\t'); end
                try
                    T = array2table(mat, 'VariableNames', hdr);
                    writetable(T, full, 'Delimiter', dlm);
                    uialert(f, sprintf('Saved: %s', full), ...
                        'Exported', 'Icon', 'success');
                catch ME
                    uialert(f, sprintf('Export failed:\n%s', ME.message), ...
                        'Error');
                end
            end

            function onClose()
                if ~isempty(pickTimer) && isvalid(pickTimer)
                    try, stop(pickTimer); catch, end
                    try, delete(pickTimer); catch, end
                    pickTimer = [];
                end
                try
                    saveSettings();
                catch
                end
                closeSVD();
                closeDiagramWindow();
                closeSyntaxHelp();
                app.targetFig = [];
                delete(f);
            end
        end
        
        % =================================================================
        % Cleanup
        % =================================================================
        % =================================================================
        % Automatic analysis window
        % =================================================================
        function openAutoWindow(app)
            %OPENAUTOWINDOW  Automatic analysis: a proposal for the
            %   corrections (background, mask, chirp, crop) and a summary
            %   of the kinetics without a kinetic model, as a report with
            %   the plots it rests on. The working data are changed only
            %   by "Apply corrections to the data".
            if isempty(app.deltaA)
                uialert(app.fig, 'Load data first.', 'No data');
                return;
            end
            if ~isempty(app.autoFig) && isvalid(app.autoFig)
                figure(app.autoFig);  return;
            end

            f = uifigure('Name', 'Auto Analysis', ...
                'Position', app.fitFigurePosition(1380, 860));
            app.applyLightTheme(f);
            app.autoFig = f;

            gl = uigridlayout(f, [1, 2]);
            gl.ColumnWidth = {340, '1x'};
            gl.Padding = [10 10 10 10];  gl.ColumnSpacing = 10;

            % ============== Left: setup ==============
            left = uipanel(gl, 'Title', 'Automatic analysis');
            left.Layout.Column = 1;
            lgl = uigridlayout(left, [18, 1]);
            lgl.RowHeight = {22, 24, 24, 28, 24, 24, 22, 30, 30, 24, 36, 22, ...
                             30, 32, 32, 32, 22, '1x'};
            lgl.Padding = [8 8 8 8];  lgl.RowSpacing = 6;

            uilabel(lgl, 'Text', '1. Corrections to propose', ...
                'FontWeight', 'bold');
            cbBg = uicheckbox(lgl, 'Value', true, ...
                'Text', 'Background (delays before the signal)', ...
                'Tooltip', ['Finds the delays that lie before the signal ' ...
                    'and subtracts their mean spectrum. Skipped when that ' ...
                    'mean is zero within the noise, or when fewer than 3 ' ...
                    'delays lie safely before the signal.']);
            cbMask = uicheckbox(lgl, 'Value', true, ...
                'Text', 'Mask channels with excess noise', ...
                'Tooltip', ['Channels whose noise is more than 3 times ' ...
                    'that of their neighbours (pump scatter, dead pixels, ' ...
                    'the ends of the probe spectrum) can be set to NaN and ' ...
                    'left out of every fit. Which of them are dropped is ' ...
                    'set by "Mask rule". Dropped channels at an end of the ' ...
                    'spectrum are cut away instead when "Crop" is ticked.']);
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [18 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Mask rule:');
            ddMask = uidropdown(row, ...
                'Items', {'Only where the signal is lost', ...
                          'Every noisy channel'}, ...
                'ItemsData', {'weak', 'all'}, 'Value', 'weak', ...
                'Tooltip', ['Which of the noisy channels are dropped.' newline ...
                    'Only where the signal is lost: a noisy channel is ' ...
                    'dropped when its signal is less than 10 times its ' ...
                    'noise, or when its noise is more than 10 times that of ' ...
                    'the neighbours. Channels that are noisier than their ' ...
                    'neighbours but show a clear signal stay in the data ' ...
                    '(the pump region with a strong bleach, the ends of ' ...
                    'the probe spectrum); the report lists them.' newline ...
                    'Every noisy channel: every channel more than 3 times ' ...
                    'noisier than its neighbours is dropped, together with ' ...
                    'the neighbours that are still 1.5 times noisier.']);
            cbMask.ValueChangedFcn = @(~,~) set(ddMask, 'Enable', ...
                app.pick(cbMask.Value, 'on', 'off'));
            cbChirp = uicheckbox(lgl, 'Value', true, ...
                'Text', 'Chirp (time zero along the spectrum)', ...
                'Tooltip', ['Locates the steepest rise of every wavelength, ' ...
                    'fits the chirp curve of the Chirp Corr window through ' ...
                    'them and shifts every wavelength so that its time ' ...
                    'zero is 0. If time zero is the same everywhere but ' ...
                    'not at 0, the delay axis is shifted instead.']);
            cbCrop = uicheckbox(lgl, 'Value', true, ...
                'Text', 'Crop what the corrections leave empty', ...
                'Tooltip', ['Removes the first and last delays that the ' ...
                    'chirp correction leaves without data at some ' ...
                    'wavelengths, and noisy channels at the ends of the ' ...
                    'spectrum.']);

            uilabel(lgl, 'Text', '2. Summary of the kinetics', ...
                'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', 60};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Largest number of components:');
            edMax = uieditfield(row, 'numeric', 'Value', 5, ...
                'Limits', [1 5], 'RoundFractionalValues', 'on', ...
                'Tooltip', ['Global fits are run with 1, 2, ... components ' ...
                    'until one more lowers the RMS by less than 5 % or ' ...
                    'brings a time constant the data do not determine. ' ...
                    'This is the largest number tried.']);
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Coherent artifact:');
            ddCA = uidropdown(row, ...
                'Items', {'Automatic', 'Off', 'On (3 terms)'}, ...
                'ItemsData', {'auto', 'off', 'on'}, 'Value', 'auto', ...
                'Tooltip', ['Terms for a signal confined to time zero (cross-' ...
                    'phase modulation, two-photon absorption, Raman of the ' ...
                    'solvent): a Gaussian of the IRF width and its first two ' ...
                    'derivatives, each with a free spectrum, as in the ' ...
                    'Global Analysis window. Automatic: used when the fits ' ...
                    'without them need a time constant below a third of the ' ...
                    'IRF width, or when the terms lower the RMS by 15 % or ' ...
                    'more. Off: never. On: always. The report says what ' ...
                    'was done and why.']);
            cbAll = uicheckbox(lgl, 'Value', false, ...
                'Text', 'Fit every number up to the largest', ...
                'Tooltip', ['Normally the fits stop at the first number of ' ...
                    'components that is rejected. Ticked: all numbers up to ' ...
                    'the largest are fitted, so that each of them can be ' ...
                    'picked under "Components used". The automatic choice ' ...
                    'stays the same. Takes longer.']);
            btnRun = uibutton(lgl, 'Text', 'Run', ...
                'BackgroundColor', [0.20 0.50 0.45], ...
                'FontColor', 'w', 'FontWeight', 'bold', ...
                'Tooltip', ['Analyse the data as they are shown in the ' ...
                    'main window. Nothing is changed.']);

            uilabel(lgl, 'Text', '3. Use the result', 'FontWeight', 'bold');
            row = uigridlayout(lgl, [1, 2]);
            row.ColumnWidth = {'fit', '1x'};
            row.Padding = [0 0 0 0];  row.ColumnSpacing = 6;
            uilabel(row, 'Text', 'Components used:');
            ddN = uidropdown(row, 'Items', {'(run first)'}, ...
                'Enable', 'off', ...
                'Tooltip', ['Number of components of the fit that section 3 ' ...
                    'of the report, the Components and Residual tabs and ' ...
                    '"Send model to Global Analysis" use. After a run it ' ...
                    'is the automatic choice. Every number in the table of ' ...
                    'fits can be picked instead, also one that the ' ...
                    'automatic choice rejected; the report then says so ' ...
                    'and why it was rejected. Nothing is fitted again.']);
            btnApply = uibutton(lgl, 'Enable', 'off', ...
                'Text', 'Apply corrections to the data', ...
                'Tooltip', ['Write the proposed corrections into the ' ...
                    'working data of the main window. Like Crop, this ' ...
                    'makes the corrected data the new baseline; "Reset ' ...
                    'all" in the main window returns to the loaded data.']);
            btnGA = uibutton(lgl, 'Enable', 'off', ...
                'Text', 'Send model to Global Analysis', ...
                'Tooltip', ['Open the Global Analysis window with the ' ...
                    'number of components in use ("Components used"), the ' ...
                    'fitted values as starting values and the same ' ...
                    'coherent-artifact terms. A stored global-fit result ' ...
                    'is discarded. Press "Run Fit" there.']);
            btnSave = uibutton(lgl, 'Enable', 'off', ...
                'Text', 'Save report...');
            lblStatus = uilabel(lgl, 'Text', 'Ready.', ...
                'FontColor', [0.2 0.2 0.2]);
            uitextarea(lgl, 'Editable', 'off', 'Value', { ...
                'Run analyses the data as they are shown in the main window and changes nothing.', ...
                '', ...
                'Report: what was found, with the numbers behind every decision, and a list of things to check by eye.', ...
                '', ...
                'Corrections: the rise of every wavelength with the fitted time zero, the data after the corrections, the background and the noise of every channel.', ...
                '', ...
                'Components: singular values, RMS against the number of components, kinetics at a few wavelengths (data and fit), DADS and EADS of the fit in use (dashed: coherent-artifact terms).', ...
                '', ...
                'Residual: what that fit leaves, and its spectrum (oscillations).', ...
                '', ...
                '"Components used" switches the report and the plots to a fit with another number of components.', ...
                '', ...
                'The report describes the data. It does not decide which states and reaction steps are behind the time constants.'});

            % ============== Right: report and plots ==============
            tabs = uitabgroup(gl);
            tabs.Layout.Column = 2;

            tabRep = uitab(tabs, 'Title', 'Report');
            gR = uigridlayout(tabRep, [1, 1]);
            gR.Padding = [6 6 6 6];
            txtRep = uitextarea(gR, 'Editable', 'off', ...
                'FontName', 'Monospaced', ...
                'Value', {'Press "Run".'});

            tabCor = uitab(tabs, 'Title', 'Corrections');
            gC = uigridlayout(tabCor, [2, 2]);
            gC.RowHeight = {'1.3x', '1x'};  gC.ColumnWidth = {'1x', '1x'};
            gC.Padding = [8 8 8 8];  gC.RowSpacing = 8;  gC.ColumnSpacing = 10;
            axRaw = uiaxes(gC);  axRaw.Layout.Row = 1;  axRaw.Layout.Column = 1;
            axCor = uiaxes(gC);  axCor.Layout.Row = 1;  axCor.Layout.Column = 2;
            axBg  = uiaxes(gC);  axBg.Layout.Row  = 2;  axBg.Layout.Column  = 1;
            axNz  = uiaxes(gC);  axNz.Layout.Row  = 2;  axNz.Layout.Column  = 2;

            tabCmp = uitab(tabs, 'Title', 'Components');
            gM = uigridlayout(tabCmp, [2, 2]);
            gM.RowHeight = {'1x', '1x'};  gM.ColumnWidth = {'1x', '1x'};
            gM.Padding = [8 8 8 8];  gM.RowSpacing = 8;  gM.ColumnSpacing = 10;
            % top left: singular values and RMS, side by side
            gTL = uigridlayout(gM, [1, 2]);
            gTL.Layout.Row = 1;  gTL.Layout.Column = 1;
            gTL.ColumnWidth = {'1x', '1x'};
            gTL.Padding = [0 0 0 0];  gTL.ColumnSpacing = 8;
            axSV  = uiaxes(gTL);  axSV.Layout.Row  = 1;  axSV.Layout.Column  = 1;
            axRms = uiaxes(gTL);  axRms.Layout.Row = 1;  axRms.Layout.Column = 2;
            % top right: kinetics at a few wavelengths, data and fit.
            % Two axes share the delay axis: the delays around time
            % zero on a linear scale, the later ones on a log scale.
            gKin = uigridlayout(gM, [2, 1]);
            gKin.Layout.Row = 1;  gKin.Layout.Column = 2;
            gKin.RowHeight = {24, '1x'};
            gKin.Padding = [0 0 0 0];  gKin.RowSpacing = 4;
            gKc = uigridlayout(gKin, [1, 3]);
            gKc.ColumnWidth = {'fit', 170, '1x'};
            gKc.Padding = [0 0 0 0];  gKc.ColumnSpacing = 6;
            uilabel(gKc, 'Text', sprintf('Kinetics at (%s):', app.xUnitTxt()));
            edKin = uieditfield(gKc, 'text', 'Value', '', ...
                'Tooltip', ['Positions on the spectral axis whose kinetics ' ...
                    'are shown, data and fit: up to 6 numbers separated ' ...
                    'by spaces or commas, for example  520 600 680. ' ...
                    'The nearest channel of the data is taken. Empty: ' ...
                    'chosen automatically (the channel with the largest ' ...
                    'signal in each fifth of the spectrum).']);
            % what the plot shows, or why it is empty (the two axes carry
            % no title: equal decorations keep their ordinates aligned)
            lblKin = uilabel(gKc, 'Text', '', 'FontColor', [0.35 0.35 0.35]);
            gKa = uigridlayout(gKin, [1, 2]);
            gKa.ColumnWidth = {'1x', '1.5x'};
            gKa.Padding = [0 0 0 0];  gKa.ColumnSpacing = 2;
            axKinL = uiaxes(gKa);  axKinL.Layout.Row = 1;  axKinL.Layout.Column = 1;
            axKinR = uiaxes(gKa);  axKinR.Layout.Row = 1;  axKinR.Layout.Column = 2;
            % bottom: spectra
            axDADS = uiaxes(gM);  axDADS.Layout.Row = 2;  axDADS.Layout.Column = 1;
            axEADS = uiaxes(gM);  axEADS.Layout.Row = 2;  axEADS.Layout.Column = 2;

            tabRes = uitab(tabs, 'Title', 'Residual');
            gS = uigridlayout(tabRes, [2, 2]);
            gS.RowHeight = {'1x', '1x'};  gS.ColumnWidth = {'1x', '1x'};
            gS.Padding = [8 8 8 8];  gS.RowSpacing = 8;  gS.ColumnSpacing = 10;
            axRes  = uiaxes(gS);  axRes.Layout.Row  = [1 2];  axRes.Layout.Column  = 1;
            axPow  = uiaxes(gS);  axPow.Layout.Row  = 1;      axPow.Layout.Column  = 2;
            axProf = uiaxes(gS);  axProf.Layout.Row = 2;      axProf.Layout.Column = 2;

            % ---- State ----
            P = [];  S = [];          % results of the last run
            L = {};                   % report lines
            sumErr = '';              % message when the summary failed
            applied = false;          % corrections written into the main data
            srcWl = [];  srcT = [];  srcD = [];   % data the run was made on
            srcDesc = '';             % their name in the report
            dlg = [];                 % progress dialog while running
            lineCol = lines(7);

            % ---- Wire up ----
            btnRun.ButtonPushedFcn   = @(~,~) safeCb(@doRun, 'Run');
            btnApply.ButtonPushedFcn = @(~,~) safeCb(@doApply, 'Apply corrections');
            btnGA.ButtonPushedFcn    = @(~,~) safeCb(@doSendGA, 'Send model');
            btnSave.ButtonPushedFcn  = @(~,~) safeCb(@doSave, 'Save report');
            ddN.ValueChangedFcn      = @(~,~) safeCb(@onPickN, 'Components used');
            edKin.ValueChangedFcn    = @(~,~) safeCb(@drawKinetics, 'Kinetics');
            emptyAxes();

            function safeCb(fcn, label)
                try
                    fcn();
                catch ME
                    % The window was closed while the callback was
                    % running: nothing left to report to.
                    if ~isvalid(f), return; end
                    closeDlg();
                    btnRun.Enable = 'on';
                    uialert(f, sprintf('"%s" failed:\n%s\n\n(In %s line %d)', ...
                        label, ME.message, ...
                        ME.stack(1).name, ME.stack(1).line), ...
                        'Callback error');
                    lblStatus.Text = sprintf('Error in %s.', label);
                end
            end

            function setProg(fr, msg)
                lblStatus.Text = msg;
                if ~isempty(dlg) && isvalid(dlg)
                    dlg.Value = max(0, min(1, fr));
                    dlg.Message = msg;
                end
                drawnow;
            end

            function closeDlg()
                if ~isempty(dlg) && isvalid(dlg)
                    try, close(dlg); catch, end
                end
                dlg = [];
            end

            function emptyAxes()
                axs = {axRaw, axCor, axBg, axNz, axSV, axRms, axDADS, ...
                       axEADS, axRes, axPow, axProf};
                for q = 1:numel(axs)
                    app.claAll(axs{q});
                    legend(axs{q}, 'off');
                    colorbar(axs{q}, 'off');
                    title(axs{q}, 'Press "Run"');
                end
                drawKinetics();      % empty, with its own message
            end

            % ============== Callbacks ==============
            function doRun()
                opt = struct( ...
                    'doBg',    cbBg.Value, ...
                    'doMask',  cbMask.Value, ...
                    'doChirp', cbChirp.Value, ...
                    'doCrop',  cbCrop.Value, ...
                    'maxComp', round(edMax.Value), ...
                    'ca',      ddCA.Value, ...
                    'tryAll',  cbAll.Value);
                if strcmp(ddMask.Value, 'all')
                    opt.maskSN = Inf;      % drop every noisy channel
                end
                srcWl = app.wavelength(:);
                srcT  = app.delay(:).';
                srcD  = app.deltaA;
                P = [];  S = [];  L = {};  sumErr = '';
                applied = false;
                btnApply.Enable = 'off';
                btnGA.Enable    = 'off';
                btnSave.Enable  = 'off';
                btnRun.Enable   = 'off';
                ddN.Items = {'(run first)'};
                ddN.Enable = 'off';
                try
                    dlg = uiprogressdlg(f, 'Title', 'Auto Analysis', ...
                        'Message', 'Corrections...', 'Value', 0.02);
                catch
                    dlg = [];
                end
                setProg(0.02, 'Corrections...');
                try
                    P = app.autoPreprocess(srcWl, srcT, srcD, opt);
                catch errP
                    % No result: do not leave the report and the plots
                    % of an earlier run on screen. (safeCb shows the
                    % message and switches "Run" back on.)
                    txtRep.Value = {'The analysis could not be run:', ...
                        errP.message};
                    emptyAxes();
                    rethrow(errP);
                end
                setProg(0.12, 'Singular values...');
                try
                    S = app.autoSummarize(P, opt, @setProg);
                catch errS
                    S = [];
                    sumErr = errS.message;
                end
                closeDlg();
                btnRun.Enable = 'on';

                srcDesc = app.dataSourceDesc;
                if isempty(srcDesc), srcDesc = 'current data'; end
                refreshReport();
                fillPickN();
                drawCorrections();
                drawComponents();
                drawResidual();
                btnApply.Enable = app.pick(P.changed, 'on', 'off');
                btnGA.Enable    = app.pick(~isempty(S), 'on', 'off');
                btnSave.Enable  = 'on';
                tabs.SelectedTab = tabRep;
                if P.changed
                    lblStatus.Text = 'Done. The data are unchanged so far.';
                else
                    lblStatus.Text = 'Done. No correction was needed.';
                end
            end

            function refreshReport()
                % Report text for the present state of P, S and "applied".
                L = app.autoReportLines(P, S, srcDesc);
                if ~isempty(sumErr)
                    L{end+1} = '';
                    L{end+1} = sprintf( ...
                        'The summary of the kinetics failed: %s', sumErr);
                end
                if applied
                    L{end+1} = '';
                    L{end+1} = ['The corrections were applied to the data of ' ...
                        'the main window.'];
                end
                txtRep.Value = L;
            end

            function fillPickN()
                % "Components used": one entry per fit in the table.
                if isempty(S)
                    ddN.Items = {'(no summary)'};
                    ddN.Enable = 'off';
                    return;
                end
                it = cell(1, numel(S.tabN));
                for q = 1:numel(it)
                    if q == S.nStar
                        it{q} = sprintf('%d  (automatic choice)', q);
                    elseif ~isempty(S.tabWhy{q})
                        it{q} = sprintf('%d  (rejected)', q);
                    else
                        it{q} = sprintf('%d', q);
                    end
                end
                ddN.Items = it;
                ddN.Value = it{S.nUse};
                ddN.Enable = app.pick(numel(it) > 1, 'on', 'off');
            end

            function onPickN()
                % Another number of components: no new fit, the stored
                % one becomes the model of the report and the plots.
                if isempty(S), return; end
                n = find(strcmp(ddN.Items, ddN.Value), 1);
                if isempty(n) || n == S.nUse, return; end
                lblStatus.Text = 'Updating...';
                drawnow;
                S = app.autoSelectModel(P, S, n);
                refreshReport();
                drawComponents();
                drawResidual();
                if n == S.nStar
                    lblStatus.Text = sprintf( ...
                        'In use: %d component(s), the automatic choice.', n);
                else
                    lblStatus.Text = sprintf( ...
                        'In use: %d component(s), set by hand.', n);
                end
            end

            function nm = fitName()
                % How the fit in use is called in the plot titles.
                if S.nUse == S.nStar
                    nm = 'the chosen fit';
                else
                    nm = sprintf('the fit with %d components', S.nUse);
                end
            end

            function ok = sameData()
                % Are the working data still the ones the run was made on?
                ok = isequal(app.wavelength(:), srcWl) && ...
                     isequal(app.delay(:).', srcT) && ...
                     isequaln(app.deltaA, srcD);
            end

            function doApply()
                if isempty(P) || applied, return; end
                if ~P.changed
                    uialert(f, 'The report proposes no correction.', ...
                        'Nothing to apply', 'Icon', 'info');
                    return;
                end
                if ~sameData()
                    uialert(f, ['The data in the main window changed after ' ...
                        'the run. Press "Run" again.'], 'Data changed');
                    return;
                end
                sel = uiconfirm(f, ['The corrections of the report are ' ...
                    'written into the working data and become the new ' ...
                    'baseline (as with Crop). The analysis windows close ' ...
                    'and their stored results are discarded. "Reset all" ' ...
                    'in the main window returns to the loaded data.'], ...
                    'Apply corrections', ...
                    'Options', {'Apply', 'Cancel'}, ...
                    'DefaultOption', 1, 'CancelOption', 2);
                if ~strcmp(sel, 'Apply'), return; end

                tK = P.tIn(P.keepT);
                % Correction windows that are still open show the data
                % as they were: close them.
                tools = {app.bgFig, app.chirpFig, app.cropFig, app.maskFig};
                for q = 1:numel(tools)
                    if ~isempty(tools{q}) && isvalid(tools{q})
                        try, close(tools{q}); catch, end
                    end
                end
                % This window causes the data change itself and stays
                % open: hide it from closeAnalysisWindows() meanwhile.
                app.autoFig = [];
                try
                    % Corrected data on the present axes, then the crop,
                    % which makes them the new baseline and resets the
                    % correction state, selection and display ranges.
                    app.deltaA = P.Dfull;
                    app.applyCropByRange(min(P.wl), max(P.wl), ...
                        tK(1), tK(end));
                    if P.chirp.shift ~= 0
                        app.applyZeroTimeShift(P.chirp.shift);
                    end
                catch errA
                    app.autoFig = f;
                    rethrow(errA);
                end
                app.autoFig = f;
                applied = true;
                btnApply.Enable = 'off';
                % The run now describes the data of the main window.
                srcWl = app.wavelength(:);
                srcT  = app.delay(:).';
                srcD  = app.deltaA;
                refreshReport();
                lblStatus.Text = 'Corrections applied to the main data.';
                figure(f);
            end

            function doSendGA()
                if isempty(S), return; end
                if P.changed && ~applied
                    uialert(f, ['The model belongs to the corrected data. ' ...
                        'Press "Apply corrections to the data" first.'], ...
                        'Apply first', 'Icon', 'info');
                    return;
                end
                if ~sameData()
                    uialert(f, ['The data in the main window changed after ' ...
                        'the run. Press "Run" again.'], 'Data changed');
                    return;
                end
                Mo = S.model;
                if ~isempty(app.gaFig) && isvalid(app.gaFig)
                    close(app.gaFig);
                end
                nT = numel(Mo.tau);
                app.gaNComp     = nT;
                app.gaTauInit   = Mo.tau(:).';
                app.gaTauFixed  = false(1, nT);
                app.gaBetaInit  = ones(1, nT);
                app.gaStretchOn = false(1, nT);
                app.gaHasInf    = Mo.hasInf;
                app.gaIRFTreat  = Mo.irf;
                app.gaT0        = Mo.t0;
                app.gaFWHM      = Mo.fwhm;
                if strcmp(Mo.irf, 'none') || Mo.irfSeed
                    % (irfSeed: a fit with fewer components than the
                    % automatic choice, made with the IRF of that choice)
                    app.gaT0Fixed   = true;
                    app.gaFWHMFixed = true;
                else
                    app.gaT0Fixed   = false;
                    app.gaFWHMFixed = Mo.fwFixed;
                end
                % Same coherent-artifact terms (width of the IRF)
                app.gaCA = app.caOptions(struct('order', Mo.caK, ...
                    'mode', 'irf', 'fwhm', NaN));
                % The summary leaves out delays without data at some
                % wavelength: fit the same delays.
                if all(S.cols)
                    app.gaFitTMin = [];  app.gaFitTMax = [];
                else
                    app.gaFitTMin = S.t(1);  app.gaFitTMax = S.t(end);
                end
                app.clearGaResults();
                app.openGlobalAnalysisWindow();
                lblStatus.Text = 'Model sent: press "Run Fit" there.';
            end

            function doSave()
                if isempty(L), return; end
                defName = 'TRSpec_auto_report.txt';
                base = app.dataStem();
                if ~isempty(base)
                    defName = ['TRSpec_' base '_auto_report.txt'];
                end
                [fn, fp] = app.putFile('*.txt', 'Save report', defName);
                if isequal(fn, 0), return; end
                fid = fopen(fullfile(fp, fn), 'w');
                if fid < 0
                    uialert(f, 'The file could not be opened for writing.', ...
                        'Save report');
                    return;
                end
                try
                    fprintf(fid, 'TRSpecAnalyzer, saved %s\n\n', ...
                        datestr(now, 'yyyy-mm-dd HH:MM')); %#ok<TNOW1,DATST>
                    for q = 1:numel(L)
                        fprintf(fid, '%s\n', L{q});
                    end
                catch errW
                    fclose(fid);
                    rethrow(errW);
                end
                fclose(fid);
                lblStatus.Text = sprintf('Report saved: %s', fn);
            end

            % ============== Plots ==============
            function heat(ax, wlv, tv, Z)
                % 2D map (wavelength x delay) on a symmetric colour scale.
                app.claAll(ax);
                [Wg, Tg] = meshgrid(wlv(:).', tv(:).');
                Cd = Z.';
                hs = surface(ax, Wg, Tg, zeros(size(Cd)), Cd, ...
                    'EdgeColor', 'none', 'FaceColor', 'flat');
                if any(isnan(Cd(:)))
                    hs.AlphaData = double(~isnan(Cd));
                    hs.AlphaDataMapping = 'none';
                    hs.FaceAlpha = 'flat';
                end
                view(ax, 2);
                ax.YDir = 'normal';
                cl = max(abs(Z(:)), [], 'omitnan');
                if isempty(cl) || ~isfinite(cl) || cl == 0, cl = 1; end
                ax.CLim = [-cl cl];
                colormap(ax, app.getColormap());
                colorbar(ax);
                ax.XLim = [min(wlv), max(wlv)];
                if tv(end) > tv(1), ax.YLim = [tv(1), tv(end)]; end
                ax.NextPlot = 'add';
                xlabel(ax, app.xAxLabel());
                ylabel(ax, sprintf('Delay time (%s)', app.tUnitAx()));
            end

            function drawCorrections()
                ch = P.chirp;
                wlI = P.wlIn;  tI = P.tIn;
                vr = P.rowOK & ~P.maskRows;
                Db = srcD - P.bg;

                % ---- data near time zero, before and after ----
                lo = tI(1);  hi = tI(end);  cRaw = NaN;
                if P.hasT0 && any(ch.good)
                    trg = ch.ridgeT(ch.good);
                    wdt = max([P.fwGuess, 2 * ch.dtZ, 0]);
                    if wdt > 0
                        cRaw = median(trg);
                        lo = max(tI(1),   min(trg) - 8 * wdt);
                        hi = min(tI(end), max(trg) + 12 * wdt);
                    end
                end
                cI = tI >= lo & tI <= hi;
                if nnz(cI) < 3, cI = true(size(tI)); end
                heat(axRaw, wlI, tI(cI), Db(:, cI));
                hLeg = gobjects(0);
                if any(ch.good)
                    hLeg(end+1) = line(axRaw, wlI(ch.good), ...
                        ch.ridgeT(ch.good), 'LineStyle', 'none', ...
                        'Marker', '.', 'MarkerSize', 8, 'Color', 'k', ...
                        'DisplayName', 'steepest rise');
                end
                if any(strcmp(ch.status, {'applied', 'found'}))
                    hLeg(end+1) = line(axRaw, wlI(vr), ch.t0(vr), ...
                        'Color', [0 0.6 0], 'LineWidth', 1.8, ...
                        'DisplayName', 'fitted time zero');
                end
                if isempty(hLeg)
                    legend(axRaw, 'off');
                else
                    legend(axRaw, hLeg, 'Location', 'best');
                end
                title(axRaw, 'Data before the corrections, near time zero');

                tO = P.t;
                lo2 = tO(1);  hi2 = tO(end);
                if isfinite(cRaw) && isfinite(P.t0Guess)
                    lo2 = max(tO(1),   P.t0Guess + (lo - cRaw));
                    hi2 = min(tO(end), P.t0Guess + (hi - cRaw));
                end
                cO = tO >= lo2 & tO <= hi2;
                if nnz(cO) < 3, cO = true(size(tO)); end
                heat(axCor, P.wl, tO(cO), P.D(:, cO));
                if isfinite(P.t0Guess)
                    yline(axCor, P.t0Guess, 'k:', 'HandleVisibility', 'off');
                end
                legend(axCor, 'off');
                if P.changed
                    title(axCor, 'After the proposed corrections');
                else
                    title(axCor, 'No correction proposed');
                end

                % ---- background ----
                app.claAll(axBg);
                if strcmp(P.bgStatus, 'applied')
                    plot(axBg, wlI, P.bg, '-', 'Color', lineCol(1, :), ...
                        'LineWidth', 1.2);
                    title(axBg, sprintf( ...
                        'Background: mean of the first %d delays', P.nBg));
                else
                    title(axBg, 'No background subtracted');
                end
                yline(axBg, 0, 'k:', 'HandleVisibility', 'off');
                grid(axBg, 'on');
                xlabel(axBg, app.xAxLabel());
                ylabel(axBg, '\DeltaA or Intensity');

                % ---- noise of every channel ----
                app.claAll(axNz);
                hold(axNz, 'on');
                hN = plot(axNz, wlI, P.noise, '-', 'Color', lineCol(1, :), ...
                    'LineWidth', 1.0, 'DisplayName', 'noise');
                hLeg = hN;
                thr = P.opt.maskK * P.noise ./ P.noiseRatio;
                if any(isfinite(thr))
                    hLeg(end+1) = plot(axNz, wlI, thr, '--', ...
                        'Color', [0.5 0.5 0.5], 'DisplayName', ...
                        sprintf('%.3g x neighbours', P.opt.maskK));
                end
                mk = P.maskRows & P.keepW;
                if any(mk)
                    hLeg(end+1) = plot(axNz, wlI(mk), P.noise(mk), 'o', ...
                        'Color', [0.85 0.2 0.2], 'DisplayName', 'masked');
                end
                kn = P.noisy & ~P.maskRows;
                if any(kn)
                    hLeg(end+1) = plot(axNz, wlI(kn), P.noise(kn), 's', ...
                        'Color', [0.2 0.55 0.3], 'MarkerSize', 5, ...
                        'DisplayName', 'noisy, kept (clear signal)');
                end
                if any(~P.keepW)
                    hLeg(end+1) = plot(axNz, wlI(~P.keepW), ...
                        P.noise(~P.keepW), 'x', ...
                        'Color', [0.3 0.3 0.3], 'DisplayName', 'cropped');
                end
                hold(axNz, 'off');
                axNz.YScale = 'log';
                grid(axNz, 'on');
                legend(axNz, hLeg, 'Location', 'best');
                xlabel(axNz, app.xAxLabel());
                ylabel(axNz, 'Noise (same unit as the data)');
                title(axNz, 'Noise of every channel');
            end

            function drawComponents()
                axs = {axSV, axRms, axDADS, axEADS};
                for q = 1:numel(axs)
                    app.claAll(axs{q});
                    legend(axs{q}, 'off');
                end
                if isempty(S)
                    for q = 1:numel(axs), title(axs{q}, 'No summary'); end
                    drawKinetics();
                    return;
                end
                tuA = app.tUnitAx();

                % ---- singular values ----
                % (These two axes are narrow. The title is kept short so
                % that the axes toolbar does not cover it; what the
                % lines and markers mean is in a legend, in the corner
                % that a falling curve leaves free.)
                kS = min(15, numel(S.sv));
                hold(axSV, 'on');
                plot(axSV, 1:kS, S.sv(1:kS), 'o-', 'Color', lineCol(1, :), ...
                    'MarkerFaceColor', lineCol(1, :), ...
                    'HandleVisibility', 'off');
                hasFl = isfinite(S.svFloor) && S.svFloor > 0;
                if hasFl
                    hFl = plot(axSV, [0.5, kS + 0.5], 5 * S.svFloor * [1 1], ...
                        'k--', 'DisplayName', '5 x noise floor');
                end
                hold(axSV, 'off');
                axSV.YScale = 'log';
                axSV.XLim = [0.5, kS + 0.5];
                grid(axSV, 'on');
                xlabel(axSV, 'Component');
                ylabel(axSV, 'Singular value');
                title(axSV, 'Singular values');
                if hasFl, legend(axSV, hFl, 'Location', 'northeast'); end

                % ---- RMS against the number of components ----
                hold(axRms, 'on');
                plot(axRms, S.tabN, S.tabRms, 'o-', 'Color', lineCol(1, :), ...
                    'HandleVisibility', 'off');
                if S.nUse ~= S.nStar
                    nmUse = 'in use';
                else
                    nmUse = 'chosen';
                end
                hR = plot(axRms, S.nUse, S.tabRms(S.nUse), 'o', ...
                    'Color', lineCol(1, :), ...
                    'MarkerFaceColor', lineCol(1, :), 'MarkerSize', 9, ...
                    'DisplayName', nmUse);
                if S.nUse ~= S.nStar
                    hR(end+1) = plot(axRms, S.nStar, S.tabRms(S.nStar), 's', ...
                        'Color', [0.3 0.3 0.3], 'MarkerSize', 14, ...
                        'LineWidth', 1.2, 'DisplayName', 'automatic choice');
                end
                if isfinite(S.rmsNoise) && S.rmsNoise > 0
                    hR(end+1) = plot(axRms, [0.5, max(S.tabN) + 0.5], ...
                        S.rmsNoise * [1 1], 'k:', ...
                        'DisplayName', 'noise estimate');
                end
                hold(axRms, 'off');
                axRms.YScale = 'log';
                axRms.XTick = S.tabN;
                axRms.XLim = [0.5, max(S.tabN) + 0.5];
                grid(axRms, 'on');
                xlabel(axRms, 'Number of components');
                if S.ca.use > 0
                    ylabel(axRms, 'RMS of the fit with artifact terms');
                else
                    ylabel(axRms, 'RMS of the global fit');
                end
                title(axRms, 'RMS of the fits');
                legend(axRms, hR, 'Location', 'northeast');

                % ---- spectra of the chosen fit ----
                Mo = S.model;
                nT = numel(Mo.tau);
                nc = size(Mo.DADS, 2);
                hD = gobjects(1, nc);  hE = gobjects(1, nc);
                hold(axDADS, 'on');  hold(axEADS, 'on');
                for k = 1:nc
                    col = lineCol(mod(k - 1, size(lineCol, 1)) + 1, :);
                    if k <= nT
                        nm = sprintf('\\tau_%d = %.3g %s', k, Mo.tau(k), tuA);
                    else
                        nm = 'non-decaying';
                    end
                    hD(k) = plot(axDADS, P.wl, Mo.DADS(:, k), '-', ...
                        'Color', col, 'LineWidth', 1.2, 'DisplayName', nm);
                    hE(k) = plot(axEADS, P.wl, Mo.EADS(:, k), '-', ...
                        'Color', col, 'LineWidth', 1.2, ...
                        'DisplayName', sprintf('EADS %d', k));
                end
                % Coherent-artifact terms: dashed, in gray
                caNm  = {'artifact: Gaussian', 'artifact: 1st derivative', ...
                         'artifact: 2nd derivative'};
                caCol = [0.15 0.15 0.15; 0.45 0.45 0.45; 0.70 0.70 0.70];
                for k = 1:Mo.caK
                    hD(end+1) = plot(axDADS, P.wl, Mo.caSpec(:, k), '--', ...
                        'Color', caCol(k, :), 'LineWidth', 1.0, ...
                        'DisplayName', caNm{k}); %#ok<AGROW>
                end
                yline(axDADS, 0, 'k:', 'HandleVisibility', 'off');
                yline(axEADS, 0, 'k:', 'HandleVisibility', 'off');
                hold(axDADS, 'off');  hold(axEADS, 'off');
                grid(axDADS, 'on');  grid(axEADS, 'on');
                legend(axDADS, hD, 'Location', 'best');
                legend(axEADS, hE, 'Location', 'best');
                xlabel(axDADS, app.xAxLabel());
                xlabel(axEADS, app.xAxLabel());
                ylabel(axDADS, 'Amplitude');
                ylabel(axEADS, 'Amplitude');
                if Mo.caK > 0
                    title(axDADS, ['DADS of ' fitName() ...
                        '  (dashed: artifact terms)']);
                else
                    title(axDADS, ['DADS of ' fitName()]);
                end
                title(axEADS, ['EADS of ' fitName() ...
                    '  (1 \rightarrow 2 \rightarrow ...)']);
                drawKinetics();
            end

            function idx = kinPick(wv, Dk)
                % Channels whose kinetics are shown: those nearest to
                % the numbers typed in "Kinetics at", or, when the field
                % is empty, the channel with the largest signal in each
                % fifth of the spectrum (strong traces from all over the
                % spectrum).
                vals = sscanf(regexprep(strtrim(edKin.Value), '[,;]', ' '), ...
                    '%f').';
                isAuto = isempty(vals);
                if isAuto
                    amp = max(abs(Dk), [], 2);
                    ed  = linspace(min(wv), max(wv), 6);
                    vals = zeros(1, 0);
                    for b = 1:5
                        inB = find(wv >= ed(b) & wv <= ed(b+1));
                        if isempty(inB), continue; end
                        [~, im] = max(amp(inB));
                        vals(end+1) = wv(inB(im)); %#ok<AGROW>
                    end
                end
                span = max(wv) - min(wv);
                idx = zeros(1, 0);
                for q = 1:numel(vals)
                    if ~isfinite(vals(q)), continue; end
                    [dW, i1] = min(abs(wv - vals(q)));
                    if dW > 0.02 * span, continue; end   % outside the data
                                                         % or in a masked range
                    if any(idx == i1), continue; end
                    if isAuto && ~isempty(idx) && ...
                            any(abs(wv(idx) - wv(i1)) < 0.03 * span)
                        continue;        % automatic picks: not side by side
                    end
                    idx(end+1) = i1; %#ok<AGROW>
                    if numel(idx) >= 6, break; end
                end
                if isAuto, idx = sort(idx); end
            end

            function drawKinetics()
                % Data (dots) and fit (lines) of the model in use at a
                % few wavelengths. The evenly sampled delays around time
                % zero go on a linear axis (left), the later, more widely
                % spaced ones on a logarithmic axis (right).
                % The two axes share the ordinate. Their plot boxes line
                % up only if both carry the same decorations above and
                % below: no title (the text goes into lblKin), tick
                % labels of one height (plain numbers on the log axis
                % too), one line of axis label.
                app.claAll(axKinL);  legend(axKinL, 'off');
                app.claAll(axKinR);  legend(axKinR, 'off');
                title(axKinL, '');  title(axKinR, '');
                axKinL.XScale = 'linear';  axKinR.XScale = 'linear';
                axKinR.YTickLabelMode = 'auto';
                axKinR.XTickMode = 'auto';  axKinR.XTickLabelMode = 'auto';
                if isempty(S)
                    gKa.ColumnWidth = {'1x', 0};
                    axKinR.Visible = 'off';  axKinL.Visible = 'on';
                    if isempty(P)
                        lblKin.Text = 'Press "Run"';
                    else
                        lblKin.Text = 'No summary';
                    end
                    return;
                end
                Mo = S.model;
                wv = S.wl(:);  tk = S.t(:).';  n = numel(tk);
                Dk = P.D(S.rows, S.cols);
                idx = kinPick(wv, Dk);

                % ---- where the linear part ends ----
                % The last delay before the spacing exceeds three times
                % the spacing just after time zero.
                tz = 0;
                if strcmp(Mo.irf, 'numerical') && isfinite(Mo.t0), tz = Mo.t0; end
                j0 = find(tk >= tz, 1);
                if isempty(j0), j0 = 1; end
                jS = [];
                if n >= 8 && j0 <= n - 2
                    d  = diff(tk);
                    d0 = median(d(j0:min(j0 + 4, n - 1)));
                    jS = find(d(j0:end) > 3 * d0, 1) + j0 - 1;
                end
                kMode = 'lin';
                if ~isempty(jS) && tk(jS) > 0 && n - jS >= 5 && ...
                        tk(end) >= 5 * tk(jS)
                    kMode = 'split';
                    if jS < 5 && tk(1) > 0, kMode = 'log'; end
                end
                switch kMode
                    case 'split'
                        iL = 1:jS;  iR = jS:n;
                        gKa.ColumnWidth = {'1x', '1.5x'};
                        axKinL.Visible = 'on';  axKinR.Visible = 'on';
                    case 'log'
                        iL = zeros(1, 0);  iR = 1:n;
                        gKa.ColumnWidth = {0, '1x'};
                        axKinL.Visible = 'off';  axKinR.Visible = 'on';
                    otherwise
                        iL = 1:n;  iR = zeros(1, 0);
                        gKa.ColumnWidth = {'1x', 0};
                        axKinL.Visible = 'on';  axKinR.Visible = 'off';
                end
                if isempty(iR), axT = axKinL; else, axT = axKinR; end
                if isempty(idx)
                    lblKin.Text = 'No channel at these positions';
                    return;
                end

                % ---- plot ----
                hK = gobjects(1, numel(idx));
                yLo = Inf;  yHi = -Inf;
                hold(axKinL, 'on');  hold(axKinR, 'on');
                for q = 1:numel(idx)
                    i1  = idx(q);
                    col = lineCol(mod(q - 1, size(lineCol, 1)) + 1, :);
                    nm  = sprintf('%.4g %s', wv(i1), app.xUnitAx());
                    if ~isempty(iL)
                        plot(axKinL, tk(iL), Dk(i1, iL), 'o', 'Color', col, ...
                            'MarkerSize', 3.5, 'HandleVisibility', 'off');
                        hK(q) = plot(axKinL, tk(iL), Mo.fit(i1, iL), '-', ...
                            'Color', col, 'LineWidth', 1.2, 'DisplayName', nm);
                    end
                    if ~isempty(iR)
                        plot(axKinR, tk(iR), Dk(i1, iR), 'o', 'Color', col, ...
                            'MarkerSize', 3.5, 'HandleVisibility', 'off');
                        hK(q) = plot(axKinR, tk(iR), Mo.fit(i1, iR), '-', ...
                            'Color', col, 'LineWidth', 1.2, 'DisplayName', nm);
                    end
                    yLo = min([yLo, Dk(i1, :), Mo.fit(i1, :)]);
                    yHi = max([yHi, Dk(i1, :), Mo.fit(i1, :)]);
                end
                if ~isempty(iL), yline(axKinL, 0, 'k:', 'HandleVisibility', 'off'); end
                if ~isempty(iR), yline(axKinR, 0, 'k:', 'HandleVisibility', 'off'); end
                hold(axKinL, 'off');  hold(axKinR, 'off');
                pad = 0.05 * (yHi - yLo);
                if ~(pad > 0), pad = 1; end
                yl = [yLo - pad, yHi + pad];
                tLab = sprintf('Delay time (%s)', app.tUnitAx());
                if ~isempty(iL)
                    axKinL.YLim = yl;
                    if tk(iL(end)) > tk(iL(1))
                        axKinL.XLim = [tk(iL(1)), tk(iL(end))];
                    end
                    grid(axKinL, 'on');
                    xlabel(axKinL, tLab);
                    ylabel(axKinL, '\DeltaA or Intensity');
                end
                if ~isempty(iR)
                    axKinR.XScale = 'log';
                    axKinR.YLim = yl;
                    axKinR.XLim = [tk(iR(1)), tk(iR(end))];
                    % decades written as plain numbers (10, 100, ...):
                    % the usual 10^n labels are taller than the labels
                    % of the linear axes and would shift the plot box
                    pw = ceil(log10(tk(iR(1)))):floor(log10(tk(iR(end))));
                    if ~isempty(pw)
                        axKinR.XTick = 10 .^ pw;
                        axKinR.XTickLabel = arrayfun(@(v) sprintf('%g', v), ...
                            10 .^ pw, 'UniformOutput', false);
                    end
                    grid(axKinR, 'on');
                    if isempty(iL)
                        xlabel(axKinR, tLab);
                        ylabel(axKinR, '\DeltaA or Intensity');
                    else
                        % same ordinate as the left axes: no labels twice
                        axKinR.YTickLabel = {};
                        xlabel(axKinR, 'log scale');
                        xlabel(axKinL, [tLab ', linear']);
                    end
                end
                % Outside the plot box: there it covers no trace, and it
                % follows the axes when the layout changes ('best' is
                % worked out once, for the size the axes has at that
                % moment, and ended up beyond the edge of the tab).
                legend(axT, hK, 'Location', 'northeastoutside');
                lblKin.Text = ['dots: data,  lines: ' fitName()];
            end

            function drawResidual()
                app.claAll(axPow);   legend(axPow, 'off');
                app.claAll(axProf);  legend(axProf, 'off');
                if isempty(S)
                    app.claAll(axRes);
                    colorbar(axRes, 'off');
                    title(axRes, 'No summary');
                    title(axPow, 'No summary');
                    title(axProf, 'No summary');
                    return;
                end
                O = S.osc;
                Rv = S.model.resid;
                tv = S.t;
                if all(isfinite(O.win))
                    cR = tv >= O.win(1) & tv <= O.win(2);
                    ttl = ['Residual of ' fitName() ', evenly sampled part'];
                else
                    cR = true(size(tv));
                    ttl = ['Residual of ' fitName()];
                end
                if nnz(cR) < 3, cR = true(size(tv)); end
                heat(axRes, S.wl, tv(cR), Rv(:, cR));
                legend(axRes, 'off');
                title(axRes, ttl);

                if isempty(O.freq)
                    switch O.status
                        case 'short'
                            msg = 'Too few evenly spaced delays';
                        case 'skipped'
                            msg = 'Oscillations not analysed';
                        otherwise
                            msg = 'No oscillation analysis';
                    end
                    title(axPow, msg);
                    title(axProf, msg);
                    return;
                end
                hold(axPow, 'on');
                hP = plot(axPow, O.freq, O.power, '-', ...
                    'Color', lineCol(1, :), 'DisplayName', 'residual');
                hP(end+1) = plot(axPow, O.freq, O.floor, '--', ...
                    'Color', [0.5 0.5 0.5], 'DisplayName', 'local floor');
                if ~isempty(O.pkIdx)
                    hP(end+1) = plot(axPow, O.pkF, O.power(O.pkIdx), 'v', ...
                        'Color', [0.85 0.2 0.2], ...
                        'MarkerFaceColor', [0.85 0.2 0.2], ...
                        'DisplayName', 'detected');
                end
                hold(axPow, 'off');
                axPow.YScale = 'log';
                axPow.XLim = [0, O.fMax];
                grid(axPow, 'on');
                legend(axPow, hP, 'Location', 'best');
                xlabel(axPow, 'Wavenumber (cm^{-1})');
                ylabel(axPow, '|FFT|^2 summed over the wavelengths (noise-weighted)');
                title(axPow, 'Spectrum of the residual');

                if isempty(O.pkIdx)
                    title(axProf, 'No oscillation detected');
                    return;
                end
                hold(axProf, 'on');
                hQ = gobjects(1, numel(O.pkIdx));
                for q = 1:numel(O.pkIdx)
                    pr = O.Pmap(:, O.pkIdx(q));
                    pm = max(pr);
                    if ~(pm > 0), pm = 1; end
                    hQ(q) = plot(axProf, S.wl, pr / pm, '-', ...
                        'Color', lineCol(mod(q - 1, size(lineCol, 1)) + 1, :), ...
                        'LineWidth', 1.2, 'DisplayName', ...
                        sprintf('%.0f cm^{-1}', O.pkF(q)));
                end
                hold(axProf, 'off');
                grid(axProf, 'on');
                legend(axProf, hQ, 'Location', 'best');
                xlabel(axProf, app.xAxLabel());
                ylabel(axProf, '|FFT|^2 (each scaled to 1)');
                title(axProf, 'Where each oscillation is strong');
            end
        end
        
        function onClose(app)
            % Close every child window (the list in closeChildWindows is
            % the single source of truth) before the main figure goes.
            app.closeChildWindows();
            delete(app.fig);
        end
    end
end
