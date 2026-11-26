"""trspec: calculation code of TRSpecAnalyzer in Python (no GUI)."""
from ._util import Result
from ._optim import fminsearch
from .kinetics import (exp_irf_conv, stretched_irf_conv, gauss_legendre,
                       build_ga_basis, build_ga_basis_ext, irf_mode_norm,
                       ca_options, ca_basis, irf_delay_step, irf_unresolved,
                       irf_width_lost, irf_too_wide)
from .chirp import (chirp_model, fit_chirp_params, delay_runs, spline_curv,
                    interp_runs, chirp_shift)
from .uncertainty import fit_param_uncertainty
from .global_analysis import fit_global_analysis
from .sequential import (compute_eads_from_dads, seq_clock_conc,
                         compute_sas_sequential, suggest_n_components)

from .single import fit_single_trace
from .dataset import (Dataset, read_data_file, write_data_file, detect_axis_swap,
                      compute_ridge, movmean, default_chirp_tmax, estimate_scan_shifts)
from .target import (parse_target_scheme, target_scheme_text, target_rate_matrix,
                     target_conc, fit_target_analysis, target_step_table)
from .lda import compute_lda, compute_lcurve
from .mcr import compute_mcr
from .coherence import (compute_coherence, compute_time_frequency, coh_spectrum, solvent_vc,
                        solvent_shifted, solvent_band, fit_solvent_coherence)
from .auto import (auto_options, auto_preprocess, auto_summarize,
                   auto_select_model, auto_report_lines)
from .bandpos import (lm_fit, lognorm_shape, band_axis, band_parabola,
                      band_position, band_model_curve, band_shift_start,
                      fit_band_shift)

__version__ = '1.1.0'      # follows the MATLAB program it is a port of
