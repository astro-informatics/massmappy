# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
import pyssht as ssht
cimport cy_mass_mapping as mm
from libc.math cimport log, exp, sqrt

#----------------------------------------------------------------------------------------------------#

cdef extern from "so3.h":

	int so3_sampling_f_size(so3_parameters_t *params);

	ctypedef struct so3_parameters_t:
		int reality
		int L0
		int L
		int N
		so3_sampling_t sampling_scheme;
		ssht_dl_method_t dl_method;

	ctypedef enum so3_sampling_t:
		SO3_SAMPLING_MW, SO3_SAMPLING_MW_SS, SO3_SAMPLING_SIZE

#----------------------------------------------------------------------------------------------------#

cdef extern from "s2let.h":

	void fill_so3_parameters(so3_parameters_t *parameters1, s2let_parameters_t *parameters2)
	int s2let_n_phi(const s2let_parameters_t *parameters);
	int s2let_n_theta(const s2let_parameters_t *parameters);
	void s2let_tiling_direction_allocate(double complex **s_elm, const s2let_parameters_t *parameters);
	void s2let_tiling_direction(double complex *s_elm, const s2let_parameters_t *parameters);

	void s2let_tiling_wavelet_allocate(double complex **psi, double **phi, const s2let_parameters_t *parameters);
	void s2letc_tiling_wavelet(double complex *psi, double *phi, const s2let_parameters_t *parameters);
	void s2let_tiling_axisym(double *kappa, double *kappa0, const s2let_parameters_t *parameters);
	void s2let_tiling_axisym_allocate(double **kappa, double **kappa0, const s2let_parameters_t *parameters);

	int s2let_n_scal(const s2let_parameters_t *parameters);
	int s2let_n_wav(const s2let_parameters_t *parameters);

	int s2let_n_wav_j(int j, const s2let_parameters_t *parameters);

	void s2let_mw_alm2map(double complex * f, const double complex * flm, int L, int spin);

	void s2let_mw_map2alm(double complex * flm, const double complex * f, int L, int spin);

	int s2let_bandlimit(int j, const s2let_parameters_t *parameters);

	void s2let_synthesis_wav2lm(
		double complex *flm,
		const double complex *f_wav,
		const double complex *f_scal,
		const s2let_parameters_t *parameters
	);

	void s2let_analysis_lm2wav(
		double complex *f_wav,
		double complex *f_scal,
		const double complex *flm,
		const s2let_parameters_t *parameters
	);

	void s2let_analysis_lm2wav_manual(
		double complex *f_wav,
		double complex *f_scal,
		const double complex *flm,
		const double *scal_l,
		const double complex *wav_lm,
		const int scal_bandlimit,
		const int *wav_bandlimits,
		int J,
		int L,
		int spin,
		int N
	);

	void s2let_synthesis_wav2lm_manual(
		double complex *flm,
		const double complex *f_wav,
		const double complex *f_scal,
		const double *scal_l,
		const double complex *wav_lm,
		const int scal_bandlimit,
		const int *wav_bandlimits,
		int J,
		int L,
		int spin,
		int N
	);

	void s2let_transform_axisym_lm_allocate_wav(
		double **wav_lm, double **scal_lm, const s2let_parameters_t *parameters);

	void s2let_transform_axisym_lm_wav(
		double *wav_lm, double *scal_lm, const s2let_parameters_t *parameters);

	void s2let_transform_axisym_lm_wav_analysis(
		double complex *f_wav_lm,
		double complex *f_scal_lm,
		const double complex *flm,
		const double *wav_lm,
		const double *scal_lm,
		const s2let_parameters_t *parameters
	);

	void s2let_transform_axisym_lm_wav_synthesis(
		double complex *flm,
		const double complex *f_wav_lm,
		const double complex *f_scal_lm,
		const double *wav_lm,
		const double *scal_lm,
		const s2let_parameters_t *parameters
	);

	int s2let_j_max(s2let_parameters_t *parameters);

	ctypedef enum ssht_dl_method_t:
		SSHT_DL_RISBO, SSHT_DL_TRAPANI
	ctypedef enum s2let_sampling_t:
		S2LET_SAMPLING_MW
	ctypedef enum s2let_wav_norm_t:
		S2LET_WAV_NORM_DEFAULT, S2LET_WAV_NORM_SPIN_LOWERED

	ctypedef struct s2let_parameters_t:
		int J_min
		double B
		int L
		int N
		int upsample
		int spin
		ssht_dl_method_t dl_method
		s2let_wav_norm_t normalisation;
		s2let_sampling_t sampling_scheme;
		int original_spin
		int reality
		int verbosity

#---------------------------------------------
# functions
def j_max(int B, int L, int J_min):
	cdef s2let_parameters_t parameters = {};
	parameters.B = B;
	parameters.L = L;
	parameters.J_min = J_min;
	return s2let_j_max(&parameters);

def scale_maskes(np.ndarray[double, ndim=2, mode="c"] mask_orig not None, int B, int L, int N, int J_min):

	cdef np.ndarray[complex, ndim=1] eb_sep_mask

	cdef int J_max=j_max(B, L, J_min)

	eb_sep_mask = np.zeros((J_max-J_min+1)*N*L*L, dtype=complex)

	return eb_sep_mask

# (J-J_min+1)*N*L*L
# [j, n, l , m]
