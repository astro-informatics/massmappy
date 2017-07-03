import numpy as np
import matplotlib.pyplot as plt
import pyssht as ssht
import cy_mass_mapping as mm
import cy_healpy_mass_mapping as hp_mm
import cy_eb_sep as eb_mm
from matplotlib import cm, colors, colorbar, gridspec
import healpy as hp

B = 2
L = 64
N = 2
J_min = 2

Method = "MW"

mask = np.ones(ssht.sample_shape(L,Method=Method))

eb_mm.scale_maskes(mask, B, L, N, J_min)

