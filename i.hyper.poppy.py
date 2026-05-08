#!/usr/bin/env python
##############################################################################
# MODULE:    i.hyper.poppy
# AUTHOR(S): Hyperspectral Opium Poppy Detection
# PURPOSE:   Find pixels in a hyperspectral 3D raster that match opium poppy
#            (Papaver somniferum L.) using built-in multi-stage canopy reference
#            spectra derived from peer-reviewed field spectrometry studies, plus
#            any user-supplied poppy spectra.  Outputs a single max-abundance
#            map: the best match score across all reference spectra for every
#            pixel.
#
#            Built-in reference spectra are derived from:
#              • Jia et al. (2011) IEEE Trans. Geosci. Remote Sens. 49(9):3414–
#                3422 — canopy-level ASD FieldSpec Pro measurements (350–2500 nm)
#                at three growth stages in NW China.  Key optimal bands:
#                  before_flowering : 438, 528, 736, 754, 1207 nm (JM ≈ 1.20)
#                  flowering        : 405, 424, 524, 760, 1974 nm (JM ≈ 1.22)
#                  harvesting       : 468, 726, 746, 982, 1220, 1689 nm (JM ≈ 1.21)
#              • Calderón et al. (2014) Precision Agric. — leaf-level integrating-
#                sphere measurements (350–800 nm) of asymptomatic poppy leaves.
#
#            Growth stages (Papaver somniferum, canopy level, NW China):
#              before_flowering — vegetative rosette; low visible reflectance,
#                typical broadleaf red-edge, NIR plateau
#              flowering        — white petals dominate visible (R550 ≈ 0.18,
#                R670 ≈ 0.15 vs 0.04 at BF); red-edge shifts slightly
#              harvesting       — maturing capsules + senescing leaves; elevated
#                mid-infrared, red-edge persists, water bands differ
#              leaf_healthy     — leaf-level asymptomatic spectrum (350–800 nm)
#
# COPYRIGHT: (C) 2026 by the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later
##############################################################################

# %module
# % description: Opium poppy (Papaver somniferum) detection — max-abundance map from built-in multi-stage spectral reference match in a hyperspectral 3D raster
# % keyword: imagery
# % keyword: hyperspectral
# % keyword: spectroscopy
# % keyword: poppy
# % keyword: opium
# % keyword: Papaver somniferum
# % keyword: target detection
# % keyword: spectral matching
# % keyword: similarity
# % keyword: SAM
# % keyword: compendium
# % keyword: classification
# %end

# %option G_OPT_R3_INPUT
# % key: input
# % required: yes
# % description: Input hyperspectral 3D raster (from i.hyper.import or i.hyper.atcorr)
# % guisection: Input
# %end

# %option G_OPT_R_OUTPUT
# % key: output
# % required: yes
# % description: Output max-abundance raster map (0=no match, 1=perfect match to best poppy spectrum)
# % guisection: Output
# %end

# %option
# % key: growth_stage
# % type: string
# % required: no
# % options: before_flowering,flowering,harvesting,leaf_healthy,all
# % answer: all
# % description: Poppy growth stage(s) to match against — before_flowering: vegetative canopy (Jia et al. 2011, optimal 438/528/736/754/1207 nm) | flowering: white-petal canopy (Jia et al. 2011, optimal 405/424/524/760/1974 nm) | harvesting: capsule+senescence canopy (Jia et al. 2011, optimal 468/726/746/982/1220/1689 nm) | leaf_healthy: asymptomatic leaf level (Calderón et al. 2014, 350–800 nm) | all: all four built-in spectra (default)
# % guisection: Reference
# %end

# %option G_OPT_F_INPUT
# % key: poppy_db
# % required: no
# % description: Optional supplementary poppy spectrum CSV (wavelength,reflectance per line) or JSON; appended to built-in spectra
# % guisection: Reference
# %end

# %option
# % key: method
# % type: string
# % required: no
# % multiple: yes
# % options: sam,sid,sid_sam,ed,sad,sca,cr_sam,cr_ed,gd1,gd2,xcorr,dtw,ssim,jsd,bhatt,mtf,cem,ensemble,consensus
# % answer: sam
# % description: Spectral similarity method(s) — sam: Spectral Angle Mapper | sid: Spectral Information Divergence | sid_sam: SID×tan(SAM) hybrid | ed: Euclidean Distance (L2) | sad: Spectral Absolute Difference (L1) | sca: Spectral Correlation Angle (Pearson r) | cr_sam: Continuum-Removed SAM | cr_ed: Continuum-Removed Euclidean Distance | gd1: 1st-Derivative Shape Matching | gd2: 2nd-Derivative Shape Matching | xcorr: Normalized Cross-Correlation (shift-tolerant) | dtw: Dynamic Time Warping | ssim: Spectral Structural Similarity Index | jsd: Jensen-Shannon Divergence | bhatt: Bhattacharyya Coefficient | mtf: Matched/Tuned Filter | cem: Constrained Energy Minimization | ensemble: Rank-based ensemble fusion | consensus: Full multi-method consensus with calibrated probability fusion (runs all base methods, see fusion_mode=)
# % guisection: Methods
# %end

# %option
# % key: fusion_mode
# % type: string
# % required: no
# % options: rank_product,fisher,stouffer,group_product,harmonic,min
# % answer: rank_product
# % description: Fusion strategy for consensus= method — rank_product: weighted geometric mean of rank fractions (default, robust) | fisher: Fisher chi-squared combined probability test (statistically rigorous) | stouffer: Stouffer weighted Z-score combination (handles diversity weights naturally) | group_product: geometric-mean AND within method-type groups then arithmetic-mean OR across groups | harmonic: harmonic mean (strictest, requires all methods to agree) | min: minimum across methods (absolute strictest)
# % guisection: Methods
# %end

# %option
# % key: agreement_threshold
# % type: double
# % required: no
# % answer: 0.80
# % description: Calibrated-probability threshold for per-pixel agreement count (consensus= only); fraction of methods that must exceed this to count as "agreeing"
# % guisection: Methods
# %end

# %option
# % key: output_prefix
# % type: string
# % required: no
# % description: Output map prefix for individual per-stage and per-method similarity maps (written alongside the main output)
# % guisection: Output
# %end

# %option
# % key: resample
# % type: string
# % required: no
# % options: linear,cubic,pchip
# % answer: linear
# % description: Interpolation method for resampling reference spectrum to sensor wavelengths
# % guisection: Processing
# %end

# %option
# % key: normalize
# % type: string
# % required: no
# % options: none,area,max,minmax,vector
# % answer: none
# % description: Spectrum normalization before matching: none=raw reflectance; area=divide by band sum; max=divide by maximum; minmax=0-1 range; vector=unit L2 norm
# % guisection: Processing
# %end

# %option
# % key: shift_window
# % type: integer
# % required: no
# % answer: 3
# % description: Maximum band-shift window for shift-tolerant methods (xcorr, dtw); 0 disables shift tolerance
# % guisection: Processing
# %end

# %option
# % key: min_wavelength
# % type: double
# % required: no
# % description: Minimum wavelength to consider (nm); defaults to sensor range
# % guisection: Processing
# %end

# %option
# % key: max_wavelength
# % type: double
# % required: no
# % description: Maximum wavelength to consider (nm); defaults to sensor range
# % guisection: Processing
# %end

# %option
# % key: coordinates
# % type: string
# % required: no
# % description: East,North coordinates for single-pixel point analysis (requires -p flag)
# % guisection: Point mode
# %end

# %flag
# % key: n
# % description: Only use bands marked valid (valid=1) in band metadata
# % guisection: Processing
# %end

# %flag
# % key: i
# % description: Info mode: print band coverage and built-in spectrum summary, then exit
# % guisection: Processing
# %end

# %flag
# % key: v
# % description: Verbose: print processing details and per-band information
# % guisection: Processing
# %end

# %flag
# % key: c
# % description: Apply convex-hull continuum removal to both reference and pixel spectra before matching
# % guisection: Processing
# %end

# %flag
# % key: p
# % description: Point mode: analyse a single pixel at coordinates= and print full per-method score table
# % guisection: Point mode
# %end

# %flag
# % key: z
# % description: Normalize both reference and pixel spectra to unit sum (probability simplex) before matching; recommended with sid, jsd, bhatt
# % guisection: Processing
# %end


from __future__ import annotations

import sys
import os
import re
import csv
import json
import ctypes
import ctypes.util
import atexit
from typing import Optional

import numpy as np
import grass.script as gs

# ---------------------------------------------------------------------------
# Method metadata
# ---------------------------------------------------------------------------

ALL_METHODS = [
    'sam', 'sid', 'sid_sam', 'ed', 'sad', 'sca',
    'cr_sam', 'cr_ed', 'gd1', 'gd2',
    'xcorr', 'dtw', 'ssim', 'jsd', 'bhatt',
    'mtf', 'cem', 'ensemble', 'consensus',
]

BASE_METHODS = [m for m in ALL_METHODS if m not in ('ensemble', 'consensus')]

METHOD_LABELS = {
    'sam':       'Spectral Angle Mapper',
    'sid':       'Spectral Information Divergence',
    'sid_sam':   'SID × tan(SAM) hybrid',
    'ed':        'Euclidean Distance (L2)',
    'sad':       'Spectral Absolute Difference (L1)',
    'sca':       'Spectral Correlation Angle (Pearson r)',
    'cr_sam':    'Continuum-Removed SAM',
    'cr_ed':     'Continuum-Removed Euclidean Distance',
    'gd1':       '1st-Derivative Shape Matching',
    'gd2':       '2nd-Derivative Shape Matching',
    'xcorr':     'Normalized Cross-Correlation',
    'dtw':       'Dynamic Time Warping',
    'ssim':      'Spectral Structural Similarity Index',
    'jsd':       'Jensen-Shannon Divergence',
    'bhatt':     'Bhattacharyya Coefficient',
    'mtf':       'Matched Tuned Filter',
    'cem':       'Constrained Energy Minimization',
    'ensemble':  'Rank-based Ensemble Fusion',
    'consensus': 'Multi-method Consensus (calibrated probability fusion)',
}

METHOD_GROUPS: dict[str, list[str]] = {
    'geometric':   ['sam', 'sca', 'cr_sam', 'gd1', 'gd2'],
    'distance':    ['ed', 'sad', 'cr_ed'],
    'information': ['sid', 'sid_sam', 'jsd'],
    'statistical': ['bhatt', 'ssim'],
    'signal':      ['xcorr', 'dtw'],
    'subpixel':    ['mtf', 'cem'],
}

FUSION_MODES = (
    'rank_product', 'fisher', 'stouffer', 'group_product', 'harmonic', 'min',
)

GLOBAL_STATS_METHODS = {'mtf', 'cem'}

# ---------------------------------------------------------------------------
# Built-in Papaver somniferum reference spectra
# ---------------------------------------------------------------------------
# All canopy spectra: Jia et al. (2011) IEEE Trans. Geosci. Remote Sens.
#   49(9):3414–3422.  ASD FieldSpec Pro (350–2500 nm), canopy level (nadir,
#   1.5 m height), NW China arid farm.  Atmospheric water bands excluded from
#   analysis: 1360–1420 nm, 1801–1965 nm, 2451–2500 nm.
#   Reflectance values digitised from Figure 3 and calibrated against Table III/V.
#
# Leaf spectrum: Calderón et al. (2014) Precision Agric.
#   Li-Cor 1800-12 integrating sphere, 350–800 nm.
# ---------------------------------------------------------------------------

# P_BF — Before Flowering (June 18; plant height 65 cm)
# Dominant features: standard broadleaf green vegetation; optimal bands at
# 438, 528 nm (visible), 736, 754 nm (red-edge), 1207 nm (water band NIR)
_P_BF_WLS = np.array([
    350, 360, 370, 380, 390, 400, 410, 420, 430, 438, 445, 450,
    460, 470, 480, 490, 500, 510, 520, 528, 535, 550, 560, 570,
    580, 590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690,
    700, 710, 720, 726, 730, 736, 740, 746, 750, 754, 760, 770,
    780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 970, 1000,
    1050, 1100, 1150, 1200, 1207, 1250, 1300, 1350, 1360, 1420,
    1450, 1500, 1550, 1600, 1650, 1689, 1700, 1750, 1800, 1870,
    1960, 1965, 2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350,
    2400, 2450,
], dtype=np.float64)

_P_BF_REF = np.array([
    0.041, 0.041, 0.041, 0.041, 0.041, 0.042, 0.041, 0.040, 0.040, 0.040,
    0.041, 0.043, 0.046, 0.049, 0.052, 0.057, 0.063, 0.072, 0.080, 0.083,
    0.087, 0.090, 0.088, 0.084, 0.078, 0.072, 0.066, 0.059, 0.053, 0.048,
    0.045, 0.042, 0.041, 0.040, 0.041, 0.048, 0.068, 0.128, 0.232, 0.282,
    0.325, 0.362, 0.385, 0.415, 0.427, 0.435, 0.442, 0.448, 0.452, 0.458,
    0.460, 0.461, 0.461, 0.459, 0.456, 0.452, 0.445, 0.436, 0.430, 0.450,
    0.455, 0.457, 0.452, 0.442, 0.438, 0.400, 0.340, 0.222, 0.158, 0.308,
    0.312, 0.308, 0.300, 0.285, 0.270, 0.252, 0.248, 0.220, 0.079, 0.057,
    0.155, 0.158, 0.205, 0.218, 0.222, 0.216, 0.200, 0.183, 0.166, 0.149,
    0.122, 0.096,
], dtype=np.float64)

# P_FL — Flowering (July 5; plant height 125 cm; white petals)
# Dominant features: large white petals elevate visible reflectance strongly;
# optimal bands at 405, 424, 524 nm (flower pigments), 760 nm (red-edge),
# 1974 nm (petal water content; JM=1.01 though small it is unique to FL stage)
_P_FL_WLS = np.array([
    350, 360, 370, 380, 390, 400, 405, 410, 420, 424, 430, 440,
    450, 460, 470, 480, 490, 500, 510, 520, 524, 530, 540, 550,
    560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 660, 670,
    680, 690, 700, 710, 720, 730, 740, 750, 760, 770, 780, 800,
    830, 860, 900, 940, 970, 1000, 1050, 1100, 1150, 1200, 1250,
    1300, 1350, 1360, 1420, 1450, 1500, 1550, 1600, 1650, 1689,
    1700, 1750, 1800, 1870, 1960, 1965, 1974, 2000, 2050, 2100,
    2150, 2200, 2250, 2300, 2350, 2400, 2450,
], dtype=np.float64)

_P_FL_REF = np.array([
    0.098, 0.100, 0.102, 0.105, 0.110, 0.115, 0.118, 0.122, 0.130, 0.135,
    0.140, 0.145, 0.148, 0.152, 0.156, 0.160, 0.164, 0.168, 0.174, 0.179,
    0.181, 0.181, 0.180, 0.178, 0.175, 0.170, 0.165, 0.160, 0.156, 0.152,
    0.149, 0.148, 0.148, 0.148, 0.148, 0.148, 0.150, 0.160, 0.182, 0.248,
    0.340, 0.408, 0.448, 0.464, 0.470, 0.474, 0.476, 0.480, 0.481, 0.480,
    0.474, 0.462, 0.448, 0.468, 0.472, 0.474, 0.469, 0.457, 0.418, 0.354,
    0.230, 0.165, 0.320, 0.324, 0.318, 0.310, 0.294, 0.278, 0.260, 0.256,
    0.228, 0.083, 0.061, 0.160, 0.162, 0.164, 0.215, 0.228, 0.232, 0.226,
    0.210, 0.192, 0.174, 0.156, 0.128, 0.100,
], dtype=np.float64)

# P_HV — Harvesting (July 14; plant height 130 cm; capsules + senescing leaves)
# Dominant features: globose green capsules + yellowing/drying leaves;
# optimal bands at 468 nm (leaf pigment change), 726, 746 nm (red-edge),
# 982, 1220 nm (water-content NIR), 1689 nm (cellulose/lignin SWIR)
_P_HV_WLS = np.array([
    350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460,
    468, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580,
    590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700,
    710, 720, 726, 730, 740, 746, 750, 760, 780, 800, 840, 880,
    920, 960, 982, 1000, 1050, 1100, 1150, 1200, 1220, 1250, 1300,
    1350, 1360, 1420, 1450, 1500, 1550, 1600, 1650, 1689, 1700,
    1750, 1800, 1870, 1960, 1965, 2000, 2050, 2100, 2150, 2200,
    2250, 2300, 2350, 2400, 2450,
], dtype=np.float64)

_P_HV_REF = np.array([
    0.052, 0.052, 0.052, 0.052, 0.053, 0.054, 0.055, 0.057, 0.059, 0.062,
    0.066, 0.071, 0.075, 0.081, 0.088, 0.092, 0.098, 0.104, 0.110, 0.116,
    0.122, 0.120, 0.116, 0.110, 0.104, 0.098, 0.091, 0.085, 0.079, 0.075,
    0.072, 0.070, 0.068, 0.068, 0.074, 0.092, 0.150, 0.250, 0.292, 0.322,
    0.358, 0.374, 0.384, 0.395, 0.412, 0.422, 0.426, 0.426, 0.422, 0.410,
    0.400, 0.416, 0.422, 0.424, 0.420, 0.410, 0.402, 0.366, 0.310, 0.202,
    0.150, 0.300, 0.304, 0.298, 0.288, 0.274, 0.260, 0.246, 0.242, 0.216,
    0.075, 0.053, 0.152, 0.155, 0.200, 0.214, 0.218, 0.212, 0.197, 0.180,
    0.164, 0.148, 0.122, 0.096,
], dtype=np.float64)

# P_LEAF — Asymptomatic (healthy) leaf, integrating sphere
# Calderón et al. (2014) Precision Agric., Fig. 1a — Li-Cor 1800-12, 350–800 nm
# Leaf-level (not canopy), useful for matching leaf-scale hyperspectral data
_P_LEAF_WLS = np.array([
    350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460,
    470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580,
    590, 600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700,
    705, 710, 715, 720, 725, 730, 740, 750, 760, 770, 780, 790, 800,
], dtype=np.float64)

_P_LEAF_REF = np.array([
    0.042, 0.042, 0.043, 0.043, 0.044, 0.045, 0.046, 0.046, 0.047, 0.048,
    0.050, 0.053, 0.057, 0.065, 0.072, 0.080, 0.087, 0.094, 0.098, 0.100,
    0.100, 0.098, 0.092, 0.086, 0.080, 0.074, 0.068, 0.060, 0.054, 0.050,
    0.046, 0.044, 0.042, 0.042, 0.046, 0.058, 0.078, 0.136, 0.252, 0.382,
    0.446, 0.470, 0.484, 0.490, 0.492, 0.493, 0.493, 0.492, 0.492,
], dtype=np.float64)

# Registry: label → (wavelengths, reflectances, citation)
_BUILTIN_SPECTRA = {
    'before_flowering': (
        _P_BF_WLS, _P_BF_REF,
        'P_BF canopy — Jia et al. (2011) IEEE TGRS 49(9), NW China, ASD FieldSpec Pro'
    ),
    'flowering': (
        _P_FL_WLS, _P_FL_REF,
        'P_FL canopy — Jia et al. (2011) IEEE TGRS 49(9), NW China, ASD FieldSpec Pro'
    ),
    'harvesting': (
        _P_HV_WLS, _P_HV_REF,
        'P_HV canopy — Jia et al. (2011) IEEE TGRS 49(9), NW China, ASD FieldSpec Pro'
    ),
    'leaf_healthy': (
        _P_LEAF_WLS, _P_LEAF_REF,
        'P_LEAF healthy — Calderón et al. (2014) Precision Agric., Li-Cor 1800-12'
    ),
}

# ---------------------------------------------------------------------------
# Temporary raster cleanup
# ---------------------------------------------------------------------------

_TMP_RASTERS: list[str] = []


def _cleanup():
    if _TMP_RASTERS:
        gs.run_command('g.remove', type='raster',
                       name=','.join(_TMP_RASTERS), flags='f', quiet=True)


atexit.register(_cleanup)


def _tmp(label: str) -> str:
    name = f"tmp_ihpoppy_{os.getpid()}_{label}"
    _TMP_RASTERS.append(name)
    return name

# ---------------------------------------------------------------------------
# libgrass_g3d: fast band extraction
# ---------------------------------------------------------------------------

_G3D_LIB = None


def _load_g3d_lib():
    global _G3D_LIB
    if _G3D_LIB is not None:
        return _G3D_LIB

    lib_dir = os.path.join(os.environ['GISBASE'], 'lib')

    import re as _re
    _vh = os.path.join(os.environ['GISBASE'], 'include', 'grass', 'version.h')
    with open(_vh) as _f:
        _m = _re.search(r'#define\s+GRASS_HEADERS_VERSION\s+"([^"]+)"', _f.read())
    headers_version = (_m.group(1).encode() if _m else b"")

    gis_lib = ctypes.CDLL(os.path.join(lib_dir, 'libgrass_gis.so'),
                          mode=ctypes.RTLD_GLOBAL)
    gis_lib.G__no_gisinit.restype = None
    gis_lib.G__no_gisinit.argtypes = [ctypes.c_char_p]
    gis_lib.G__no_gisinit(headers_version)

    lib = ctypes.CDLL(os.path.join(lib_dir, 'libgrass_g3d.so'))
    lib.Rast3d_extract_z_slice.restype = ctypes.c_int
    lib.Rast3d_extract_z_slice.argtypes = [
        ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
    ]
    _G3D_LIB = lib
    return lib


def extract_band(raster3d: str, band_num: int) -> str:
    lib = _load_g3d_lib()
    z = band_num - 1
    base = raster3d.replace('@', '_').replace('#', '_').replace('.', '_')
    tmp_name = _tmp(f"band_{base}_{band_num}")
    name3d, mapset3d = (raster3d.split('@') + [''])[:2]
    ret = lib.Rast3d_extract_z_slice(
        name3d.encode(),
        mapset3d.encode() if mapset3d else None,
        ctypes.c_int(z),
        tmp_name.encode(),
    )
    if ret != 0:
        gs.fatal(f"Rast3d_extract_z_slice failed for band {band_num} of {raster3d}")
    return tmp_name

# ---------------------------------------------------------------------------
# Band metadata helpers
# ---------------------------------------------------------------------------


def _parse_wl_from_r3info(raster3d: str) -> dict[int, tuple]:
    try:
        info_text = gs.read_command('r3.info', flags='h', map=raster3d)
    except Exception:
        try:
            info_text = gs.read_command('r3.info', map=raster3d)
        except Exception:
            return {}
    pat = re.compile(r'Band\s+(\d+):\s+([\d.]+)\s+nm[,\s]+FWHM:\s+([\d.]+)\s+nm')
    bands: dict[int, tuple] = {}
    for line in info_text.split('\n'):
        m = pat.search(line)
        if m:
            bands[int(m.group(1))] = (float(m.group(2)), float(m.group(3)))
    return bands


def _convert_wl_nm(wl: float, unit: str) -> float:
    u = unit.lower().strip()
    if u in ('nm', 'nanometer', 'nanometers'):
        return wl
    if u in ('um', 'µm', 'micrometer', 'micron', 'microns'):
        return wl * 1000.0
    if u in ('m', 'meter', 'meters'):
        return wl * 1e9
    gs.warning(f"Unknown wavelength unit '{unit}'; assuming nm")
    return wl


def get_band_info(raster3d: str, only_valid: bool = False,
                  min_wl: Optional[float] = None,
                  max_wl: Optional[float] = None) -> list[dict]:
    info = gs.raster3d_info(raster3d)
    depths = int(info['depths'])

    base = raster3d.split('@')[0]
    mapset = raster3d.split('@')[1] if '@' in raster3d else None
    slices_exist = bool(gs.find_file(f"{base}#1", element='cell',
                                     mapset=mapset).get('name'))

    bands: list[dict] = []

    if slices_exist:
        for i in range(1, depths + 1):
            band_name = f"{raster3d}#{i}"
            wl = fwhm = None
            valid = True
            unit = 'nm'
            try:
                result = gs.read_command('r.info', map=band_name, flags='h')
                for line in result.split('\n'):
                    line = line.strip()
                    if line.startswith('wavelength='):
                        wl = float(line.split('=')[1])
                    elif line.startswith('FWHM='):
                        fwhm = float(line.split('=')[1])
                    elif line.startswith('valid='):
                        valid = int(line.split('=')[1]) == 1
                    elif line.startswith('unit='):
                        unit = line.split('=')[1].strip()
            except Exception:
                pass
            if wl is None:
                continue
            wl_nm = _convert_wl_nm(wl, unit)
            if min_wl is not None and wl_nm < min_wl:
                continue
            if max_wl is not None and wl_nm > max_wl:
                continue
            if only_valid and not valid:
                continue
            bands.append({'band': i, 'wavelength': wl_nm,
                          'fwhm': fwhm or 10.0, 'valid': valid,
                          'map_name': band_name})
    else:
        wl_dict = _parse_wl_from_r3info(raster3d)
        for i in range(1, depths + 1):
            if i not in wl_dict:
                continue
            wl_nm, fwhm = wl_dict[i]
            if min_wl is not None and wl_nm < min_wl:
                continue
            if max_wl is not None and wl_nm > max_wl:
                continue
            if only_valid:
                continue
            bands.append({'band': i, 'wavelength': wl_nm,
                          'fwhm': fwhm, 'valid': True,
                          'map_name': None})

    if not bands:
        gs.fatal(f"No wavelength metadata found in '{raster3d}'. "
                 "Import data with i.hyper.import or ensure Band N/FWHM lines "
                 "are in r3.info history.")

    bands.sort(key=lambda b: b['wavelength'])
    return bands

# ---------------------------------------------------------------------------
# Reference spectrum I/O
# ---------------------------------------------------------------------------


def get_builtin_spectra(growth_stage: str) -> list[dict]:
    """Return built-in Papaver somniferum spectra for the requested stage(s).

    :param str growth_stage: one of ``before_flowering``, ``flowering``,
        ``harvesting``, ``leaf_healthy``, or ``all``
    :return: list of spectrum dicts with keys ``wavelengths``, ``reflectances``,
        ``stage``, ``citation``
    :rtype: list[dict]
    """
    if growth_stage == 'all':
        keys = list(_BUILTIN_SPECTRA.keys())
    elif growth_stage in _BUILTIN_SPECTRA:
        keys = [growth_stage]
    else:
        gs.fatal(f"Unknown growth_stage '{growth_stage}'. "
                 f"Valid: {', '.join(list(_BUILTIN_SPECTRA) + ['all'])}")

    spectra = []
    for key in keys:
        wls, refs, citation = _BUILTIN_SPECTRA[key]
        spectra.append({
            'wavelengths':  wls,
            'reflectances': refs,
            'stage':        key,
            'citation':     citation,
        })
    return spectra


def load_poppy_db(path: str) -> list[dict]:
    """Load supplementary poppy spectra from a CSV or JSON file.

    **CSV format**: two columns ``wavelength,reflectance``.  Non-numeric rows
    are skipped.  Multiple spectra can be separated by blank lines or all
    concatenated as a single spectrum.

    **JSON formats**: same as the inline reference parser — array of pairs
    ``[[wl,r],...]`` or parallel arrays ``{"wavelengths":[...],"reflectances":[...]}``.

    :param str path: path to the supplementary spectrum file
    :return: list of spectrum dicts (each has ``wavelengths``, ``reflectances``,
        ``stage='external'``, ``citation=path``)
    :rtype: list[dict]
    """
    if not os.path.isfile(path):
        gs.fatal(f"poppy_db file not found: {path}")

    # Try JSON first
    try:
        with open(path) as f:
            obj = json.load(f)
        if isinstance(obj, list) and obj and isinstance(obj[0], (list, tuple)):
            pairs = [(float(r[0]), float(r[1])) for r in obj]
            pairs.sort()
            wls = np.array([p[0] for p in pairs], dtype=np.float64)
            refs = np.array([p[1] for p in pairs], dtype=np.float64)
            return [{'wavelengths': wls, 'reflectances': refs,
                     'stage': 'external', 'citation': path}]
        if isinstance(obj, dict):
            wls_j = obj.get('wavelengths') or obj.get('wavelength') or obj.get('wl')
            refs_j = obj.get('reflectances') or obj.get('reflectance') or obj.get('r')
            if wls_j and refs_j:
                pairs = sorted(zip(map(float, wls_j), map(float, refs_j)))
                wls = np.array([p[0] for p in pairs], dtype=np.float64)
                refs = np.array([p[1] for p in pairs], dtype=np.float64)
                return [{'wavelengths': wls, 'reflectances': refs,
                         'stage': 'external', 'citation': path}]
    except (json.JSONDecodeError, KeyError, TypeError, IndexError, ValueError):
        pass

    # CSV: collect all numeric rows as one spectrum
    pairs = []
    with open(path, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or len(row) < 2:
                continue
            try:
                pairs.append((float(row[0]), float(row[1])))
            except ValueError:
                continue
    if len(pairs) < 2:
        gs.fatal(f"Could not parse ≥2 wavelength,reflectance pairs from {path}")
    pairs.sort()
    wls = np.array([p[0] for p in pairs], dtype=np.float64)
    refs = np.array([p[1] for p in pairs], dtype=np.float64)
    return [{'wavelengths': wls, 'reflectances': refs,
             'stage': 'external', 'citation': path}]

# ---------------------------------------------------------------------------
# Wavelength LUT
# ---------------------------------------------------------------------------


class WavelengthLUT:
    """Precomputed resampling lookup table between two wavelength grids."""

    def __init__(self, src_wls, dst_wls, fill: str | float = 'edge') -> None:
        src_wls = np.asarray(src_wls, dtype=np.float64)
        dst_wls = np.asarray(dst_wls, dtype=np.float64)

        if src_wls.ndim != 1 or dst_wls.ndim != 1:
            raise ValueError("WavelengthLUT: both wavelength arrays must be 1-D")
        if src_wls.size < 2 or dst_wls.size < 2:
            raise ValueError("WavelengthLUT: need at least 2 points in each grid")
        if np.any(np.diff(src_wls) <= 0) or np.any(np.diff(dst_wls) <= 0):
            raise ValueError("WavelengthLUT: wavelength arrays must be strictly increasing")

        self.src_wls = src_wls
        self.dst_wls = dst_wls
        self.fill = fill
        n_src = len(src_wls)

        ridx = np.searchsorted(src_wls, dst_wls, side='right')
        ridx = np.clip(ridx, 1, n_src - 1)
        lidx = ridx - 1

        dx = src_wls[ridx] - src_wls[lidx]
        alpha = np.where(dx > 0, (dst_wls - src_wls[lidx]) / dx, 0.0)
        alpha = np.clip(alpha, 0.0, 1.0)

        self.left_idx  = lidx
        self.right_idx = ridx
        self.alpha     = alpha

        self.valid_dst = ((dst_wls >= src_wls[0]) & (dst_wls <= src_wls[-1]))
        self.valid_src = ((src_wls >= dst_wls[0]) & (src_wls <= dst_wls[-1]))

        wl_lo = max(src_wls[0], dst_wls[0])
        wl_hi = min(src_wls[-1], dst_wls[-1])
        self.overlap_lo: float = float(wl_lo)
        self.overlap_hi: float = float(wl_hi)
        self.has_overlap: bool = bool(wl_lo < wl_hi)

        self.overlap_dst_idx = np.where(self.valid_dst)[0]
        self.overlap_src_idx = np.where(self.valid_src)[0]
        self.overlap_dst_wls = dst_wls[self.overlap_dst_idx]
        self.overlap_src_wls = src_wls[self.overlap_src_idx]

    def apply(self, src_vals: np.ndarray) -> np.ndarray:
        src_vals = np.asarray(src_vals, dtype=np.float64)
        result = (src_vals[self.left_idx] * (1.0 - self.alpha) +
                  src_vals[self.right_idx] * self.alpha)
        self._apply_fill(result)
        return result

    def apply_cube(self, cube: np.ndarray) -> np.ndarray:
        left_vals  = cube[self.left_idx]
        right_vals = cube[self.right_idx]
        alpha = self.alpha
        for _ in range(cube.ndim - 1):
            alpha = alpha[..., np.newaxis]
        result = left_vals * (1.0 - alpha) + right_vals * alpha
        self._apply_fill_nd(result)
        return result

    def _apply_fill(self, result: np.ndarray) -> None:
        oor = ~self.valid_dst
        if not np.any(oor):
            return
        if self.fill == 'edge':
            pass
        elif self.fill == 'nan':
            result[oor] = np.nan
        else:
            result[oor] = float(self.fill)

    def _apply_fill_nd(self, result: np.ndarray) -> None:
        oor = ~self.valid_dst
        if not np.any(oor):
            return
        if self.fill == 'edge':
            return
        fill_val = np.nan if self.fill == 'nan' else float(self.fill)
        result[oor] = fill_val

    def restrict_to_overlap(self, src_vals=None, dst_bands=None):
        ovl_src_wls  = self.overlap_src_wls
        ovl_dst_wls  = self.overlap_dst_wls
        ovl_src_vals = (src_vals[self.overlap_src_idx]
                        if src_vals is not None else None)
        ovl_dst_bands = ([dst_bands[i] for i in self.overlap_dst_idx]
                         if dst_bands is not None else None)
        sub_lut = WavelengthLUT(ovl_src_wls, ovl_dst_wls, fill=self.fill)
        return ovl_src_wls, ovl_src_vals, ovl_dst_wls, ovl_dst_bands, sub_lut

    def coverage_report(self) -> str:
        n_dst = len(self.dst_wls)
        n_src = len(self.src_wls)
        n_cov = int(self.valid_dst.sum())
        n_vis = int(self.valid_src.sum())
        pct_cov = 100.0 * n_cov / n_dst if n_dst else 0.0
        pct_vis = 100.0 * n_vis / n_src if n_src else 0.0
        return (
            f"src [{self.src_wls[0]:.1f}–{self.src_wls[-1]:.1f} nm, n={n_src}] → "
            f"dst [{self.dst_wls[0]:.1f}–{self.dst_wls[-1]:.1f} nm, n={n_dst}]  |  "
            f"overlap [{self.overlap_lo:.1f}–{self.overlap_hi:.1f} nm]  |  "
            f"dst covered {n_cov}/{n_dst} ({pct_cov:.1f}%)  "
            f"src visible {n_vis}/{n_src} ({pct_vis:.1f}%)"
        )

    def __repr__(self) -> str:
        return f"WavelengthLUT({self.coverage_report()})"


def resample_reference(ref_wls, ref_vals, sensor_wls,
                       method: str = 'linear',
                       lut: Optional['WavelengthLUT'] = None) -> np.ndarray:
    if method == 'linear':
        if lut is None:
            lut = WavelengthLUT(ref_wls, sensor_wls)
        return lut.apply(ref_vals)

    from scipy.interpolate import CubicSpline, PchipInterpolator

    clipped = np.clip(sensor_wls, ref_wls[0], ref_wls[-1])
    if method == 'cubic':
        interp = CubicSpline(ref_wls, ref_vals, bc_type='not-a-knot',
                             extrapolate=False)
        result = interp(clipped)
    elif method == 'pchip':
        interp = PchipInterpolator(ref_wls, ref_vals, extrapolate=False)
        result = interp(clipped)
    else:
        gs.fatal(f"Unknown resample method: {method}")

    bad = ~np.isfinite(result)
    if np.any(bad):
        result[bad] = np.interp(clipped[bad], ref_wls, ref_vals)

    if lut is None:
        lut = WavelengthLUT(ref_wls, sensor_wls)
    oor = ~lut.valid_dst
    if np.any(oor):
        if lut.fill == 'nan':
            result[oor] = np.nan
        elif lut.fill != 'edge':
            result[oor] = float(lut.fill)

    return result.astype(np.float64)

# ---------------------------------------------------------------------------
# Cube I/O
# ---------------------------------------------------------------------------


def load_cube(bands: list[dict], raster3d: str,
              verbose: bool = False) -> np.ndarray:
    import grass.script.array as garray

    n = len(bands)
    cube = None

    for idx, b in enumerate(bands):
        if verbose:
            gs.verbose(f"  Loading band {b['band']} ({b['wavelength']:.1f} nm) ...")
        gs.percent(idx, n, 5)

        map_name = b.get('map_name')
        if map_name is None:
            map_name = extract_band(raster3d, b['band'])
            b['map_name'] = map_name

        arr = garray.array()
        arr.read(map_name)
        band_data = np.asarray(arr, dtype=np.float32)

        if cube is None:
            rows, cols = band_data.shape
            cube = np.empty((n, rows, cols), dtype=np.float32)

        cube[idx] = band_data

    gs.percent(n, n, 5)
    return cube


def read_pixel_spectrum(bands: list[dict], raster3d: str,
                        east: float, north: float) -> np.ndarray:
    coords = f"{east},{north}"
    spec = []
    for b in bands:
        map_name = b.get('map_name')
        if map_name is None:
            map_name = extract_band(raster3d, b['band'])
            b['map_name'] = map_name
        result = gs.read_command('r.what', map=map_name, coordinates=coords)
        lines = result.strip().split('\n')
        val = np.nan
        for line in lines:
            parts = line.split('|')
            if len(parts) >= 4:
                try:
                    val = float(parts[3])
                except ValueError:
                    val = np.nan
                break
        spec.append(val)
    return np.array(spec, dtype=np.float64)


def write_raster(data: np.ndarray, name: str, overwrite: bool = True) -> None:
    import grass.script.array as garray
    arr = garray.array()
    arr[:] = data.astype(np.float64)
    arr.write(name, overwrite=overwrite)

# ---------------------------------------------------------------------------
# Spectrum preprocessing
# ---------------------------------------------------------------------------


def normalize_spectrum(spec: np.ndarray, method: str) -> np.ndarray:
    if method == 'none':
        return spec
    s = spec.copy()
    if method == 'area':
        total = np.sum(s)
        if total > 0:
            s /= total
    elif method == 'max':
        m = np.max(s)
        if m > 0:
            s /= m
    elif method == 'minmax':
        lo, hi = np.min(s), np.max(s)
        if hi > lo:
            s = (s - lo) / (hi - lo)
        else:
            s = np.zeros_like(s)
    elif method == 'vector':
        n = np.linalg.norm(s)
        if n > 0:
            s /= n
    else:
        gs.fatal(f"Unknown normalize method: {method}")
    return s


def normalize_cube(cube: np.ndarray, method: str) -> np.ndarray:
    if method == 'none':
        return cube
    c = cube.copy()
    if method == 'area':
        total = np.sum(c, axis=0, keepdims=True)
        total = np.where(total > 0, total, 1.0)
        c /= total
    elif method == 'max':
        mx = np.max(c, axis=0, keepdims=True)
        mx = np.where(mx > 0, mx, 1.0)
        c /= mx
    elif method == 'minmax':
        lo = np.min(c, axis=0, keepdims=True)
        hi = np.max(c, axis=0, keepdims=True)
        rng = np.where(hi > lo, hi - lo, 1.0)
        c = (c - lo) / rng
    elif method == 'vector':
        norms = np.linalg.norm(c, axis=0, keepdims=True)
        norms = np.where(norms > 0, norms, 1.0)
        c /= norms
    return c


def to_prob_simplex(spec: np.ndarray) -> np.ndarray:
    s = np.maximum(spec, 1e-12)
    return s / s.sum()


def to_prob_simplex_cube(cube: np.ndarray) -> np.ndarray:
    c = np.maximum(cube, 1e-12)
    totals = c.sum(axis=0, keepdims=True)
    return c / totals

# ---------------------------------------------------------------------------
# Continuum removal (upper convex hull)
# ---------------------------------------------------------------------------


def _upper_hull(wavelengths: np.ndarray, reflectances: np.ndarray) -> np.ndarray:
    pts = list(zip(wavelengths.tolist(), reflectances.tolist()))
    hull: list[tuple] = []
    for p in pts:
        while len(hull) >= 2:
            o, a, b = hull[-2], hull[-1], p
            cross = (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
            if cross >= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    hull_wl = np.array([h[0] for h in hull])
    hull_rf = np.array([h[1] for h in hull])
    return np.interp(wavelengths, hull_wl, hull_rf)


def continuum_remove(spec: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
    continuum = _upper_hull(wavelengths, spec)
    with np.errstate(invalid='ignore', divide='ignore'):
        cr = np.where(continuum > 0, spec / continuum, 1.0)
    return np.clip(cr, 0.0, 1.0)


def continuum_remove_cube(cube: np.ndarray, wavelengths: np.ndarray) -> np.ndarray:
    n_bands, rows, cols = cube.shape
    n_pix = rows * cols
    flat = cube.reshape(n_bands, n_pix).T
    result = np.empty_like(flat)
    for i in range(n_pix):
        spec = flat[i]
        continuum = _upper_hull(wavelengths, spec)
        with np.errstate(invalid='ignore', divide='ignore'):
            cr = np.where(continuum > 0, spec / continuum, 1.0)
        result[i] = np.clip(cr, 0.0, 1.0)
    return result.T.reshape(n_bands, rows, cols)

# ---------------------------------------------------------------------------
# Spectral similarity / distance methods (all return similarity in [0, 1])
# ---------------------------------------------------------------------------

EPS = 1e-12


def _safe_cube(cube: np.ndarray) -> np.ndarray:
    return np.where(np.isfinite(cube), cube, 0.0)


def match_sam(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    dot = np.einsum('i,ijk->jk', ref, cube)
    norm_ref = np.linalg.norm(ref) + EPS
    norm_cube = np.sqrt(np.einsum('ijk,ijk->jk', cube, cube)) + EPS
    cos_a = np.clip(dot / (norm_ref * norm_cube), -1.0, 1.0)
    theta = np.arccos(cos_a)
    return np.clip(1.0 - 2.0 * theta / np.pi, 0.0, 1.0)


def match_sid(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    p = ref / (ref.sum() + EPS) + EPS
    q = cube / (cube.sum(axis=0, keepdims=True) + EPS) + EPS
    p3 = p[:, np.newaxis, np.newaxis]
    kld_pq = np.sum(p3 * np.log(p3 / q), axis=0)
    kld_qp = np.sum(q * np.log(q / p3), axis=0)
    sid = kld_pq + kld_qp
    return np.clip(np.exp(-np.abs(sid)), 0.0, 1.0)


def match_sid_sam(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    dot = np.einsum('i,ijk->jk', ref, cube)
    norm_ref = np.linalg.norm(ref) + EPS
    norm_cube = np.sqrt(np.einsum('ijk,ijk->jk', cube, cube)) + EPS
    theta = np.arccos(np.clip(dot / (norm_ref * norm_cube), -1.0, 1.0))
    tan_theta = np.tan(np.clip(theta, 0.0, np.pi / 2 - 0.001))
    p = ref / (ref.sum() + EPS) + EPS
    q = cube / (cube.sum(axis=0, keepdims=True) + EPS) + EPS
    p3 = p[:, np.newaxis, np.newaxis]
    sid = np.sum(p3 * np.log(p3 / q), axis=0) + np.sum(q * np.log(q / p3), axis=0)
    score = np.abs(sid) * tan_theta
    return np.clip(np.exp(-score), 0.0, 1.0)


def match_ed(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    diff = cube - ref[:, np.newaxis, np.newaxis]
    n = cube.shape[0]
    dist = np.sqrt(np.sum(diff ** 2, axis=0)) / np.sqrt(n)
    return 1.0 / (1.0 + dist)


def match_sad(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    diff = np.abs(cube - ref[:, np.newaxis, np.newaxis])
    sad = np.mean(diff, axis=0)
    return 1.0 / (1.0 + sad)


def match_sca(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    mu_ref = ref.mean()
    mu_cube = cube.mean(axis=0)
    ref_c = ref - mu_ref
    cube_c = cube - mu_cube[np.newaxis]
    cov = np.sum(ref_c[:, np.newaxis, np.newaxis] * cube_c, axis=0)
    std_ref = (ref_c ** 2).sum() ** 0.5 + EPS
    std_cube = np.sqrt((cube_c ** 2).sum(axis=0)) + EPS
    r = cov / (std_ref * std_cube)
    return np.clip((r + 1.0) / 2.0, 0.0, 1.0)


def match_cr_sam(cube: np.ndarray, ref: np.ndarray,
                 wavelengths: np.ndarray) -> np.ndarray:
    ref_cr = continuum_remove(ref, wavelengths)
    cube_cr = continuum_remove_cube(cube, wavelengths)
    return match_sam(cube_cr, ref_cr)


def match_cr_ed(cube: np.ndarray, ref: np.ndarray,
                wavelengths: np.ndarray) -> np.ndarray:
    ref_cr = continuum_remove(ref, wavelengths)
    cube_cr = continuum_remove_cube(cube, wavelengths)
    return match_ed(cube_cr, ref_cr)


def match_gd1(cube: np.ndarray, ref: np.ndarray,
              wavelengths: np.ndarray) -> np.ndarray:
    dx = np.diff(wavelengths)
    ref_d = np.diff(ref) / np.where(dx > 0, dx, 1.0)
    cube_d = np.diff(cube, axis=0) / np.where(dx > 0, dx, 1.0)[:, np.newaxis, np.newaxis]
    return match_sam(cube_d, ref_d)


def match_gd2(cube: np.ndarray, ref: np.ndarray,
              wavelengths: np.ndarray) -> np.ndarray:
    dx = np.diff(wavelengths)
    ref_d1 = np.diff(ref) / np.where(dx > 0, dx, 1.0)
    cube_d1 = np.diff(cube, axis=0) / np.where(dx > 0, dx, 1.0)[:, np.newaxis, np.newaxis]
    dx2 = (dx[:-1] + dx[1:]) / 2.0
    ref_d2 = np.diff(ref_d1) / np.where(dx2 > 0, dx2, 1.0)
    cube_d2 = np.diff(cube_d1, axis=0) / np.where(dx2 > 0, dx2, 1.0)[:, np.newaxis, np.newaxis]
    return match_sam(cube_d2, ref_d2)


def match_xcorr(cube: np.ndarray, ref: np.ndarray,
                max_lag: int = 3) -> np.ndarray:
    cube = _safe_cube(cube)
    n, rows, cols = cube.shape
    mu_ref = ref.mean()
    mu_cube = cube.mean(axis=0)
    ref_c = ref - mu_ref
    cube_c = cube - mu_cube[np.newaxis]
    std_ref = np.sqrt((ref_c ** 2).mean()) + EPS
    std_cube = np.sqrt((cube_c ** 2).mean(axis=0)) + EPS

    best_ncc = np.full((rows, cols), -1.0)
    for lag in range(-max_lag, max_lag + 1):
        if lag == 0:
            i_ref = slice(None); i_pix = slice(None)
        elif lag > 0:
            i_ref = slice(None, n - lag); i_pix = slice(lag, None)
        else:
            i_ref = slice(-lag, None); i_pix = slice(None, n + lag)
        n_eff = len(range(n)[i_ref])
        if n_eff < 2:
            continue
        ncc = np.sum(ref_c[i_ref, np.newaxis, np.newaxis] * cube_c[i_pix],
                     axis=0) / (n_eff * std_ref * std_cube)
        best_ncc = np.maximum(best_ncc, ncc)
    return np.clip((best_ncc + 1.0) / 2.0, 0.0, 1.0)


def match_dtw(cube: np.ndarray, ref: np.ndarray,
              window: int = 3) -> np.ndarray:
    cube = _safe_cube(cube)
    n_bands, rows, cols = cube.shape
    n_pix = rows * cols
    flat = cube.reshape(n_bands, n_pix).T.astype(np.float32)
    m = len(ref)
    ref_f = ref.astype(np.float32)

    CHUNK = 512
    distances = np.empty(n_pix, dtype=np.float32)

    for start in range(0, n_pix, CHUNK):
        end = min(start + CHUNK, n_pix)
        chunk = flat[start:end]
        k = end - start

        INF = np.float32(1e30)
        prev = np.full((k, m + 1), INF, dtype=np.float32)
        prev[:, 0] = 0.0

        for i in range(1, n_bands + 1):
            curr = np.full((k, m + 1), INF, dtype=np.float32)
            j_lo = max(1, i - window)
            j_hi = min(m, i + window)
            for j in range(j_lo, j_hi + 1):
                cost = (chunk[:, i - 1] - ref_f[j - 1]) ** 2
                best = np.minimum(
                    np.minimum(prev[:, j], curr[:, j - 1]),
                    prev[:, j - 1]
                )
                curr[:, j] = cost + best
            prev = curr

        distances[start:end] = np.sqrt(prev[:, m] / max(n_bands, m))

    dist_map = distances.reshape(rows, cols)
    return np.clip(np.exp(-dist_map.astype(np.float64)), 0.0, 1.0)


def match_ssim(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    C1 = (0.01) ** 2
    C2 = (0.03) ** 2
    C3 = C2 / 2.0

    mu_r = ref.mean()
    mu_p = cube.mean(axis=0)
    sigma_r = ref.std() + EPS
    sigma_p = cube.std(axis=0) + EPS
    sigma_rp = np.mean(
        (ref - mu_r)[:, np.newaxis, np.newaxis] * (cube - mu_p[np.newaxis]),
        axis=0
    )

    luminance = (2 * mu_r * mu_p + C1) / (mu_r ** 2 + mu_p ** 2 + C1)
    contrast  = (2 * sigma_r * sigma_p + C2) / (sigma_r ** 2 + sigma_p ** 2 + C2)
    structure = (sigma_rp + C3) / (sigma_r * sigma_p + C3)
    ssim = luminance * contrast * structure
    return np.clip((ssim + 1.0) / 2.0, 0.0, 1.0)


def match_jsd(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    p = (ref / (ref.sum() + EPS) + EPS)[:, np.newaxis, np.newaxis]
    q = cube / (cube.sum(axis=0, keepdims=True) + EPS) + EPS
    m = (p + q) / 2.0
    kld_pm = np.sum(p * np.log(p / m), axis=0)
    kld_qm = np.sum(q * np.log(q / m), axis=0)
    jsd = (kld_pm + kld_qm) / 2.0
    return np.clip(np.exp(-jsd), 0.0, 1.0)


def match_bhatt(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    p = np.maximum(ref / (ref.sum() + EPS), 0.0)
    q = np.maximum(cube / (cube.sum(axis=0, keepdims=True) + EPS), 0.0)
    bc = np.sum(np.sqrt(p[:, np.newaxis, np.newaxis] * q), axis=0)
    return np.clip(bc, 0.0, 1.0)


def match_mtf(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    n, rows, cols = cube.shape
    flat = cube.reshape(n, rows * cols)
    mu = flat.mean(axis=1)
    centered = flat - mu[:, np.newaxis]
    cov = (centered @ centered.T) / (flat.shape[1] - 1)
    cov += np.eye(n) * 1e-6 * np.trace(cov) / n

    d = ref - mu
    try:
        cov_inv = np.linalg.pinv(cov)
    except np.linalg.LinAlgError:
        gs.warning("MTF: covariance matrix is singular; falling back to SAM")
        return match_sam(cube, ref)

    w = cov_inv @ d
    denom = d @ w
    if abs(denom) < EPS:
        return np.zeros((rows, cols), dtype=np.float32)
    w /= denom

    scores = (w @ centered).reshape(rows, cols)
    lo, hi = np.percentile(scores, [1, 99])
    rng = hi - lo if hi > lo else 1.0
    return np.clip((scores - lo) / rng, 0.0, 1.0)


def match_cem(cube: np.ndarray, ref: np.ndarray) -> np.ndarray:
    cube = _safe_cube(cube)
    n, rows, cols = cube.shape
    flat = cube.reshape(n, rows * cols)
    cov = (flat @ flat.T) / flat.shape[1]
    cov += np.eye(n) * 1e-6 * np.trace(cov) / n

    try:
        cov_inv = np.linalg.pinv(cov)
    except np.linalg.LinAlgError:
        gs.warning("CEM: covariance matrix is singular; falling back to SAM")
        return match_sam(cube, ref)

    d = ref.astype(np.float64)
    w = cov_inv @ d
    denom = d @ w
    if abs(denom) < EPS:
        return np.zeros((rows, cols), dtype=np.float32)
    w /= denom

    scores = (w @ flat).reshape(rows, cols)
    lo, hi = np.percentile(scores, [1, 99])
    rng = hi - lo if hi > lo else 1.0
    return np.clip((scores - lo) / rng, 0.0, 1.0)


def match_ensemble(score_maps: dict[str, np.ndarray]) -> np.ndarray:
    maps = {k: v for k, v in score_maps.items() if k != 'ensemble'}
    if not maps:
        gs.fatal("Ensemble requires at least one other method to be computed")
    rows, cols = next(iter(maps.values())).shape
    n_pix = rows * cols
    rank_sum = np.zeros(n_pix, dtype=np.float64)
    for name, smap in maps.items():
        flat = smap.ravel().astype(np.float64)
        order = np.argsort(flat)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(n_pix)
        rank_sum += ranks
    rank_sum /= (len(maps) * n_pix)
    return rank_sum.reshape(rows, cols)

# ---------------------------------------------------------------------------
# Consensus hotspot analysis
# ---------------------------------------------------------------------------


def _norm_ppf(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, EPS, 1.0 - EPS)
    sign = np.where(p < 0.5, -1.0, 1.0)
    t = np.sqrt(-2.0 * np.log(np.where(p < 0.5, p, 1.0 - p)))
    num = 2.515517 + 0.802853 * t + 0.010328 * t ** 2
    den = 1.0 + 1.432788 * t + 0.189269 * t ** 2 + 0.001308 * t ** 3
    return sign * (t - num / den)


def _norm_cdf(z: np.ndarray) -> np.ndarray:
    return 0.5 * (1.0 + np.erf(z / np.sqrt(2.0)))


def _chi2_sf(x: np.ndarray, df: int) -> np.ndarray:
    try:
        from scipy.special import chdtrc
        return np.asarray(chdtrc(df, x), dtype=np.float64)
    except ImportError:
        z = (x - df) / np.sqrt(2.0 * df)
        return 1.0 - _norm_cdf(z)


def calibrate_scores(score_maps: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    calibrated: dict[str, np.ndarray] = {}
    for name, smap in score_maps.items():
        if name in ('ensemble', 'consensus'):
            continue
        flat = smap.ravel().astype(np.float64)
        n = len(flat)
        order = np.argsort(flat)
        ranks = np.empty(n, dtype=np.float64)
        ranks[order] = np.arange(1, n + 1, dtype=np.float64)
        calibrated[name] = (ranks / n).reshape(smap.shape)
    return calibrated


def compute_diversity_weights(calibrated: dict[str, np.ndarray]) -> dict[str, float]:
    names = list(calibrated.keys())
    k = len(names)
    if k == 1:
        return {names[0]: 1.0}

    flat = np.array([calibrated[nm].ravel() for nm in names], dtype=np.float64)
    centered = flat - flat.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centered, axis=1, keepdims=True) + EPS
    normed = centered / norms
    corr = normed @ normed.T

    np.fill_diagonal(corr, 0.0)
    mean_abs = np.abs(corr).sum(axis=1) / max(k - 1, 1)

    raw = 1.0 / (mean_abs + 0.15)
    raw *= k / raw.sum()
    return {nm: float(w) for nm, w in zip(names, raw)}


def fuse_probabilities(calibrated: dict[str, np.ndarray],
                       weights: dict[str, float],
                       mode: str) -> np.ndarray:
    names = [n for n in calibrated if n not in ('ensemble', 'consensus')]
    k = len(names)
    if k == 0:
        gs.fatal("consensus: no base method score maps available for fusion")

    rows, cols = next(c for n, c in calibrated.items()
                      if n not in ('ensemble', 'consensus')).shape
    P = np.array([calibrated[nm].ravel() for nm in names], dtype=np.float64)
    P = np.clip(P, EPS, 1.0 - EPS)
    W = np.array([weights.get(nm, 1.0) for nm in names], dtype=np.float64)
    W /= W.sum()

    if mode == 'rank_product':
        result = np.exp((W[:, np.newaxis] * np.log(P)).sum(axis=0))
    elif mode == 'fisher':
        Q = np.clip(1.0 - P, EPS, 1.0 - EPS)
        T = -2.0 * (W[:, np.newaxis] * np.log(Q)).sum(axis=0)
        eff_df = max(int(round(2.0 * W.sum() ** 2 / (W ** 2).sum())), 2)
        result = 1.0 - np.clip(_chi2_sf(T, eff_df), 0.0, 1.0)
    elif mode == 'stouffer':
        try:
            from scipy.special import ndtri, ndtr
            Z = ndtri(P)
            Z_comb = (W[:, np.newaxis] * Z).sum(axis=0) / np.sqrt((W ** 2).sum())
            result = ndtr(Z_comb)
        except ImportError:
            Z = _norm_ppf(P)
            Z_comb = (W[:, np.newaxis] * Z).sum(axis=0) / np.sqrt((W ** 2).sum())
            result = _norm_cdf(Z_comb)
    elif mode == 'group_product':
        name_to_idx = {nm: i for i, nm in enumerate(names)}
        active_groups: dict[str, list[int]] = {}
        for nm in names:
            placed = False
            for gname, gmembers in METHOD_GROUPS.items():
                if nm in gmembers:
                    active_groups.setdefault(gname, []).append(name_to_idx[nm])
                    placed = True
                    break
            if not placed:
                active_groups.setdefault('_other', []).append(name_to_idx[nm])
        group_results: list[np.ndarray] = []
        group_wts: list[float] = []
        for gname, idxs in active_groups.items():
            g_P = P[idxs]
            g_W = W[idxs]
            g_W = g_W / g_W.sum()
            log_g = np.log(np.clip(g_P, EPS, 1.0))
            g_score = np.exp((g_W[:, np.newaxis] * log_g).sum(axis=0))
            group_results.append(g_score)
            group_wts.append(float(len(idxs)))
        gw = np.array(group_wts)
        gw /= gw.sum()
        result = sum(w * s for w, s in zip(gw, group_results))
    elif mode == 'harmonic':
        result = 1.0 / (W[:, np.newaxis] / P).sum(axis=0)
    elif mode == 'min':
        result = P.min(axis=0)
    else:
        gs.fatal(f"Unknown fusion mode '{mode}'. Valid: {', '.join(FUSION_MODES)}")

    return np.clip(result.reshape(rows, cols), 0.0, 1.0)


def compute_consensus_stats(calibrated: dict[str, np.ndarray],
                            probability: np.ndarray,
                            agreement_threshold: float = 0.80) -> dict[str, np.ndarray]:
    names = [n for n in calibrated if n not in ('ensemble', 'consensus')]
    k = len(names)
    if k == 0:
        return {}

    rows, cols = probability.shape
    P = np.array([calibrated[nm].ravel() for nm in names], dtype=np.float64)

    above = (P > agreement_threshold).sum(axis=0) / k

    col_sum = P.sum(axis=0, keepdims=True) + EPS
    P_norm = P / col_sum
    log_P = np.where(P_norm > EPS, np.log(P_norm + EPS), 0.0)
    raw_h = -(P_norm * log_P).sum(axis=0)
    max_h = np.log(k) if k > 1 else 1.0
    entropy = raw_h / max_h

    spread = P.std(axis=0)

    prob_flat = probability.ravel()
    p75_prob = np.percentile(prob_flat, 75)
    p75_ent  = np.percentile(entropy, 75)
    conflict = ((prob_flat > p75_prob) & (entropy > p75_ent)).astype(np.float32)

    return {
        'agreement': above.reshape(rows, cols).astype(np.float32),
        'entropy':   entropy.reshape(rows, cols).astype(np.float32),
        'conflict':  conflict.reshape(rows, cols),
        'spread':    spread.reshape(rows, cols).astype(np.float32),
    }


def run_consensus_analysis(cube, ref_proc, wls, shift_win,
                           existing_score_maps,
                           fusion_mode='rank_product',
                           agreement_threshold=0.80,
                           skip_slow=False,
                           verbose=False) -> dict:
    SLOW = {'cr_sam', 'cr_ed', 'dtw'}
    score_maps = dict(existing_score_maps)
    to_run = [m for m in BASE_METHODS
              if m not in score_maps and (not skip_slow or m not in SLOW)]

    if verbose:
        gs.verbose(f"consensus: {len(to_run)} method(s) to compute: "
                   f"{', '.join(to_run) if to_run else '(none — all cached)'}")

    for idx, mth in enumerate(to_run):
        if verbose:
            gs.verbose(f"  [{idx + 1}/{len(to_run)}] {METHOD_LABELS.get(mth, mth)}")
        score_maps[mth] = compute_method(mth, cube, ref_proc, wls,
                                         shift_win, score_maps)

    if verbose:
        gs.verbose("consensus: calibrating scores (empirical CDF rank transform)")
    calibrated = calibrate_scores(score_maps)

    if len(calibrated) == 0:
        gs.fatal("consensus: no base-method score maps available for fusion")

    if verbose:
        gs.verbose(f"consensus: computing diversity weights "
                   f"({len(calibrated)} methods)")
    weights = compute_diversity_weights(calibrated)

    if verbose:
        for nm, w in sorted(weights.items(), key=lambda x: -x[1]):
            corr_note = " (down-weighted: correlated)" if w < 0.75 else ""
            gs.verbose(f"    {nm:<12s} w={w:.3f}{corr_note}")

    if verbose:
        gs.verbose(f"consensus: fusing with mode='{fusion_mode}'")
    probability = fuse_probabilities(calibrated, weights, fusion_mode)
    stats = compute_consensus_stats(calibrated, probability, agreement_threshold)

    p_flat = probability.ravel()
    n_pix = p_flat.size
    for thr in (0.9, 0.8, 0.7, 0.5):
        n_hot = int((p_flat >= thr).sum())
        if n_hot > 0:
            gs.verbose(f"consensus: {n_hot} pixels ({100.*n_hot/n_pix:.2f}%) ≥ {thr:.1f}")
            break

    return {
        'probability': probability,
        'agreement':   stats.get('agreement'),
        'entropy':     stats.get('entropy'),
        'conflict':    stats.get('conflict'),
        'spread':      stats.get('spread'),
        'weights':     weights,
        'score_maps':  score_maps,
        'calibrated':  calibrated,
    }

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


def compute_method(method: str, cube: np.ndarray, ref: np.ndarray,
                   wavelengths: np.ndarray,
                   shift_window: int,
                   score_maps: dict) -> np.ndarray:
    if method == 'sam':
        return match_sam(cube, ref)
    if method == 'sid':
        return match_sid(cube, ref)
    if method == 'sid_sam':
        return match_sid_sam(cube, ref)
    if method == 'ed':
        return match_ed(cube, ref)
    if method == 'sad':
        return match_sad(cube, ref)
    if method == 'sca':
        return match_sca(cube, ref)
    if method == 'cr_sam':
        return match_cr_sam(cube, ref, wavelengths)
    if method == 'cr_ed':
        return match_cr_ed(cube, ref, wavelengths)
    if method == 'gd1':
        return match_gd1(cube, ref, wavelengths)
    if method == 'gd2':
        return match_gd2(cube, ref, wavelengths)
    if method == 'xcorr':
        return match_xcorr(cube, ref, max_lag=shift_window)
    if method == 'dtw':
        return match_dtw(cube, ref, window=shift_window)
    if method == 'ssim':
        return match_ssim(cube, ref)
    if method == 'jsd':
        return match_jsd(cube, ref)
    if method == 'bhatt':
        return match_bhatt(cube, ref)
    if method == 'mtf':
        return match_mtf(cube, ref)
    if method == 'cem':
        return match_cem(cube, ref)
    if method == 'ensemble':
        return match_ensemble(score_maps)
    gs.fatal(f"Unknown method: {method}")

# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def set_similarity_colors(raster: str) -> None:
    try:
        gs.run_command('r.colors', map=raster, color='bcyr', quiet=True)
    except Exception:
        pass


def set_raster_metadata(raster: str, source: str,
                        method: str, description: str) -> None:
    try:
        gs.run_command('r.support', map=raster,
                       title='Papaver somniferum similarity map',
                       history=f'i.hyper.poppy: method={method}; {description}',
                       description=description, quiet=True)
    except Exception:
        pass


def point_analysis(pixel_spec: np.ndarray, ref: np.ndarray,
                   wavelengths: np.ndarray,
                   methods: list[str],
                   shift_window: int,
                   flag_c: bool,
                   flag_z: bool = False) -> dict[str, float]:
    cube1 = pixel_spec[:, np.newaxis, np.newaxis].astype(np.float64)
    if flag_c:
        cube1 = continuum_remove_cube(cube1, wavelengths)
    scores: dict[str, float] = {}
    for mth in methods:
        if mth in ('ensemble', 'consensus'):
            continue
        smap = compute_method(mth, cube1, ref, wavelengths, shift_window, scores)
        scores[mth] = float(smap[0, 0])
    if 'ensemble' in methods and len(scores) >= 1:
        scores['ensemble'] = float(match_ensemble(
            {k: np.array([[v]]) for k, v in scores.items()}
        )[0, 0])
    return scores

# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------


def main(options: dict, flags: dict) -> int:
    """Entry point called by the GRASS parser after ``gs.parser()`` returns.

    Orchestrates the Papaver somniferum detection pipeline:

    1. Load built-in poppy spectra for selected growth stage(s).
    2. Optionally append user-supplied spectra from poppy_db=.
    3. Read band metadata from the 3D raster.
    4. Load the hyperspectral cube.
    5. For each spectrum: resample, preprocess, score every pixel.
    6. Max-pool across all spectra → single abundance map (best match wins).
    7. Write output; optionally write per-stage maps.
    """
    raster3d    = options['input']
    output      = options['output']
    growth_stg  = options.get('growth_stage', 'all')
    poppy_db    = options.get('poppy_db') or ''
    methods_str = options.get('method', 'sam')
    out_prefix  = options.get('output_prefix') or ''
    resample_m  = options.get('resample', 'linear')
    normalize_m = options.get('normalize', 'none')
    shift_win   = int(options.get('shift_window', 3))
    fusion_mode = options.get('fusion_mode', 'rank_product')
    agr_thresh  = float(options.get('agreement_threshold', 0.80))
    min_wl      = float(options['min_wavelength']) if options.get('min_wavelength') else None
    max_wl      = float(options['max_wavelength']) if options.get('max_wavelength') else None
    coords_str  = options.get('coordinates') or ''

    flag_n = flags.get('n', False)
    flag_i = flags.get('i', False)
    flag_v = flags.get('v', False)
    flag_c = flags.get('c', False)
    flag_p = flags.get('p', False)
    flag_z = flags.get('z', False)

    methods = [m.strip() for m in methods_str.split(',') if m.strip()]
    unknown = [m for m in methods if m not in ALL_METHODS]
    if unknown:
        gs.fatal(f"Unknown method(s): {', '.join(unknown)}. "
                 f"Valid: {', '.join(ALL_METHODS)}")
    if fusion_mode not in FUSION_MODES:
        gs.fatal(f"Unknown fusion_mode '{fusion_mode}'. "
                 f"Valid: {', '.join(FUSION_MODES)}")

    do_consensus = 'consensus' in methods
    if do_consensus and 'ensemble' in methods:
        gs.warning("Both 'consensus' and 'ensemble' requested; "
                   "'ensemble' will be skipped — consensus subsumes it.")
        methods = [m for m in methods if m != 'ensemble']

    if flag_p and not coords_str:
        gs.fatal("Point mode (-p) requires coordinates= to be set")
    if 'ensemble' in methods and len(methods) < 2:
        gs.fatal("ensemble requires at least one additional method")

    # ------------------------------------------------------------------
    # Step 1: Assemble poppy spectra
    # ------------------------------------------------------------------
    all_spectra = get_builtin_spectra(growth_stg)
    if poppy_db:
        extra = load_poppy_db(poppy_db)
        all_spectra.extend(extra)
        gs.message(f"Appended {len(extra)} external spectrum/spectra from {poppy_db}")

    n_spectra = len(all_spectra)
    gs.message(f"Poppy reference spectra: {n_spectra} "
               f"({growth_stg}" +
               (f" + {len(all_spectra) - len(get_builtin_spectra(growth_stg))} external"
                if poppy_db else "") + ")")
    for sp in all_spectra:
        wl0, wl1 = sp['wavelengths'][0], sp['wavelengths'][-1]
        gs.message(f"  {sp['stage']:20s}  {wl0:.0f}–{wl1:.0f} nm  "
                   f"({sp['citation']})")

    # ------------------------------------------------------------------
    # Step 2: Band metadata
    # ------------------------------------------------------------------
    gs.message(f"Scanning hyperspectral bands in: {raster3d}")
    bands = get_band_info(raster3d, only_valid=flag_n, min_wl=min_wl, max_wl=max_wl)
    wls = np.array([b['wavelength'] for b in bands], dtype=np.float64)
    gs.message(f"Found {len(bands)} usable bands: {wls[0]:.1f} – {wls[-1]:.1f} nm")

    if flag_i:
        gs.message("Built-in poppy spectrum key wavelengths:")
        jia_notes = {
            'before_flowering': '438, 528, 736, 754, 1207 nm (Jia et al. 2011, JM≈1.20)',
            'flowering':        '405, 424, 524, 760, 1974 nm (Jia et al. 2011, JM≈1.22)',
            'harvesting':       '468, 726, 746, 982, 1220, 1689 nm (Jia et al. 2011, JM≈1.21)',
            'leaf_healthy':     '550, 670 nm green/red wells, 710–730 nm red-edge (Calderón 2014)',
        }
        for stage, note in jia_notes.items():
            gs.message(f"  {stage:20s}  {note}")
        return 0

    # ------------------------------------------------------------------
    # Step 3: Load full cube (shared across all spectra)
    # ------------------------------------------------------------------
    gs.message("Loading hyperspectral cube into memory...")
    cube = load_cube(bands, raster3d, verbose=flag_v)

    if flag_z:
        cube = to_prob_simplex_cube(cube)
    if normalize_m != 'none':
        cube = normalize_cube(cube, normalize_m)
    if flag_c and 'cr_sam' not in methods and 'cr_ed' not in methods:
        gs.message("Applying continuum removal to all pixels...")
        cube = continuum_remove_cube(cube, wls)

    cube = cube.astype(np.float64)

    abundance_map = np.full(cube.shape[1:], -np.inf, dtype=np.float64)
    best_stage_map = np.zeros(cube.shape[1:], dtype=np.int32)

    if flag_p:
        east, north = (float(v) for v in coords_str.split(','))
        gs.message(f"Point mode: extracting spectrum at E={east}, N={north}")
        spec_pixel = read_pixel_spectrum(bands, raster3d, east, north)
        if np.all(np.isnan(spec_pixel)):
            gs.fatal("All band values are null at the given coordinates")
        if flag_z:
            spec_pixel = to_prob_simplex(spec_pixel)
        if normalize_m != 'none':
            spec_pixel = normalize_spectrum(spec_pixel, normalize_m)
    else:
        spec_pixel = None

    # ------------------------------------------------------------------
    # Step 4: Score every pixel against every spectrum
    # ------------------------------------------------------------------
    gs.message(f"Running poppy search: {n_spectra} spectrum/spectra × "
               f"{len(methods)} method(s)...")

    primary_method = ('ensemble' if 'ensemble' in methods
                      else ('consensus' if do_consensus else methods[0]))

    per_stage_maps: dict[str, np.ndarray] = {}

    for sp_idx, spectrum in enumerate(all_spectra):
        ref_wls  = spectrum['wavelengths']
        ref_vals = spectrum['reflectances']
        stage    = spectrum['stage']
        sp_label = f"{stage}"

        gs.percent(sp_idx, n_spectra, 5)

        lut = WavelengthLUT(ref_wls, wls, fill='edge')
        if not lut.has_overlap:
            if flag_v:
                gs.verbose(f"  Spectrum '{stage}': no overlap with sensor — skipped")
            continue

        if flag_v:
            gs.verbose(f"  [{sp_idx+1}/{n_spectra}] {stage}: {lut.coverage_report()}")

        ref_resampled = resample_reference(ref_wls, ref_vals, wls,
                                           method=resample_m, lut=lut)
        ref_proc = ref_resampled.copy()
        if flag_z:
            ref_proc = to_prob_simplex(ref_proc)
        if normalize_m != 'none':
            ref_proc = normalize_spectrum(ref_proc, normalize_m)
        if flag_c and 'cr_sam' not in methods and 'cr_ed' not in methods:
            ref_proc = continuum_remove(ref_proc, wls)

        if flag_p:
            scores = point_analysis(spec_pixel, ref_proc, wls, methods,
                                    shift_win, flag_c, flag_z=False)
            score_val = scores.get(primary_method, float('nan'))
            gs.message(f"  [{sp_idx+1:2d}/{n_spectra}] {sp_label:20s}  "
                       f"{primary_method}={score_val:.4f}")
            continue

        if do_consensus:
            consensus_result = run_consensus_analysis(
                cube, ref_proc, wls, shift_win,
                existing_score_maps={},
                fusion_mode=fusion_mode,
                agreement_threshold=agr_thresh,
                skip_slow=False,
                verbose=flag_v,
            )
            sp_map = consensus_result['probability']
        else:
            score_maps: dict[str, np.ndarray] = {}
            non_ensemble = [m for m in methods if m != 'ensemble']
            do_ensemble  = 'ensemble' in methods
            for mth in non_ensemble:
                score_maps[mth] = compute_method(mth, cube, ref_proc, wls,
                                                 shift_win, score_maps)
            if do_ensemble:
                score_maps['ensemble'] = match_ensemble(score_maps)
            sp_map = score_maps['ensemble' if do_ensemble else non_ensemble[0]]

        better = sp_map > abundance_map
        abundance_map = np.where(better, sp_map, abundance_map)
        best_stage_map = np.where(better, sp_idx, best_stage_map)
        per_stage_maps[stage] = sp_map

        if flag_v:
            valid = sp_map[np.isfinite(sp_map)]
            if valid.size > 0:
                gs.verbose(f"  [{sp_idx+1:2d}] {stage:20s}  "
                           f"max={valid.max():.4f} mean={valid.mean():.4f}")

    gs.percent(n_spectra, n_spectra, 5)

    if flag_p:
        sep = "-" * 64
        gs.message(sep)
        gs.message(f"i.hyper.poppy — Point analysis  E={coords_str.split(',')[0]}  "
                   f"N={coords_str.split(',')[1]}")
        gs.message(f"  (scores shown above per growth stage, method={primary_method})")
        gs.message(sep)
        return 0

    abundance_map = np.where(np.isfinite(abundance_map), abundance_map, np.nan)

    # ------------------------------------------------------------------
    # Step 5: Write output
    # ------------------------------------------------------------------
    stage_desc = (growth_stg if growth_stg != 'all'
                  else ', '.join(sp['stage'] for sp in all_spectra))
    ref_desc = f"Papaver somniferum ({n_spectra} spectra, stage={stage_desc})"

    gs.message(f"Writing max-abundance map: {output}")
    write_raster(abundance_map, output)
    set_similarity_colors(output)
    set_raster_metadata(output, raster3d, primary_method, ref_desc)

    if out_prefix:
        # Per-stage score maps
        for stage, smap in per_stage_maps.items():
            name = f"{out_prefix}_{stage}"
            gs.message(f"  Writing stage map: {name}")
            write_raster(smap, name)
            set_similarity_colors(name)
            try:
                gs.run_command('r.support', map=name,
                               title=f'Poppy similarity — {stage}',
                               description=f'method={primary_method}; {stage}',
                               quiet=True)
            except Exception:
                pass

        # Best-stage index map
        idx_name = f"{out_prefix}_best_stage_idx"
        gs.message(f"  Writing best-stage index map: {idx_name}")
        write_raster(best_stage_map.astype(np.float64), idx_name)
        try:
            stage_labels = [sp['stage'] for sp in all_spectra]
            gs.run_command('r.colors', map=idx_name, color='random', quiet=True)
            cats = '\n'.join(f"{i}:{s}" for i, s in enumerate(stage_labels))
            gs.run_command('r.support', map=idx_name,
                           title='Best-matching poppy growth stage index',
                           description=f'0={stage_labels[0]}, '
                                       + ', '.join(f'{i}={s}' for i, s
                                                   in enumerate(stage_labels)),
                           quiet=True)
        except Exception:
            pass

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    s = abundance_map[np.isfinite(abundance_map)]
    sep = "=" * 60
    gs.message(" ")
    gs.message(sep)
    gs.message("i.hyper.poppy completed successfully.")
    gs.message(f"  Spectra matched  : {n_spectra} (stage={growth_stg})")
    gs.message(f"  Method           : {primary_method}")
    gs.message(f"  Output map       : {output}")
    if s.size > 0:
        gs.message(f"  Score range      : [{s.min():.4f}, {s.max():.4f}]")
        gs.message(f"  Score mean       : {s.mean():.4f}")
        for thr in (0.9, 0.8, 0.7, 0.5):
            n_above = int((s >= thr).sum())
            pct = 100.0 * n_above / s.size
            if n_above > 0:
                gs.message(f"  Pixels >= {thr:.1f}    : {n_above} ({pct:.2f}%)")
                break
    if out_prefix:
        gs.message(f"  Stage maps prefix: {out_prefix}_<stage>")
        gs.message(f"  Best-stage index : {out_prefix}_best_stage_idx")
    gs.message("")
    gs.message("  Key spectral references:")
    gs.message("    Jia et al. (2011) IEEE TGRS 49(9):3414–3422")
    gs.message("      — canopy-level field spectrometry, 3 growth stages")
    gs.message("    Calderón et al. (2014) Precision Agric.")
    gs.message("      — leaf-level integrating-sphere, healthy vs DM-infected")
    gs.message("    Wang (2013) UNB Tech. Rep. 286 — EO-1 Hyperion, Afghanistan")
    gs.message("    Bennington (2008) Cranfield PhD — IKONOS multi-temporal, Afghanistan")
    gs.message(sep)

    return 0


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main(options, flags))
