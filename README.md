# i.hyper.poppy

**GRASS GIS module — Opium poppy (*Papaver somniferum*) detection from
hyperspectral imagery using built-in multi-stage spectral references**

Part of the [i.hyper](../README.md) module family for VNIR-SWIR hyperspectral
data processing in GRASS GIS.

---

## Overview

`i.hyper.poppy` finds pixels in a hyperspectral 3D raster that best match
*Papaver somniferum* (opium poppy) canopy or leaf spectra.  For each pixel a
similarity score in **[0, 1]** (0 = no match, 1 = perfect match) is computed
using one or more of **19 similarity methods** drawn from remote sensing,
signal analysis, information theory, morphological mathematics, and
subpixel-detection theory.

Unlike a generic spectral matcher, the module ships **four built-in reference
spectra** compiled from the peer-reviewed literature so no external spectral
library is required:

| Stage | Sensor / scale | Optimal bands (nm) | JM dist. | Source |
|---|---|---|---|---|
| `before_flowering` | ASD FieldSpec Pro, canopy | 438, 528, 736, 754, 1207 | ≈ 1.20 | Jia et al. 2011 |
| `flowering` | ASD FieldSpec Pro, canopy | 405, 424, 524, 760, 1974 | ≈ 1.22 | Jia et al. 2011 |
| `harvesting` | ASD FieldSpec Pro, canopy | 468, 726, 746, 982, 1220, 1689 | ≈ 1.21 | Jia et al. 2011 |
| `leaf_healthy` | Li-Cor integrating sphere, leaf | 550 (green peak), 670 (red well), 710–730 (red-edge) | — | Calderón et al. 2014 |

The flowering stage provides the strongest spectral contrast: white petals raise
visible-range reflectance by ~2.7× relative to before-flowering, making
400–530 nm the most discriminative region during bloom.

The module scores every pixel against all selected stage spectra and outputs a
**max-abundance map** — the best match score across all stages.  With
`output_prefix=` it also writes one per-stage similarity map plus a
best-matching-stage index raster.

A user-supplied spectrum CSV or JSON can be appended to the built-in set via
`poppy_db=`.

The optional **multi-method consensus pipeline** (`method=consensus`) runs all
17 base methods simultaneously, calibrates their scores to per-pixel
probabilities via empirical CDF rank transform, down-weights correlated methods
using inter-method diversity analysis, and fuses the result into a single
hotspot probability map — plus four diagnostic maps (agreement, entropy,
conflict, spread).

Input 3D rasters are produced by
[`i.hyper.import`](../i.hyper.import) or
[`i.hyper.atcorr`](../i.hyper.atcorr).

---

## Built-in spectra

### Canopy spectra — Jia et al. (2011)

Field spectrometry study conducted on an official poppy farm in NW China
(~1900 m altitude, arid valley).  Measurements with an ASD FieldSpec Pro
portable spectrometer (350–2500 nm, 1.4 nm sampling) taken at canopy nadir
from 1.5 m height, averaged over 300 spectra per growth period.

Atmospheric water-absorption bands excluded from analysis:
1360–1420 nm, 1801–1965 nm, 2451–2500 nm.

The three growth periods studied:

| Code | Date | Height | Competing species |
|---|---|---|---|
| P_BF (before flowering) | 2008-06-18 | 65 cm | Wheat, barley, alfalfa |
| P_FL (flowering) | 2008-07-05 | 125 cm | Wheat, barley, alfalfa |
| P_HV (harvesting) | 2008-07-14 | 130 cm | Wheat, barley, alfalfa |

Optimal waveband selection used a three-level synthetic method combining
Mann–Whitney U-test, Jeffries–Matusita distance, and CART classification trees.
The selected wavebands achieved **100 % poppy discrimination accuracy** (PA)
across all three growth periods at canopy level (Table IV, Jia et al. 2011).

### Leaf spectrum — Calderón et al. (2014)

Leaf-level reflectance and transmittance measured with a Li-Cor 1800-12
integrating sphere (350–1000 nm, 0.5 nm sampling) on asymptomatic (healthy)
*P. somniferum* cv. Nigrum leaves.  Useful for matching high-spatial-resolution
airborne or UAV hyperspectral data where individual leaf signatures are resolved.

---

## Similarity methods

| Key | Name | Category |
|-----|------|----------|
| `sam` | Spectral Angle Mapper | Geometric |
| `sid` | Spectral Information Divergence | Information theory |
| `sid_sam` | SID × tan(SAM) hybrid | Combined |
| `ed` | Euclidean Distance (L2) | Distance |
| `sad` | Spectral Absolute Difference (L1) | Distance |
| `sca` | Spectral Correlation Angle (Pearson *r*) | Statistical |
| `cr_sam` | Continuum-Removed SAM | Morphological |
| `cr_ed` | Continuum-Removed Euclidean Distance | Morphological |
| `gd1` | 1st-Derivative Shape Matching | Signal analysis |
| `gd2` | 2nd-Derivative Shape Matching | Signal analysis |
| `xcorr` | Normalized Cross-Correlation | Signal analysis |
| `dtw` | Dynamic Time Warping | Signal analysis |
| `ssim` | Spectral Structural Similarity Index | Signal analysis |
| `jsd` | Jensen-Shannon Divergence | Information theory |
| `bhatt` | Bhattacharyya Coefficient | Statistical |
| `mtf` | Matched Tuned Filter | Subpixel detection |
| `cem` | Constrained Energy Minimization | Subpixel detection |
| `ensemble` | Rank-based Borda-count fusion | Meta |
| `consensus` | Multi-method calibrated probability fusion | Meta |

---

## Consensus analysis (`method=consensus`)

When `consensus` is requested the module executes a four-step pipeline:

1. **Compute all base methods** — all 17 methods except `ensemble`/`consensus`
   are run against the cube.

2. **Empirical CDF calibration** — each score map is rank-transformed to a
   uniform probability in (0, 1].  This removes the scale bias that would
   otherwise allow high-range methods (SAM in [0.8–1.0]) to dominate
   over low-range methods (ED in [0.0–0.2]).

3. **Diversity weighting** — the full *k × k* Pearson correlation matrix is
   built across calibrated maps.  Each method receives weight
   `w ∝ 1 / mean |r_ij|`, so correlated clusters (e.g. SAM + SCA) do not
   over-count their shared evidence.

4. **Fusion** — six modes available (controlled by `fusion_mode=`):

| Mode | Description |
|------|-------------|
| `rank_product` | Weighted geometric mean of rank fractions *(default)* |
| `fisher` | Fisher χ² combined probability test |
| `stouffer` | Stouffer weighted Z-score combination |
| `group_product` | AND within method-type groups, OR across groups |
| `harmonic` | Harmonic mean — strictest, all methods must agree |
| `min` | Minimum across methods — absolute strictest |

### Consensus output maps (requires `output_prefix=`)

| Map | Content |
|-----|---------|
| `{prefix}_consensus_agreement` | Fraction of methods voting above threshold [0, 1] |
| `{prefix}_consensus_entropy` | Method agreement entropy: 0 = unanimous, 1 = maximal conflict |
| `{prefix}_consensus_conflict` | High-probability pixels where methods disagree |
| `{prefix}_consensus_spread` | Std dev of calibrated scores per pixel |

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input=` | — | Input hyperspectral 3D raster (from `i.hyper.import` / `i.hyper.atcorr`) |
| `output=` | — | Output max-abundance raster (FCELL, [0, 1]) |
| `growth_stage=` | `all` | Growth stage(s): `before_flowering`, `flowering`, `harvesting`, `leaf_healthy`, or `all` |
| `poppy_db=` | — | Optional supplementary spectrum CSV or JSON file (appended to built-in set) |
| `method=` | `sam` | Comma-separated similarity method(s) — see table above |
| `fusion_mode=` | `rank_product` | Fusion strategy for `method=consensus` |
| `agreement_threshold=` | `0.80` | Calibrated-probability threshold for per-pixel agreement count (consensus only) |
| `output_prefix=` | — | Prefix for per-stage similarity maps and best-stage index map |
| `resample=` | `linear` | Interpolation for resampling reference to sensor grid: `linear`, `cubic`, `pchip` |
| `normalize=` | `none` | Spectrum normalisation before matching: `none`, `area`, `max`, `minmax`, `vector` |
| `shift_window=` | `3` | Maximum band-shift for shift-tolerant methods (`xcorr`, `dtw`); 0 disables |
| `min_wavelength=` | — | Lower wavelength limit (nm) |
| `max_wavelength=` | — | Upper wavelength limit (nm) |
| `coordinates=` | — | `east,north` for point-mode analysis (requires `-p`) |

## Flags

| Flag | Effect |
|------|--------|
| `-n` | Only use bands marked `valid=1` in metadata |
| `-i` | Info mode: print built-in spectrum summary and key bands, then exit |
| `-v` | Verbose: show per-stage LUT coverage and score range |
| `-c` | Apply convex-hull continuum removal before matching |
| `-p` | Point mode: print score table for one pixel at `coordinates=` |
| `-z` | Normalize spectra to probability simplex (sum-to-one); recommended with `sid`, `jsd`, `bhatt` |

---

## Output

| Output | GRASS type | Range | Content |
|--------|-----------|-------|---------|
| `output=` | FCELL | 0–1 | Max similarity across all stages |
| `{prefix}_before_flowering` | FCELL | 0–1 | Before-flowering stage similarity |
| `{prefix}_flowering` | FCELL | 0–1 | Flowering stage similarity |
| `{prefix}_harvesting` | FCELL | 0–1 | Harvesting stage similarity |
| `{prefix}_leaf_healthy` | FCELL | 0–1 | Healthy-leaf similarity |
| `{prefix}_best_stage_idx` | FCELL | 0–3 | Index of best-matching stage |

All output maps share a blue → yellow → red colour ramp (0 = blue, 1 = red).

---

## Quick examples

```bash
# Default: match all four built-in stages with SAM
i.hyper.poppy input=scene_atcorr output=poppy_abundance

# Flowering stage only — most discriminative period
i.hyper.poppy input=scene_atcorr output=poppy_fl \
  growth_stage=flowering method=sam

# All stages, ensemble of six methods, per-stage maps
i.hyper.poppy input=scene_atcorr output=poppy_best \
  growth_stage=all \
  method=sam,cr_sam,gd1,jsd,bhatt,ensemble \
  output_prefix=poppy

# Full consensus with group-product fusion
i.hyper.poppy input=scene_atcorr output=poppy_consensus \
  growth_stage=all method=consensus fusion_mode=group_product \
  output_prefix=poppy_cons

# Point inspection at one pixel — compare all stages
i.hyper.poppy input=scene_atcorr output=_ \
  growth_stage=all method=sam,sid,bhatt,dtw \
  coordinates="452300,4325100" -p -v -i

# Add your own field-measured spectrum alongside the built-in ones
i.hyper.poppy input=scene_atcorr output=poppy_custom \
  growth_stage=all poppy_db=my_field_spectra.csv \
  method=sam,cr_sam,ensemble output_prefix=poppy_custom

# Restrict to visible + red-edge where flowering is most discriminative
i.hyper.poppy input=scene_atcorr output=poppy_vis \
  growth_stage=flowering \
  min_wavelength=400 max_wavelength=800 \
  method=sam,sca,gd1,ensemble -c
```

---

## Performance notes

| Method class | Speed | Note |
|---|---|---|
| `sam`, `sid`, `ed`, `sad`, `sca`, `jsd`, `bhatt`, `ssim`, `xcorr`, `gd1`, `gd2` | Fast | Fully vectorized over all pixels |
| `mtf`, `cem` | Fast | One k×k covariance inversion, then linear |
| `dtw` | Moderate | Chunked Sakoe-Chiba rolling-window |
| `cr_sam`, `cr_ed` | Slow | Per-pixel Graham-scan convex hull |
| `consensus` | Slow (first run) | Runs all 17 base methods |

With `growth_stage=all` the module runs the selected method(s) four times
(once per spectrum) and max-pools the results.  Total time ≈ 4× a single-stage
run; negligible for all fast methods.

---

## References

- **Jia K., Wu B., Tian Y., Li Q. & Du X. (2011)**. Spectral discrimination of
  opium poppy using field spectrometry. *IEEE Trans. Geosci. Remote Sens.*
  49(9): 3414–3422. DOI: 10.1109/TGRS.2011.2126582
  — *Source for before_flowering, flowering, harvesting canopy spectra and
  optimal band selection.*

- **Calderón R., Montes-Borrego M., Landa B.B., Navas-Cortés J.A. &
  Zarco-Tejada P.J. (2014)**. Detection of downy mildew of opium poppy using
  high-resolution multi-spectral and thermal imagery acquired with an unmanned
  aerial vehicle. *Precision Agric.* DOI: 10.1007/s11119-014-9360-y
  — *Source for healthy leaf integrating-sphere spectrum.*

- **Wang J. (2013)**. Unsupervised detection of opium poppy fields in
  Afghanistan from EO-1 Hyperion data. UNB Technical Report No. 286,
  University of New Brunswick, Fredericton.
  — *Hyperspectral endmember extraction and MESMA/MTMF methodology context.*

- **Bennington A.L. (2008)**. Application of multi-spectral remote sensing
  for crop discrimination in Afghanistan. PhD Thesis, Cranfield University.
  — *Multi-temporal spectral separability analysis, Afghanistan context.*

- **Iqbal F., Lucieer A. & Barry K. (2018)**. Poppy crop capsule volume
  estimation using UAS remote sensing and random forest regression.
  *Int. J. Appl. Earth Obs. Geoinf.* 73: 362–373.
  — *UAS-scale precision agriculture context for Tasmanian pharmaceutical poppy.*

### Similarity method references

- Kruse *et al.* (1993) — SAM. *Remote Sens. Environ.* 44, 145–163.
- Chang C.I. (2000) — SID and SID-SAM. *IEEE Trans. Inf. Theory* 46(5).
- Clark *et al.* (1987) — Continuum removal. *J. Geophys. Res.* 92(B12).
- Reed & Yu (1990) — Matched Tuned Filter. *IEEE Trans. ASSP* 38(10).
- Chang & Heinz (2000) — CEM. *IEEE Trans. GRSS* 38(3).
- Wang *et al.* (2004) — SSIM. *IEEE Trans. Image Process.* 13(4).
- Sakoe & Chiba (1978) — DTW. *IEEE Trans. ASSP* 26(1).
- Fisher R.A. (1932) — Combined probability test. Oliver & Boyd.
- Stouffer S.A. *et al.* (1949) — Measurement and Prediction. Princeton UP.

---

## See also

[`i.hyper.import`](../i.hyper.import) ·
[`i.hyper.atcorr`](../i.hyper.atcorr) ·
[`i.hyper.sleuth`](../i.hyper.sleuth) ·
[`i.hyper.cannabis`](../i.hyper.cannabis) ·
[`i.hyper.continuum`](../i.hyper.continuum) ·
[`i.hyper.spectroscopy`](../i.hyper.spectroscopy)
