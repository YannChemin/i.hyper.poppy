## DESCRIPTION

*i.hyper.poppy* finds the pixels in a hyperspectral 3D raster that most
closely match *Papaver somniferum* (opium poppy), using built-in canopy-level
and leaf-level reference spectra compiled from peer-reviewed literature. For
each pixel the similarity of every selected reference spectrum is computed and
the **maximum** is written to **output**, so the result is a max-abundance map
in **[0, 1]** (0 = no match, 1 = perfect match to the best-matching spectrum).

Reference spectra are resampled onto the sensor's own wavelength grid before
matching, and any preprocessing (continuum removal, normalization, simplex
projection) is applied identically to the reference and to every pixel.

### Built-in reference spectra

| **growth_stage** | Measurement | Optimal discrimination bands |
|---|---|---|
| **before_flowering** | Vegetative rosette canopy, June, NW China, 65 cm; ASD FieldSpec Pro 350-2500 nm | 438, 528, 736, 754, 1207 nm (JM ~ 1.20) |
| **flowering** | White-petal canopy, July, NW China, 125 cm; elevated visible reflectance from the large white petals | 405, 424, 524, 760, 1974 nm (JM ~ 1.22) |
| **harvesting** | Mature capsules and senescing leaves, mid-July, NW China, 130 cm | 468, 726, 746, 982, 1220, 1689 nm (JM ~ 1.21) |
| **leaf_healthy** | Asymptomatic leaf-level spectrum, 350-800 nm; Li-Cor 1800-12 integrating sphere | — |

**growth_stage=all** (the default) matches against all four. JM is the
Jeffries-Matusita separability of poppy against the surrounding crops in the
source study, not a property of the scene being analysed.

**poppy_db=** appends a supplementary spectrum to the built-in set, as a CSV
of `wavelength,reflectance` lines or as JSON.

### Similarity methods

**method=** selects one or more of 19 measures: **sam** (Spectral Angle
Mapper, the default), **sid**, **sid_sam**, **ed**, **sad**, **sca**,
**cr_sam**, **cr_ed**, **gd1**, **gd2**, **xcorr**, **dtw**, **ssim**,
**jsd**, **bhatt**, **mtf**, **cem**, **ensemble** and **consensus**. The
scale-invariant ones (**sam**, **sca**) are insensitive to illumination
magnitude; the divergence measures (**sid**, **jsd**, **bhatt**) expect
spectra on the probability simplex and should be combined with the **-z**
flag.

**method=consensus** runs all base methods, calibrates each score map to a
per-pixel probability through an empirical CDF rank transform, and fuses the
result according to **fusion_mode=**: **rank_product** (the default,
weighted geometric mean of rank fractions), **fisher** (chi-squared combined
probability), **stouffer** (weighted Z-score), **group_product**
(geometric mean within method groups, arithmetic mean across them),
**harmonic** (strict: every method must agree) or **min** (strictest).
**agreement_threshold=** sets the calibrated probability above which a method
counts as agreeing for the per-pixel agreement count.

### Parameters

- **input** — hyperspectral 3D raster, from *i.hyper.import* or
  *i.hyper.atcorr*
- **output** — max-abundance raster
- **output_prefix** — also write the individual per-stage and per-method
  similarity maps
- **resample** — **linear** (default), **cubic** or **pchip** interpolation of
  the reference onto the sensor bands
- **normalize** — **none** (default), **area**, **max**, **minmax** or
  **vector**
- **shift_window** — maximum band shift for the shift-tolerant methods
  (**xcorr**, **dtw**); 0 disables shift tolerance, default 3
- **min_wavelength**, **max_wavelength** — restrict matching to a wavelength
  range; both default to the sensor's own range
- **coordinates** — east,north for point mode

### Flags

- **-n** — use only bands marked `valid=1` in the band metadata
- **-i** — print band coverage and a summary of the built-in spectra, then
  exit
- **-v** — print processing details and per-band information
- **-c** — convex-hull continuum removal of reference and pixel spectra before
  matching
- **-p** — point mode: analyse the single pixel at **coordinates=** and print
  the full per-method score table
- **-z** — normalize reference and pixel spectra to unit sum (probability
  simplex); recommended with **sid**, **jsd** and **bhatt**

## NOTES

### Input requirements

Spectra must be surface reflectance. Apply atmospheric correction with
*i.hyper.atcorr* or *i.hyper.smac* first; matching radiance against
reflectance references is meaningless.

The module reads wavelength metadata written by *i.hyper.import*, either from
the `r3.info` history or from per-band `r.support` fields.

### Point mode

Point mode reads the pixel's whole depth column with a single *r3.what* call
rather than extracting one 2D raster per band, which on a large cube is the
difference between under a second and tens of minutes. It falls back to
per-band extraction if *r3.what* cannot be used.

### Interpretation

A high score means the pixel's spectral shape resembles one of the reference
spectra; it is not by itself an identification. The reference spectra come
from one growing region and one set of instruments, so the absolute score
level is not transferable between scenes. Calibrate a threshold against known
ground for the scene at hand, and prefer **method=consensus** where no such
ground is available, since the agreement count makes disagreement between
methods visible instead of hiding it.

## EXAMPLES

Match all four built-in spectra with the default SAM:

```sh
i.hyper.poppy input=cube output=poppy_max
```

Full consensus run, keeping the per-method maps:

```sh
i.hyper.poppy input=cube output=poppy_consensus method=consensus \
    fusion_mode=stouffer output_prefix=poppy_m
```

Flowering stage only, continuum-removed, over the VNIR alone:

```sh
i.hyper.poppy input=cube output=poppy_flower growth_stage=flowering \
    method=cr_sam -c max_wavelength=1000
```

Inspect one pixel's scores:

```sh
i.hyper.poppy -p input=cube output=dummy coordinates=722821.9,4824854.6
```

## REFERENCES

- Bennington A.L. (2008). *Application of multi-spectral remote sensing for
  crop discrimination in Afghanistan.* PhD thesis, Cranfield University.
- Calderón R., Montes-Borrego M., Landa B.B., Navas-Cortés J.A. &
  Zarco-Tejada P.J. (2014). Detection of downy mildew of opium poppy using
  high-resolution multi-spectral and thermal imagery acquired with an
  unmanned aerial vehicle. *Precision Agriculture*.
- Jia K., Wu B., Tian Y., Li Q. & Du X. (2011). Spectral discrimination of
  opium poppy using field spectrometry. *IEEE Transactions on Geoscience and
  Remote Sensing* 49(9):3414-3422.
- Kruse F.A. et al. (1993). The Spectral Image Processing System (SIPS):
  interactive visualization and analysis of imaging spectrometer data.
  *Remote Sensing of Environment* 44(2-3):145-163.
- Wang J. (2013). *Unsupervised detection of opium poppy fields in Afghanistan
  from EO-1 Hyperion data.* UNB Technical Report No. 286.

## SEE ALSO

*[i.hyper.atcorr](i.hyper.atcorr.html)*,
*[i.hyper.cannabis](i.hyper.cannabis.html)*,
*[i.hyper.import](i.hyper.import.html)*,
*[i.hyper.sleuth](i.hyper.sleuth.html)*,
*[i.hyper.spectroscopy](i.hyper.spectroscopy.html)*,
*[r3.info](https://grass.osgeo.org/grass-stable/manuals/r3.info.html)*,
*[r3.what](https://grass.osgeo.org/grass-stable/manuals/r3.what.html)*

## AUTHORS

Yann Chemin, GRASS Development Team.

Part of the *i.hyper* module family for VNIR-SWIR hyperspectral data
processing in GRASS GIS.
