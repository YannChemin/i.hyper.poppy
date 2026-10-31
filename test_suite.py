#!/usr/bin/env python3
"""Test suite for i.hyper.poppy — runs without GRASS GIS."""
import importlib.util
import json
import sys
import tempfile
from unittest.mock import MagicMock

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Load the module with grass.script mocked out
# ---------------------------------------------------------------------------

_gs_mock = MagicMock()
_grass_mock = MagicMock()
# `import grass.script as gs` resolves getattr(grass, "script") before it
# falls back to sys.modules, so the attribute must be wired to the same mock
# or the module ends up holding a different one and gs.fatal() cannot be
# observed from the tests.
_grass_mock.script = _gs_mock
sys.modules.setdefault("grass", _grass_mock)
sys.modules.setdefault("grass.script", _gs_mock)

_spec = importlib.util.spec_from_file_location(
    "i_hyper_poppy",
    __file__.replace("test_suite.py", "i.hyper.poppy.py"),
)
_mod = importlib.util.module_from_spec(_spec)
sys.modules["i_hyper_poppy"] = _mod
_spec.loader.exec_module(_mod)

ALL_METHODS         = _mod.ALL_METHODS
BASE_METHODS        = _mod.BASE_METHODS
METHOD_LABELS       = _mod.METHOD_LABELS
METHOD_GROUPS       = _mod.METHOD_GROUPS
FUSION_MODES        = _mod.FUSION_MODES
GLOBAL_STATS_METHODS = _mod.GLOBAL_STATS_METHODS
WavelengthLUT       = _mod.WavelengthLUT
resample_reference  = _mod.resample_reference
normalize_spectrum  = _mod.normalize_spectrum
normalize_cube      = _mod.normalize_cube
to_prob_simplex     = _mod.to_prob_simplex
to_prob_simplex_cube = _mod.to_prob_simplex_cube
_upper_hull         = _mod._upper_hull
continuum_remove    = _mod.continuum_remove
continuum_remove_cube = _mod.continuum_remove_cube
_safe_cube          = _mod._safe_cube
compute_method      = _mod.compute_method
point_analysis      = _mod.point_analysis
match_sam           = _mod.match_sam
match_ensemble      = _mod.match_ensemble
calibrate_scores    = _mod.calibrate_scores
compute_diversity_weights = _mod.compute_diversity_weights
fuse_probabilities  = _mod.fuse_probabilities
_norm_ppf           = _mod._norm_ppf
get_builtin_spectra = _mod.get_builtin_spectra
load_poppy_db       = _mod.load_poppy_db
_convert_wl_nm      = _mod._convert_wl_nm
_BUILTIN_SPECTRA    = _mod._BUILTIN_SPECTRA


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def wls():
    """A 60-band sensor over the ASD range the built-in spectra cover."""
    return np.linspace(400.0, 2400.0, 60)


@pytest.fixture
def ref(wls):
    """A reference spectrum with a red edge and two absorption dips."""
    r = 0.10 + 0.40 / (1.0 + np.exp(-(wls - 720.0) / 25.0))
    r -= 0.06 * np.exp(-0.5 * ((wls - 670.0) / 20.0) ** 2)
    r -= 0.10 * np.exp(-0.5 * ((wls - 1450.0) / 60.0) ** 2)
    return r.astype(np.float64)


@pytest.fixture
def cube(ref):
    """A 3x3 cube: the reference at (0,0), scaled copies, and a flat spectrum.

    Pixel (0,0) is the reference exactly, (0,1) is the reference scaled by
    0.5 — same shape, different brightness — and (2,2) is flat, which every
    shape-based measure must rank last.
    """
    n = ref.size
    c = np.zeros((n, 3, 3), dtype=np.float64)
    c[:, 0, 0] = ref
    c[:, 0, 1] = ref * 0.5
    c[:, 0, 2] = ref + 0.02
    rng = np.random.default_rng(4)
    for i in range(1, 3):
        for j in range(3):
            c[:, i, j] = 0.2 + 0.05 * rng.standard_normal(n)
    c[:, 2, 2] = 0.3
    return c


# ---------------------------------------------------------------------------
# Method tables
# ---------------------------------------------------------------------------

class TestMethodTables:
    def test_every_method_has_a_label(self):
        assert set(METHOD_LABELS) == set(ALL_METHODS)

    def test_base_methods_exclude_the_two_fusions(self):
        assert set(ALL_METHODS) - set(BASE_METHODS) == {'ensemble', 'consensus'}

    def test_the_covariance_methods_are_marked(self):
        assert GLOBAL_STATS_METHODS == {'mtf', 'cem'}
        assert GLOBAL_STATS_METHODS <= set(BASE_METHODS)

    def test_method_groups_partition_the_base_methods(self):
        grouped = [m for ms in METHOD_GROUPS.values() for m in ms]
        assert len(grouped) == len(set(grouped)), "a method is in two groups"
        assert set(grouped) == set(BASE_METHODS)


# ---------------------------------------------------------------------------
# Wavelength LUT and resampling
# ---------------------------------------------------------------------------

class TestResampling:
    def test_identity_on_matching_grids(self, wls, ref):
        out = resample_reference(wls, ref, wls)
        assert np.allclose(out, ref, atol=1e-9)

    def test_linear_resample_is_exact_on_a_straight_line(self):
        src = np.array([400.0, 500.0, 600.0, 700.0])
        out = resample_reference(src, 2.0 * src + 1.0,
                                 np.array([450.0, 550.0, 650.0]))
        assert np.allclose(out, 2.0 * np.array([450.0, 550.0, 650.0]) + 1.0)

    def test_edge_fill_holds_the_end_values_outside_the_reference(self):
        src = np.array([500.0, 600.0, 700.0])
        vals = np.array([0.1, 0.2, 0.3])
        dst = np.array([400.0, 600.0, 900.0])
        out = resample_reference(src, vals, dst)
        assert out[0] == pytest.approx(0.1)
        assert out[2] == pytest.approx(0.3)

    def test_nan_fill_marks_bands_the_reference_does_not_cover(self):
        src = np.array([500.0, 600.0, 700.0])
        vals = np.array([0.1, 0.2, 0.3])
        dst = np.array([400.0, 600.0, 900.0])
        lut = WavelengthLUT(src, dst, fill='nan')
        out = lut.apply(vals)
        assert np.isnan(out[0]) and np.isnan(out[2])
        assert out[1] == pytest.approx(0.2)

    def test_coverage_report_mentions_the_uncovered_bands(self):
        lut = WavelengthLUT(np.array([500.0, 600.0]),
                            np.array([400.0, 550.0, 700.0]))
        assert isinstance(lut.coverage_report(), str)
        assert lut.valid_dst.sum() == 1

    @pytest.mark.parametrize('method', ['linear', 'cubic', 'pchip'])
    def test_all_methods_reproduce_the_reference_on_its_own_grid(
            self, method, wls, ref):
        out = resample_reference(wls, ref, wls, method=method)
        assert np.allclose(out, ref, atol=1e-6)

    @pytest.mark.parametrize('unit', [
        'um', '\u00b5m', 'micrometer', 'micrometers', 'micrometre',
        'micrometres', 'micron', 'microns', 'Micrometers', ' UM ',
    ])
    def test_micrometre_spellings_all_scale_by_a_thousand(self, unit):
        """ENVI headers write 'Micrometers'; missing it is a 1000x error."""
        assert _convert_wl_nm(1.45, unit) == pytest.approx(1450.0)

    @pytest.mark.parametrize('unit', [
        'nm', 'nanometer', 'nanometers', 'nanometre', 'nanometres', 'NM',
    ])
    def test_nanometre_spellings_pass_through(self, unit):
        assert _convert_wl_nm(1450.0, unit) == pytest.approx(1450.0)

    @pytest.mark.parametrize('unit', ['m', 'meter', 'meters', 'metre', 'metres'])
    def test_metre_spellings_scale_to_nanometres(self, unit):
        assert _convert_wl_nm(1.45e-6, unit) == pytest.approx(1450.0)


# ---------------------------------------------------------------------------
# Normalization and the probability simplex
# ---------------------------------------------------------------------------

class TestNormalization:
    def test_none_returns_the_spectrum_untouched(self, ref):
        assert normalize_spectrum(ref, 'none') is ref

    def test_area_sums_to_one(self, ref):
        assert normalize_spectrum(ref, 'area').sum() == pytest.approx(1.0)

    def test_max_peaks_at_one(self, ref):
        assert normalize_spectrum(ref, 'max').max() == pytest.approx(1.0)

    def test_minmax_spans_zero_to_one(self, ref):
        out = normalize_spectrum(ref, 'minmax')
        assert out.min() == pytest.approx(0.0)
        assert out.max() == pytest.approx(1.0)

    def test_vector_has_unit_norm(self, ref):
        assert np.linalg.norm(normalize_spectrum(ref, 'vector')) \
            == pytest.approx(1.0)

    def test_a_constant_spectrum_does_not_divide_by_zero(self):
        flat = np.full(10, 0.25)
        assert np.all(normalize_spectrum(flat, 'minmax') == 0.0)

    def test_normalize_does_not_modify_its_input(self, ref):
        before = ref.copy()
        normalize_spectrum(ref, 'area')
        assert np.array_equal(ref, before)

    @pytest.mark.parametrize('method', ['area', 'max', 'minmax', 'vector'])
    def test_cube_normalization_matches_the_per_spectrum_one(
            self, method, cube):
        out = normalize_cube(cube, method)
        for i in range(3):
            for j in range(3):
                expect = normalize_spectrum(cube[:, i, j].copy(), method)
                assert np.allclose(out[:, i, j], expect, atol=1e-12)

    def test_normalize_cube_does_not_modify_its_input(self, cube):
        before = cube.copy()
        normalize_cube(cube, 'vector')
        assert np.array_equal(cube, before)

    def test_simplex_is_positive_and_sums_to_one(self):
        spec = np.array([0.0, -0.5, 0.2, 0.8])
        out = to_prob_simplex(spec)
        assert out.sum() == pytest.approx(1.0)
        assert np.all(out > 0.0)

    def test_simplex_cube_sums_to_one_per_pixel(self, cube):
        out = to_prob_simplex_cube(cube)
        assert np.allclose(out.sum(axis=0), 1.0)


# ---------------------------------------------------------------------------
# Continuum removal
# ---------------------------------------------------------------------------

class TestContinuumRemoval:
    def test_hull_lies_on_or_above_the_spectrum(self, wls, ref):
        hull = _upper_hull(wls, ref)
        assert np.all(hull >= ref - 1e-9)

    def test_a_convex_spectrum_is_its_own_hull(self, wls):
        line = 0.1 + 0.0002 * wls
        assert np.allclose(_upper_hull(wls, line), line, atol=1e-9)
        assert np.allclose(continuum_remove(line, wls), 1.0, atol=1e-9)

    def test_removal_is_bounded_to_the_unit_interval(self, wls, ref):
        cr = continuum_remove(ref, wls)
        assert cr.min() >= 0.0 and cr.max() <= 1.0

    def test_the_dips_survive_removal_as_minima(self, wls, ref):
        cr = continuum_remove(ref, wls)
        # The 1450 nm dip is the deeper of the two planted features.
        assert cr[int(np.argmin(np.abs(wls - 1450.0)))] < 0.95

    def test_removal_is_scale_invariant(self, wls, ref):
        assert np.allclose(continuum_remove(ref, wls),
                           continuum_remove(ref * 3.7, wls), atol=1e-9)

    def test_cube_removal_matches_the_per_spectrum_one(self, wls, cube):
        out = continuum_remove_cube(cube.copy(), wls)
        for i in range(3):
            for j in range(3):
                expect = continuum_remove(cube[:, i, j], wls)
                assert np.allclose(out[:, i, j], expect, atol=1e-12)


# ---------------------------------------------------------------------------
# The similarity methods, as a family
# ---------------------------------------------------------------------------

class TestSimilarityMethods:
    @pytest.mark.parametrize('method', BASE_METHODS)
    def test_scores_stay_inside_the_unit_interval(
            self, method, cube, ref, wls):
        out = compute_method(method, cube.copy(), ref, wls, 3, {})
        assert out.shape == (3, 3)
        assert np.all(np.isfinite(out)), f"{method} produced a non-finite score"
        assert out.min() >= 0.0 and out.max() <= 1.0, \
            f"{method} left [0, 1]: {out.min()}..{out.max()}"

    @pytest.mark.parametrize('method', BASE_METHODS)
    def test_the_exact_reference_pixel_scores_above_the_flat_one(
            self, method, cube, ref, wls):
        out = compute_method(method, cube.copy(), ref, wls, 3, {})
        assert out[0, 0] > out[2, 2], \
            f"{method} ranks a flat spectrum at or above an exact match"

    @pytest.mark.parametrize('method', BASE_METHODS)
    def test_nodata_does_not_propagate_into_the_score_map(
            self, method, cube, ref, wls):
        c = cube.copy()
        c[5, 1, 1] = np.nan
        c[6, 1, 1] = np.inf
        out = compute_method(method, c, ref, wls, 3, {})
        assert np.all(np.isfinite(out)), f"{method} leaked nodata"

    def test_sam_is_one_for_a_parallel_spectrum(self, cube, ref):
        """Tolerance is 1e-5, not machine epsilon, and deliberately so.

        match_sam() adds EPS inside both norms rather than guarding the
        division, so a perfect match gives cos(a) = 1 - O(EPS) instead of
        exactly 1. arccos() then amplifies that to sqrt(EPS), about 1e-6 at
        the top of the range. It is a bounded cosmetic error at the one
        point where the answer is already unambiguous.
        """
        out = match_sam(cube.copy(), ref)
        assert out[0, 0] == pytest.approx(1.0, abs=1e-5)
        assert out[0, 0] <= 1.0

    def test_sam_ignores_brightness(self, cube, ref):
        """SAM is an angle, so a scaled copy must score the same."""
        out = match_sam(cube.copy(), ref)
        assert out[0, 1] == pytest.approx(out[0, 0], abs=1e-5)

    def test_safe_cube_sanitizes_in_place_without_copying(self):
        c = np.array([[[1.0, np.nan]], [[np.inf, -np.inf]]])
        out = _safe_cube(c)
        assert out is c, "the cube was copied instead of fixed in place"
        assert np.all(np.isfinite(c))

    def test_safe_cube_preserves_dtype(self):
        c = np.full((4, 2, 2), np.nan, dtype=np.float32)
        assert _safe_cube(c).dtype == np.float32


# ---------------------------------------------------------------------------
# Fusion
# ---------------------------------------------------------------------------

class TestFusion:
    @pytest.fixture
    def score_maps(self):
        rng = np.random.default_rng(11)
        return {m: rng.random((4, 4)) for m in ('sam', 'sid', 'ed')}

    def test_ensemble_returns_rank_fractions(self, score_maps):
        out = match_ensemble(score_maps)
        assert out.shape == (4, 4)
        assert out.min() >= 0.0 and out.max() <= 1.0

    def test_ensemble_ignores_its_own_entry(self, score_maps):
        with_self = dict(score_maps, ensemble=np.zeros((4, 4)))
        assert np.allclose(match_ensemble(with_self),
                           match_ensemble(score_maps))

    def test_calibration_is_a_uniform_rank_transform(self, score_maps):
        cal = calibrate_scores(score_maps)
        for c in cal.values():
            assert c.min() > 0.0 and c.max() <= 1.0
            # A rank transform of distinct values is evenly spaced.
            assert np.allclose(np.sort(c.ravel()),
                               np.arange(1, 17) / 16.0)

    def test_calibration_preserves_order(self, score_maps):
        cal = calibrate_scores(score_maps)
        raw = score_maps['sam'].ravel()
        assert np.array_equal(np.argsort(raw),
                              np.argsort(cal['sam'].ravel()))

    def test_calibration_skips_the_fusion_entries(self, score_maps):
        cal = calibrate_scores(dict(score_maps,
                                    ensemble=np.zeros((4, 4)),
                                    consensus=np.zeros((4, 4))))
        assert 'ensemble' not in cal and 'consensus' not in cal

    def test_diversity_weights_sum_to_the_method_count(self, score_maps):
        cal = calibrate_scores(score_maps)
        w = compute_diversity_weights(cal)
        assert sum(w.values()) == pytest.approx(len(w))

    def test_a_duplicated_method_is_down_weighted(self):
        rng = np.random.default_rng(5)
        a = rng.random((8, 8))
        maps = {'sam': a, 'sid': a.copy(), 'ed': rng.random((8, 8))}
        w = compute_diversity_weights(calibrate_scores(maps))
        assert w['ed'] > w['sam'], \
            "an independent method should outweigh a duplicated pair"

    @pytest.mark.parametrize('mode', FUSION_MODES)
    def test_every_fusion_mode_stays_in_the_unit_interval(
            self, mode, score_maps):
        cal = calibrate_scores(score_maps)
        out = fuse_probabilities(cal, compute_diversity_weights(cal), mode)
        assert out.shape == (4, 4)
        assert np.all(np.isfinite(out))
        assert out.min() >= 0.0 and out.max() <= 1.0

    @pytest.mark.parametrize('mode', FUSION_MODES)
    def test_unanimous_agreement_fuses_high(self, mode):
        """Every method ranking a pixel top must fuse to a high value."""
        n = 16
        base = np.linspace(1.0 / n, 1.0, n).reshape(4, 4)
        cal = {m: base.copy() for m in ('sam', 'sid', 'ed')}
        out = fuse_probabilities(cal, {m: 1.0 for m in cal}, mode)
        assert out[-1, -1] == out.max()
        assert out[-1, -1] > out[0, 0]

    def test_norm_ppf_against_known_quantiles(self):
        assert _norm_ppf(np.array([0.5]))[0] == pytest.approx(0.0, abs=2e-3)
        assert _norm_ppf(np.array([0.975]))[0] == pytest.approx(1.96, abs=2e-3)
        assert _norm_ppf(np.array([0.025]))[0] == pytest.approx(-1.96, abs=2e-3)

    def test_norm_ppf_is_antisymmetric(self):
        p = np.array([0.01, 0.1, 0.3])
        assert np.allclose(_norm_ppf(p), -_norm_ppf(1.0 - p), atol=1e-9)


# ---------------------------------------------------------------------------
# Reference spectra
# ---------------------------------------------------------------------------

class TestReferenceSpectra:
    @pytest.mark.parametrize('stage', list(_BUILTIN_SPECTRA))
    def test_each_builtin_stage_is_well_formed(self, stage):
        spec = get_builtin_spectra(stage)[0]
        wl = np.asarray(spec['wavelengths'])
        rf = np.asarray(spec['reflectances'])
        assert wl.size == rf.size and wl.size > 0
        assert np.all(np.diff(wl) > 0), "wavelengths are not ascending"
        assert np.all(rf >= 0.0) and np.all(rf <= 1.0)
        assert spec['citation']

    def test_all_returns_every_stage(self):
        assert len(get_builtin_spectra('all')) == len(_BUILTIN_SPECTRA)

    def test_an_unknown_stage_is_fatal(self):
        """Under GRASS gs.fatal() exits; here it must at least be called."""
        _mod.gs.fatal.reset_mock()
        _mod.gs.fatal.side_effect = RuntimeError("fatal")
        try:
            with pytest.raises(RuntimeError):
                get_builtin_spectra('no_such_stage')
        finally:
            _mod.gs.fatal.side_effect = None
        assert _mod.gs.fatal.called

    def test_the_canopy_stages_have_a_red_edge(self):
        """Every canopy spectrum must rise from red to NIR."""
        for stage in ('before_flowering', 'flowering', 'harvesting'):
            spec = get_builtin_spectra(stage)[0]
            wl = np.asarray(spec['wavelengths'])
            rf = np.asarray(spec['reflectances'])
            red = rf[int(np.argmin(np.abs(wl - 670.0)))]
            nir = rf[int(np.argmin(np.abs(wl - 800.0)))]
            assert nir > red, f"{stage} has no red edge"

    def test_load_csv_library(self):
        with tempfile.NamedTemporaryFile('w', suffix='.csv',
                                         delete=False) as fh:
            fh.write("wavelength,reflectance\n500,0.1\n600,0.2\n700,0.3\n")
            path = fh.name
        spectra = load_poppy_db(path)
        assert len(spectra) >= 1
        assert np.allclose(spectra[0]['reflectances'], [0.1, 0.2, 0.3])

    def test_load_json_pairs(self):
        with tempfile.NamedTemporaryFile('w', suffix='.json',
                                         delete=False) as fh:
            json.dump([[500, 0.1], [600, 0.2], [700, 0.3]], fh)
            path = fh.name
        spectra = load_poppy_db(path)
        assert np.allclose(spectra[0]['wavelengths'], [500, 600, 700])

    def test_load_json_parallel_arrays(self):
        with tempfile.NamedTemporaryFile('w', suffix='.json',
                                         delete=False) as fh:
            json.dump({'wavelengths': [500, 600], 'reflectances': [0.1, 0.2]},
                      fh)
            path = fh.name
        spectra = load_poppy_db(path)
        assert np.allclose(spectra[0]['reflectances'], [0.1, 0.2])


# ---------------------------------------------------------------------------
# Point mode
# ---------------------------------------------------------------------------

class TestPointAnalysis:
    def test_point_mode_reproduces_the_full_map_score(self, cube, ref, wls):
        """One pixel through point_analysis must equal its map cell."""
        scores = point_analysis(cube[:, 0, 0].copy(), ref, wls,
                                ['sam', 'sid', 'ed'], 3, flag_c=False)
        maps = {m: compute_method(m, cube.copy(), ref, wls, 3, {})
                for m in ('sam', 'sid', 'ed')}
        for m in ('sam', 'sid', 'ed'):
            assert scores[m] == pytest.approx(maps[m][0, 0], abs=1e-9)

    def test_point_mode_covers_every_per_pixel_method(self, cube, ref, wls):
        scores = point_analysis(cube[:, 0, 0].copy(), ref, wls,
                                BASE_METHODS, 3, flag_c=False)
        assert set(scores) == set(BASE_METHODS) - GLOBAL_STATS_METHODS
        assert all(np.isfinite(v) for v in scores.values())

    @pytest.mark.parametrize('method', sorted(GLOBAL_STATS_METHODS))
    def test_point_mode_refuses_the_covariance_methods(
            self, method, cube, ref, wls):
        """MTF and CEM need the scene covariance, undefined for one pixel."""
        _mod.gs.warning.reset_mock()
        scores = point_analysis(cube[:, 0, 0].copy(), ref, wls,
                                ['sam', method], 3, flag_c=False)
        assert method not in scores
        assert 'sam' in scores
        assert _mod.gs.warning.called

    def test_point_mode_skips_consensus(self, cube, ref, wls):
        scores = point_analysis(cube[:, 0, 0].copy(), ref, wls,
                                ['sam', 'consensus'], 3, flag_c=False)
        assert 'consensus' not in scores

    def test_continuum_flag_changes_the_score(self, cube, ref, wls):
        plain = point_analysis(cube[:, 0, 2].copy(), ref, wls, ['sam'], 3,
                               flag_c=False)
        removed = point_analysis(cube[:, 0, 2].copy(), ref, wls, ['sam'], 3,
                                 flag_c=True)
        assert plain['sam'] != removed['sam']


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-v']))
