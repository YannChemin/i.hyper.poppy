"""Tests that libras3d is functional for i.hyper.poppy standalone mode."""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_ras3d_common import (
    WYVERN_PATH, TANAGER_PATH,
    skip_without_ras3d, skip_without_wyvern, skip_without_tanager,
    open_cube_checked, assert_band_valid, install_ras3d_shim, make_wl_sidecar,
)


@skip_without_ras3d
def test_shim_installs():
    """ras3d shim provides grass.script + grass.script.array API."""
    install_ras3d_shim()
    import grass.script as gs
    assert hasattr(gs, 'raster3d_info')
    assert hasattr(gs, 'parser')


@skip_without_ras3d
@skip_without_wyvern
def test_open_wyvern_geotiff():
    """Wyvern 23-band GeoTIFF opens with correct dimensions."""
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    assert r['depths'] == 23
    assert r['rows']   == 7825
    assert r['cols']   == 6003
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_tanager
def test_open_tanager_hdf5():
    """Tanager 426-band HDF5 opens with correct dimensions."""
    import ras3d
    h, r = open_cube_checked(TANAGER_PATH)
    assert r['depths'] == 426
    assert r['rows']   == 732
    assert r['cols']   == 607
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_wyvern
def test_read_all_bands_wyvern():
    """read_all_bands() returns float32 cube [23, 7825, 6003]."""
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    cube = ras3d.read_all_bands(h)
    assert cube.shape == (r['depths'], r['rows'], r['cols'])
    assert cube.dtype == np.float32
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_tanager
def test_read_all_bands_tanager():
    """read_all_bands() on Tanager HDF5 returns valid data."""
    import ras3d
    h, r = open_cube_checked(TANAGER_PATH)
    cube = ras3d.read_all_bands(h)
    assert cube.shape == (r['depths'], r['rows'], r['cols'])
    assert_band_valid(cube[0],   'Tanager band 0')
    assert_band_valid(cube[425], 'Tanager band 425')
    ras3d.close_cube(h)


@skip_without_ras3d
@skip_without_wyvern
def test_extract_band_ras3d(tmp_path):
    """extract_band() in ras3d mode populates band cache and writes GeoTIFF."""
    install_ras3d_shim()
    os.environ['RAS3D_OUTDIR'] = str(tmp_path)
    sys.path.insert(0, '/home/yann/dev/i.hyper.poppy')
    import i_hyper_poppy
    name = i_hyper_poppy.extract_band(WYVERN_PATH, 5)
    from ras3d_grass_shim import get_band_cache
    assert name in get_band_cache()
    assert_band_valid(get_band_cache()[name], 'poppy extract_band 5')


@skip_without_ras3d
@skip_without_wyvern
def test_wavelength_sidecar():
    """get_band_info() reads .wl.json sidecar in ras3d mode."""
    install_ras3d_shim()
    import ras3d
    h, r = open_cube_checked(WYVERN_PATH)
    sidecar, _ = make_wl_sidecar(WYVERN_PATH, r['depths'])
    ras3d.close_cube(h)
    sys.path.insert(0, '/home/yann/dev/i.hyper.poppy')
    import i_hyper_poppy
    bands = i_hyper_poppy.get_band_info(WYVERN_PATH)
    assert len(bands) == r['depths']
    os.unlink(sidecar)


@skip_without_ras3d
@skip_without_wyvern
def test_get_region_from_shim():
    """gs.raster3d_info() shim returns correct metadata."""
    install_ras3d_shim()
    import grass.script as gs
    info = gs.raster3d_info(WYVERN_PATH)
    assert info['depths'] == 23
    assert info['rows']   == 7825
    assert info['cols']   == 6003
