# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_all

def _try_collect(pkg):
    """collect_all fails if the package isn't installed — return empty lists then."""
    try:
        return collect_all(pkg)
    except Exception:
        return [], [], []

# Mandatory packages
datas_bio, binaries_bio, hiddens_bio = collect_all('Bio')
datas_pd,  binaries_pd,  hiddens_pd  = collect_all('pandas')

# Optional packages (wrapped so the spec doesn't break if they're absent)
datas_mpl, binaries_mpl, hiddens_mpl = _try_collect('matplotlib')
datas_pg,  binaries_pg,  hiddens_pg  = _try_collect('pygame')

a = Analysis(
    ['gui/matrix_app.py'],
    pathex=[
        'scripts',   # makes print_matrix + fetch_group_order importable
        'src',       # makes xlms importable
    ],
    binaries=binaries_bio + binaries_pd + binaries_mpl + binaries_pg,
    datas=datas_bio + datas_pd + datas_mpl + datas_pg,
    hiddenimports=[
        # Project-local modules found only via sys.path at runtime
        'print_matrix',
        'fetch_group_order',
        'xlms',
        'xlms.io',
        'xlms.models',
        # Deferred / conditional imports inside print_matrix.py
        'scipy',
        'scipy.cluster',
        'scipy.cluster.hierarchy',
        'scipy.spatial',
        'scipy.spatial.distance',
        # matplotlib — listed explicitly because matrix_app.py wraps the import
        # in try/except, which PyInstaller's static analyser skips
        'matplotlib',
        'matplotlib.cm',
        'matplotlib.colors',
        'matplotlib.pyplot',
        'matplotlib.backends.backend_agg',
    ] + hiddens_bio + hiddens_pd + hiddens_mpl + hiddens_pg,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['tkinter.test', 'test', 'unittest'],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='XL-MS Matrix Viewer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,      # no terminal window — GUI app
    onefile=True,       # single .exe
    icon=None,          # set to an .ico path if you have one
)
