from pathlib import Path

source = (Path(__file__).resolve().parents[1] / "desmond_gmxmmpbsa_wx.py").read_text()
assert 'GUI_VERSION = "1.0.6"' in source
assert 'ENGINE_VERSION = "9.1.0"' in source
assert 'modules = ["numpy", "scipy", "mdtraj", "networkx", "parmed"]' in source
print("GUI v9 constants and SciPy environment check passed")
