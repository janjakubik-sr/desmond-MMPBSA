from pathlib import Path

SOURCE = Path(__file__).resolve().parents[1] / "desmond_gmxmmpbsa_wx.py"
text = SOURCE.read_text(encoding="utf-8")

assert 'result_actions = wx.Panel(panel)' in text
assert 'self.parse_result_button = wx.Button(result_actions, label="Parse result")' in text
assert 'self.open_result_button = wx.Button(result_actions, label="Open result file")' in text
assert 'self.parse_result_button = wx.Button(panel' not in text
assert 'self.open_result_button = wx.Button(panel' not in text

# Audit the other nested panel/sizer blocks as well.
assert 'self.ligand_resname = wx.TextCtrl(ligand_row' in text
assert 'self.detect_ligand_button = wx.Button(ligand_row' in text
assert 'self.ligand_charge = wx.TextCtrl(parameter_row' in text
assert 'self.stride = wx.SpinCtrl(parameter_row' in text
assert 'self.solvent_model = wx.Choice(parameter_row' in text
assert 'self.mpi_processes = wx.SpinCtrl(' in text

print("wx parent/sizer regression test passed")
