# Desmond gmx_MMPBSA wx GUI 1.0.2

This maintenance release fixes a wxWidgets parent/sizer assertion in the Results tab.

## Fixed

The **Parse result** and **Open result file** buttons were created with the notebook page as their parent, but were placed in a sizer owned by a child `result_actions` panel. wxPython 4.2.3 correctly rejected that inconsistent ownership.

The buttons are now created with `result_actions` as their parent:

```python
result_actions = wx.Panel(panel)
self.parse_result_button = wx.Button(result_actions, label="Parse result")
self.open_result_button = wx.Button(result_actions, label="Open result file")
```

The scientific engine remains version 7.1.0 and is unchanged.
