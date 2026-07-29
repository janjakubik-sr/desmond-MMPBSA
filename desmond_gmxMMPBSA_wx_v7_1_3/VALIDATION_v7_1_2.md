# Validation of wx GUI 1.0.2

- Python source compiled successfully with `py_compile`.
- All direct `SetFont(wx.FontInfo(...))` calls remain corrected from GUI 1.0.1.
- The Results-tab action buttons now have the same parent as the sizer-containing window.
- A static ownership audit found no other nested-panel/sizer parent mismatch in the GUI construction methods.
- The scientific engine and molecular-modelling files are byte-identical to release 7.1.1.

Interactive rendering could not be executed in the build container because wxPython and a graphical display are not installed there.
