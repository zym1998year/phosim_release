@echo off
rem Windows entry point for the compact PhoSim GUI launcher (tools/phosim_gui2.py).
rem Run from a Python environment that has tkinter (e.g. the conda phosim-build env).
python "%~dp0tools\phosim_gui2.py" %*
