@echo off
rem Windows entry point for PhoSim (equivalent of the POSIX `phosim` shell wrapper).
rem %~dp0 resolves this script's directory so it can be invoked from anywhere.
python "%~dp0phosim.py" %*
