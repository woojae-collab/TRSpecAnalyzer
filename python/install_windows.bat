@echo off
setlocal
cd /d "%~dp0"
echo ============================================================
echo  TRSpecAnalyzer, Python version - installation
echo ============================================================
echo.

set "PY="
py -3 -c "import sys" >nul 2>nul
if not errorlevel 1 set "PY=py -3"
if defined PY goto :havepy
python -c "import sys" >nul 2>nul
if not errorlevel 1 set "PY=python"
if defined PY goto :havepy

echo Python was not found on this computer.
echo.
echo 1. Download Python 3.13, 64-bit, from https://www.python.org/downloads/
echo 2. In the installer, tick "Add python.exe to PATH", then "Install Now".
echo 3. Run this file again.
echo.
pause
exit /b 1

:havepy
echo Python found:
%PY% --version
%PY% -c "import sys; sys.exit(0 if sys.version_info[:2] >= (3, 10) and sys.maxsize > 2**32 else 1)"
if errorlevel 1 goto :oldpy

if exist ".venv\Scripts\python.exe" goto :packages
echo.
echo Creating a private environment in the folder .venv ...
%PY% -m venv .venv
if errorlevel 1 goto :fail

:packages
echo.
echo Installing numpy, scipy, matplotlib and PySide6.
echo This needs an internet connection and takes a few minutes ...
echo.
".venv\Scripts\python.exe" -m pip install --upgrade pip
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto :fail
".venv\Scripts\python.exe" -c "import numpy, scipy, matplotlib, PySide6; print('Packages are in place.')"
if errorlevel 1 goto :fail
echo.
echo ============================================================
echo  Installation finished.
echo  Start the program with run_windows.bat
echo ============================================================
pause
exit /b 0

:oldpy
echo.
echo This Python cannot be used: a 64-bit Python 3.10 or later is needed.
echo Install Python 3.13, 64-bit, from https://www.python.org/downloads/
echo and run this file again.
pause
exit /b 1

:fail
echo.
echo ============================================================
echo  Installation did not finish. Please copy the text above
echo  and send it along.
echo ============================================================
pause
exit /b 1
