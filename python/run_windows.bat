@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" goto :noinst
echo TRSpecAnalyzer, Python version. This window shows messages of the
echo program; it closes when the program is closed.
".venv\Scripts\python.exe" -m trspecgui %*
if errorlevel 1 goto :err
exit /b 0

:noinst
echo The program is not installed yet. Run install_windows.bat first.
pause
exit /b 1

:err
echo.
echo The program ended with an error. Please copy the text above and
echo send it along.
pause
exit /b 1
