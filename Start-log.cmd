@echo off
rem Diagnostic launcher: same as Start.cmd, but writes all output to start-log.txt.
setlocal
cd /d "%~dp0"
echo Writing diagnostics to start-log.txt ...
echo If the window shows "Local UI", the start worked: close it with Ctrl+C.
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0start.ps1" %* > "%~dp0start-log.txt" 2>&1
echo.
echo Done. Tell Claude: "прочитай start-log.txt"
pause
