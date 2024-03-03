@echo off
setlocal enabledelayedexpansion
cls; 1>&2
chcp 65001 1>&2
for /f "usebackq tokens=*" %%i IN (`cd`) DO @set current_dir=%%i
if /i "%current_dir%" == "%systemroot%" (
  set mafft_working_dir="%~dp0"
) else (
  set mafft_working_dir="%current_dir%"
)
pushd "%mafft_working_dir%"

echo; 1>&2
echo Preparing environment to run MAFFT on Windows. 1>&2
echo This may take a while, if real-time scanning by anti-virus software is on. 1>&2

set ROOTDIR=%~d0%~p0
set PATH=/usr/bin/:%PATH%
set MAFFT_BINARIES=/usr/lib/mafft
set TMPDIR=%TMP%
set MAFFT_TMPDIR=%TMPDIR%

REM set TMPDIR=%ROOTDIR%/tmp
REM set MAFFT_TMPDIR=%TMPDIR%
REM If you do not have write permission for Windows temporary folder
REM (typically C:\Users\username\AppData\Local\Temp\), then
REM uncomment (remove REM) the above two lines to use an alternative 
REM temporary folder.

"%ROOTDIR%\usr\bin\bash" "/usr/bin/mafft" %*

popd
exit /b
