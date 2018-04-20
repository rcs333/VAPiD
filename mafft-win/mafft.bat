@echo off
setlocal enabledelayedexpansion

echo; 1>&2
echo Preparing environment to run MAFFT on Windows. 1>&2
echo This may take a while, if real-time scanning by anti-virus software is on. 1>&2

set ROOTDIR=%~d0%~p0
set ROOTDIRWBS=!ROOTDIR:^\\=\\\!
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

"%ROOTDIR%\usr\bin\bash" "%ROOTDIRWBS%\usr\bin\mafft" %*

exit /b
