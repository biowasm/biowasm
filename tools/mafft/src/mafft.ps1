Write-Host;
Write-Host Preparing environment to run MAFFT on Windows.
Write-Host This may take a while, if real-time scanning by anti-virus software is on.

Set-Item Env:Path "/usr/bin;$Env:Path"
Set-Item Env:MAFFT_BINARIES "/usr/lib/mafft"
Set-Item Env:TMPDIR "$Env:TMP"
Set-Item Env:MAFFT_TMPDIR "$Env:TMP"
Set-Item Env:mafft_working_dir "$PWD"

#Set-Item Env:TMPDIR "/tmp"
#Set-Item Env:MAFFT_TMPDIR "/tmp"
# If you do not have write permission for standard temporary folder
# (typically C:\Users\username\AppData\Local\Temp\), then
# uncomment (remove #) the above two lines to use an alternative 
# temporary folder.

#$ROOTDIR=$PSScriptRoot # not supported by powershell versions <= 2
$ROOTDIR=Split-Path -Parent $MyInvocation.MyCommand.Path
$proc = Start-Process -Wait -NoNewWindow -PassThru -FilePath "$ROOTDIR\usr\bin\bash.exe" -ArgumentList "'/usr/bin/mafft' $args"
exit $proc.ExitCode
