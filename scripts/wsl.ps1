[CmdletBinding(PositionalBinding=$false)]
param(
    [Parameter(Position=0,Mandatory=$true)]
    [ValidateSet('build','prepare','run','check','plot','geometry','study','viewer')][string]$Task,
    [string]$Distribution = $env:POREBRINE_WSL_DISTRO,
    [string]$User = $env:POREBRINE_WSL_USER,
    [Parameter(ValueFromRemainingArguments=$true)][string[]]$TaskArgs
)
$ProjectRoot = Split-Path $PSScriptRoot -Parent
$WslOptions = @()
if ($Distribution) { $WslOptions += @('--distribution', $Distribution) }
if ($User) { $WslOptions += @('--user', $User) }
if ($Task -in @('build','run','check','geometry','study')) {
    $Script = if ($Task -eq 'geometry') { 'build_geometry' } else { $Task }
    & wsl.exe @WslOptions --cd $ProjectRoot -- bash "scripts/$Script.sh" @TaskArgs
} elseif ($Task -eq 'viewer') {
    & wsl.exe @WslOptions --cd $ProjectRoot -- python3 -B scripts/experiments/export_reactive_viewer.py @TaskArgs
} else {
    & wsl.exe @WslOptions --cd $ProjectRoot -- python3 -B "scripts/$Task.py" @TaskArgs
}
exit $LASTEXITCODE
