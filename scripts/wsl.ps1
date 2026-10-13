param(
    [Parameter(Position=0,Mandatory=$true)]
    [ValidateSet('build','prepare','run','check','plot')][string]$Task,
    [Parameter(ValueFromRemainingArguments=$true)][string[]]$TaskArgs
)
$ProjectRoot = Split-Path $PSScriptRoot -Parent
if ($Task -in @('build','run','check')) {
    & wsl.exe -d Debian -u wang --cd $ProjectRoot -- bash "scripts/$Task.sh" @TaskArgs
} else {
    & wsl.exe -d Debian -u wang --cd $ProjectRoot -- python3 -B "scripts/$Task.py" @TaskArgs
}
exit $LASTEXITCODE
