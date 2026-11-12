param([switch]$External)
$ErrorActionPreference = 'Stop'
. (Join-Path $PSScriptRoot 'scripts/python-runtime.ps1')
$validationPython = Resolve-GenomicPython -ProjectRoot $PSScriptRoot
$env:PYTHONUTF8 = '1'
$env:PYTHONUNBUFFERED = '1'
$validationScript = if ($External) { 'scripts/validate_external.py' } else { 'scripts/validate_controls.py' }
& $validationPython (Join-Path $PSScriptRoot $validationScript) @args
exit $LASTEXITCODE
