$ErrorActionPreference = 'Stop'
$projectRoot = $PSScriptRoot
. (Join-Path $projectRoot 'scripts/python-runtime.ps1')
$pythonPath = Resolve-GenomicPython -ProjectRoot $projectRoot
$env:PYTHONUTF8 = '1'
$loadedCredential = $false
$secretPath = Join-Path $projectRoot '.secrets\gigachat.dpapi'
if ((('--llm' -in $args) -or ('llm-check' -in $args)) -and -not $env:GIGACHAT_CREDENTIALS -and (Test-Path -LiteralPath $secretPath)) {
    $secureKey = Get-Content -LiteralPath $secretPath -Raw | ConvertTo-SecureString
    try {
        $env:GIGACHAT_CREDENTIALS = [System.Net.NetworkCredential]::new('', $secureKey).Password
        $loadedCredential = $true
    } finally {
        $secureKey.Dispose()
    }
}
Push-Location -LiteralPath $projectRoot
try {
    # Explicit project path also supports python.org's isolated embeddable distribution.
    & $pythonPath -I -c "import sys,runpy; sys.path.insert(0,sys.argv.pop(1)); runpy.run_module('genomic_detective',run_name='__main__')" $projectRoot @args
    $resultCode = $LASTEXITCODE
} finally {
    if ($loadedCredential) { Remove-Item Env:\GIGACHAT_CREDENTIALS -ErrorAction SilentlyContinue }
    Pop-Location
}
exit $resultCode
