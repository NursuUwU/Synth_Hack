$ErrorActionPreference = 'Stop'
$projectRoot = $PSScriptRoot
$candidates = @(
    $env:GENOMIC_PYTHON,
    (Join-Path $projectRoot '.venv\Scripts\python.exe'),
    (Join-Path $env:USERPROFILE '.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe')
)
$pythonPath = $candidates | Where-Object { $_ -and (Test-Path -LiteralPath $_ -PathType Leaf) } | Select-Object -First 1
if (-not $pythonPath) {
    $command = Get-Command python -ErrorAction SilentlyContinue
    if ($command -and $command.Source -notlike '*WindowsApps*') { $pythonPath = $command.Source }
}
if (-not $pythonPath) { throw 'Python 3.11+ required. Set GENOMIC_PYTHON to the executable path.' }
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
    & $pythonPath -m genomic_detective @args
    $resultCode = $LASTEXITCODE
} finally {
    if ($loadedCredential) { Remove-Item Env:\GIGACHAT_CREDENTIALS -ErrorAction SilentlyContinue }
    Pop-Location
}
exit $resultCode
