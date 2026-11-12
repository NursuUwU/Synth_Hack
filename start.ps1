$ErrorActionPreference = 'Stop'
try {
    . (Join-Path $PSScriptRoot 'scripts/python-runtime.ps1')
    $env:GENOMIC_PYTHON = Resolve-GenomicPython -ProjectRoot $PSScriptRoot -AllowDownload
    & (Join-Path $PSScriptRoot 'run.ps1') serve --open @args
} catch {
    Write-Host "START FAILED: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
}
