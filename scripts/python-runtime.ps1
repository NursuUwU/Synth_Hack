# Shared by run.ps1 and the double-click launcher. No user-specific runtime paths.
function Test-GenomicPython {
    param([string]$Executable)
    if (-not $Executable -or -not (Test-Path -LiteralPath $Executable -PathType Leaf)) { return $false }
    try {
        $null = & $Executable -I -c "import sys, sqlite3, ssl; sys.exit(0 if sys.version_info >= (3,11) else 1)" 2>$null
        return $LASTEXITCODE -eq 0
    } catch { return $false }
}

function Install-GenomicPython {
    param([string]$ProjectRoot)
    $settings = Get-Content -LiteralPath (Join-Path $ProjectRoot 'config/python-runtime.json') -Raw | ConvertFrom-Json
    $architecture = if ($env:PROCESSOR_ARCHITEW6432) { $env:PROCESSOR_ARCHITEW6432 } else { $env:PROCESSOR_ARCHITECTURE }
    $kind = switch ($architecture.ToUpperInvariant()) {
        'AMD64' { 'amd64' }
        'ARM64' { 'arm64' }
        'X86' { 'win32' }
        default { throw "Unsupported Windows architecture: $architecture. Install Python 3.11+ manually." }
    }
    $download = $settings.windows.$kind
    $runtimeRoot = [IO.Path]::GetFullPath((Join-Path $ProjectRoot '.runtime'))
    $target = [IO.Path]::GetFullPath((Join-Path $runtimeRoot "python-$($settings.version)-$kind"))
    $executable = Join-Path $target 'python.exe'
    if (Test-GenomicPython $executable) { return $executable }
    if (Test-Path -LiteralPath $target) { throw "Incomplete portable Python in $target. Rename that folder and start again." }
    New-Item -ItemType Directory -Path $runtimeRoot -Force | Out-Null
    $identifier = [Guid]::NewGuid().ToString('N')
    $bundled = Join-Path $ProjectRoot ('resources/' + [IO.Path]::GetFileName(([Uri]$download.url).AbsolutePath))
    $useBundled = Test-Path -LiteralPath $bundled -PathType Leaf
    $archive = if ($useBundled) { $bundled } else { Join-Path $runtimeRoot "python-$identifier.zip" }
    $staging = [IO.Path]::GetFullPath((Join-Path $runtimeRoot "install-$identifier"))
    if ($useBundled) {
        Write-Host "Preparing included portable Python $($settings.version)..."
    } else {
        Write-Host "Downloading official portable Python $($settings.version) from python.org..."
        [Net.ServicePointManager]::SecurityProtocol = [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12
        Invoke-WebRequest -UseBasicParsing -Uri $download.url -OutFile $archive -TimeoutSec 180
    }
    if ((Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLowerInvariant() -ne $download.sha256) {
        throw 'Python download checksum mismatch. The downloaded program will not be run.'
    }
    Expand-Archive -LiteralPath $archive -DestinationPath $staging
    # Verify both absolute paths remain inside the project runtime directory before moving.
    $prefix = $runtimeRoot.TrimEnd([IO.Path]::DirectorySeparatorChar) + [IO.Path]::DirectorySeparatorChar
    if (-not $target.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase) -or
        -not $staging.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) { throw 'Unsafe runtime destination' }
    Move-Item -LiteralPath $staging -Destination $target
    if (-not $useBundled) { Remove-Item -LiteralPath $archive }
    if (-not (Test-GenomicPython $executable)) { throw 'Portable Python could not start. Install Python 3.11+ from python.org.' }
    return $executable
}

function Resolve-GenomicPython {
    param([string]$ProjectRoot, [switch]$AllowDownload)
    if ($env:GENOMIC_PYTHON) {
        if (Test-GenomicPython $env:GENOMIC_PYTHON) { return $env:GENOMIC_PYTHON }
        throw 'GENOMIC_PYTHON must point to a working Python 3.11+ with sqlite3 and ssl.'
    }
    $candidates = @((Join-Path $ProjectRoot '.venv/Scripts/python.exe'))
    $runtimeRoot = Join-Path $ProjectRoot '.runtime'
    if (Test-Path -LiteralPath $runtimeRoot) {
        $candidates += @(Get-ChildItem -LiteralPath $runtimeRoot -Filter 'python-*' -Directory |
            ForEach-Object { Join-Path $_.FullName 'python.exe' })
    }
    foreach ($name in @('python', 'python3')) {
        $command = Get-Command $name -ErrorAction SilentlyContinue
        if ($command -and $command.Source -notlike '*WindowsApps*') { $candidates += $command.Source }
    }
    $launcher = Get-Command py -ErrorAction SilentlyContinue
    if ($launcher) {
        try {
            $located = & $launcher.Source -3 -c 'import sys; print(sys.executable)' 2>$null
            if ($LASTEXITCODE -eq 0 -and $located) { $candidates += [string]$located }
        } catch { }
    }
    if ($env:LOCALAPPDATA) {
        $installed = Join-Path $env:LOCALAPPDATA 'Programs/Python'
        if (Test-Path -LiteralPath $installed) {
            $candidates += @(Get-ChildItem -LiteralPath $installed -Directory |
                ForEach-Object { Join-Path $_.FullName 'python.exe' })
        }
    }
    foreach ($candidate in $candidates) {
        if (Test-GenomicPython $candidate) { return $candidate }
    }
    if ($AllowDownload) { return Install-GenomicPython -ProjectRoot $ProjectRoot }
    throw 'Python 3.11+ not found. Double-click Start.cmd to download a project-local Python, or set GENOMIC_PYTHON.'
}
