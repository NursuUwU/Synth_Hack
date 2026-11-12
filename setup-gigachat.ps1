$ErrorActionPreference = 'Stop'
$secretDirectory = Join-Path $PSScriptRoot '.secrets'
$secretPath = Join-Path $secretDirectory 'gigachat.dpapi'
Write-Host 'Paste the Authorization Key from your GigaChat project. Input is hidden.'
$secureKey = Read-Host 'GigaChat Authorization Key' -AsSecureString
if ($secureKey.Length -eq 0) { throw 'The key is empty; nothing was saved.' }
New-Item -ItemType Directory -Path $secretDirectory -Force | Out-Null
# On Windows, without -Key, ConvertFrom-SecureString uses the current user's DPAPI.
$encrypted = ConvertFrom-SecureString -SecureString $secureKey
Set-Content -LiteralPath $secretPath -Value $encrypted -Encoding utf8 -NoNewline
$secureKey.Dispose()
Write-Host 'Saved encrypted credentials for this Windows user. The key was not displayed.'
Write-Host 'Next: .\run.ps1 llm-check --synthetic'
