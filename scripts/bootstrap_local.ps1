param(
  [string]$PythonExe = "python",
  [string]$VenvDir = ".venv"
)

$ErrorActionPreference = "Stop"

& $PythonExe -m venv $VenvDir
& "$VenvDir\\Scripts\\python.exe" -m pip install --upgrade pip
& "$VenvDir\\Scripts\\pip.exe" install -e .

New-Item -ItemType Directory -Force projects | Out-Null
New-Item -ItemType Directory -Force logs | Out-Null

Write-Host ""
Write-Host "Bootstrap complete."
Write-Host "Activate with:"
Write-Host "  $VenvDir\\Scripts\\Activate.ps1"
Write-Host ""
Write-Host "Then you can run:"
Write-Host "  binderdesign-cuilab --help"
Write-Host "  binderdesign-cuilab run-project --help"
Write-Host "  binderdesign-cuilab table --help"
