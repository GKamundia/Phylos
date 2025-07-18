# Full Auspice Build Script
# This requires Node.js and creates a complete self-contained dashboard

Write-Host "🚀 Building complete Auspice dashboard..." -ForegroundColor Green

# Check if Node.js is installed
try {
    $nodeVersion = node --version
    Write-Host "✅ Node.js is installed: $nodeVersion" -ForegroundColor Green
} catch {
    Write-Host "❌ Node.js is not installed. Please install Node.js first." -ForegroundColor Red
    Write-Host "Download from: https://nodejs.org/" -ForegroundColor Yellow
    exit 1
}

# Check if auspice is installed
try {
    $auspiceVersion = auspice --version
    Write-Host "✅ Auspice is installed: $auspiceVersion" -ForegroundColor Green
} catch {
    Write-Host "📦 Installing Auspice..." -ForegroundColor Yellow
    npm install -g auspice
}

# Check if JSON files exist
if (-not (Test-Path "auspice") -or -not (Get-ChildItem "auspice\*.json" -ErrorAction SilentlyContinue)) {
    Write-Host "❌ No Auspice JSON files found in .\auspice\" -ForegroundColor Red
    exit 1
}

Write-Host "✅ Found Auspice JSON files:" -ForegroundColor Green
Get-ChildItem "auspice\*.json" | Select-Object Name, Length

# Validate JSON files
Write-Host "🔍 Validating JSON files..." -ForegroundColor Cyan
Get-ChildItem "auspice\*.json" | ForEach-Object {
    try {
        Get-Content $_.FullName -Raw | ConvertFrom-Json | Out-Null
        Write-Host "✅ $($_.Name) is valid" -ForegroundColor Green
    } catch {
        Write-Host "❌ $($_.Name) is invalid JSON" -ForegroundColor Red
        exit 1
    }
}

# Build static site
Write-Host "🏗️  Building static Auspice site..." -ForegroundColor Cyan
if (Test-Path "build") {
    Remove-Item "build" -Recurse -Force
}
New-Item -ItemType Directory -Path "build" -Force | Out-Null

# Use auspice build command
auspice build --extend .\auspice --outputDir .\build

Write-Host "✅ Build complete! Static site generated in .\build\" -ForegroundColor Green

# List generated files
Write-Host "📂 Generated files:" -ForegroundColor Cyan
Get-ChildItem "build" -Recurse | Select-Object Name, Length

# Test locally option
$choice = Read-Host "🌐 Would you like to test the dashboard locally? (y/n)"
if ($choice -eq "y" -or $choice -eq "Y") {
    Write-Host "🌍 Starting local server at http://localhost:8000" -ForegroundColor Cyan
    Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
    Set-Location build
    python -m http.server 8000
}
