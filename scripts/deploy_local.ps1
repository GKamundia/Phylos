# Local deployment script for testing Auspice dashboard (Windows PowerShell)
# Run this script to build and test your dashboard locally before deploying

Write-Host "🚀 Building RVF Nextstrain Dashboard locally..." -ForegroundColor Green

# Check if auspice is installed
try {
    auspice --version | Out-Null
    Write-Host "✅ Auspice is installed" -ForegroundColor Green
}
catch {
    Write-Host "❌ Auspice is not installed. Installing..." -ForegroundColor Red
    npm install -g auspice
}

# Check if JSON files exist
if (-not (Test-Path "auspice") -or -not (Get-ChildItem "auspice\*.json" -ErrorAction SilentlyContinue)) {
    Write-Host "❌ No Auspice JSON files found in .\auspice\" -ForegroundColor Red
    Write-Host "Please ensure your pipeline has generated the JSON files:" -ForegroundColor Yellow
    Write-Host "  - rvf_L.json"
    Write-Host "  - rvf_M.json" 
    Write-Host "  - rvf_S.json"
    exit 1
}

Write-Host "✅ Found Auspice JSON files:" -ForegroundColor Green
Get-ChildItem "auspice\*.json" | Select-Object Name, Length, LastWriteTime

# Validate JSON files
Write-Host "🔍 Validating JSON files..." -ForegroundColor Cyan
Get-ChildItem "auspice\*.json" | ForEach-Object {
    try {
        Get-Content $_.FullName -Raw | ConvertFrom-Json | Out-Null
        Write-Host "✅ $($_.Name) is valid" -ForegroundColor Green
    }
    catch {
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

auspice build --extend .\auspice --outputDir .\build

# Create index page
Write-Host "📄 Creating index page..." -ForegroundColor Cyan

$htmlContent = @"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RVF Nextstrain Dashboard</title>
    <style>
        body { 
            font-family: Arial, sans-serif;
            max-width: 800px; 
            margin: 2rem auto; 
            padding: 0 1rem;
            line-height: 1.6;
        }
        .dashboard-link {
            display: block;
            padding: 1rem;
            margin: 0.5rem 0;
            background: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            text-decoration: none;
            color: #495057;
        }
        .dashboard-link:hover {
            background: #e9ecef;
        }
        h1 { color: #2c3e50; }
        .segment-description { font-size: 0.9em; color: #6c757d; margin-top: 0.5rem; }
    </style>
</head>
<body>
    <h1>Rift Valley Fever Nextstrain Dashboard</h1>
    <p>Genomic surveillance and phylogenetic analysis of Rift Valley Fever virus segments.</p>
    
    <a href="./rvf_L" class="dashboard-link">
        <strong>RVF Segment L (Large)</strong>
        <div class="segment-description">RNA-dependent RNA polymerase gene</div>
    </a>
    
    <a href="./rvf_M" class="dashboard-link">
        <strong>RVF Segment M (Medium)</strong>
        <div class="segment-description">Glycoprotein genes and NSm protein</div>
    </a>
    
    <a href="./rvf_S" class="dashboard-link">
        <strong>RVF Segment S (Small)</strong>
        <div class="segment-description">Nucleocapsid protein and NSs protein genes</div>
    </a>
    
    <hr style="margin: 2rem 0;">
    <p><small>Built with Nextstrain | Source Code on GitHub</small></p>
</body>
</html>
"@

$htmlContent | Out-File -FilePath "build\index.html" -Encoding utf8

Write-Host "✅ Build complete! Static site generated in .\build\" -ForegroundColor Green

# Option to serve locally
$choice = Read-Host "🌐 Would you like to serve the dashboard locally for testing? (y/n)"
if ($choice -eq "y" -or $choice -eq "Y") {
    Write-Host "🌍 Starting local server at http://localhost:8000" -ForegroundColor Cyan
    Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
    Set-Location build
    python -m http.server 8000
}
