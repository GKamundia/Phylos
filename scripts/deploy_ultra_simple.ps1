# Ultra-Simple Deployment (No PowerShell HTML issues)
Write-Host "🚀 Creating simple deployment..." -ForegroundColor Green

# Check if JSON files exist
if (-not (Test-Path "auspice")) {
    Write-Host "❌ No auspice folder found" -ForegroundColor Red
    exit 1
}

$jsonFiles = Get-ChildItem "auspice\*.json" -ErrorAction SilentlyContinue
if (-not $jsonFiles) {
    Write-Host "❌ No JSON files found in auspice folder" -ForegroundColor Red
    exit 1
}

Write-Host "✅ Found JSON files:" -ForegroundColor Green
$jsonFiles | ForEach-Object { Write-Host "  - $($_.Name)" }

# Create deployment directory
if (Test-Path "deploy") {
    Remove-Item "deploy" -Recurse -Force
}
New-Item -ItemType Directory -Path "deploy" -Force | Out-Null

# Copy JSON files
Write-Host "📁 Copying JSON files..." -ForegroundColor Cyan
Copy-Item "auspice\*.json" -Destination "deploy\" -Force

# Create simple index.html using separate file
$htmlFile = "deploy\index.html"

# Write HTML content line by line to avoid PowerShell parsing issues
@"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RVF Nextstrain Dashboard</title>
    <style>
        body { font-family: Arial, sans-serif; max-width: 800px; margin: 2rem auto; padding: 1rem; background: #f5f5f5; }
        .container { background: white; padding: 2rem; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .link { display: block; padding: 1rem; margin: 1rem 0; background: #e3f2fd; border: 2px solid #2196f3; border-radius: 8px; text-decoration: none; color: #0d47a1; }
        .link:hover { background: #bbdefb; }
        h1 { color: #1565c0; text-align: center; }
        .desc { font-size: 0.9em; color: #666; margin-top: 0.5rem; font-style: italic; }
        .info { background: #fff3cd; border: 1px solid #ffeaa7; border-radius: 4px; padding: 1rem; margin: 1rem 0; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🦠 Rift Valley Fever Nextstrain Dashboard</h1>
        
        <div class="info">
            <strong>How to view:</strong> Click links to open in Nextstrain web app, or download JSON files for local viewing.
        </div>
        
        <p><strong>Genomic surveillance and phylogenetic analysis of Rift Valley Fever virus segments.</strong></p>
        
"@ | Out-File -FilePath $htmlFile -Encoding utf8 -NoNewline

# Add download links for each JSON file
foreach ($file in $jsonFiles) {
    $segment = $file.BaseName -replace "rvf_", ""
    $segmentName = switch ($segment) {
        "L" { "Large" }
        "M" { "Medium" }  
        "S" { "Small" }
        default { $segment }
    }
    
    $linkHtml = @"
        
        <a href="$($file.Name)" class="link" download>
            <strong>🧬 RVF Segment $segment ($segmentName)</strong>
            <div class="desc">Download JSON file for Nextstrain visualization</div>
        </a>
"@
    Add-Content -Path $htmlFile -Value $linkHtml -Encoding utf8
}

# Add footer
$footer = @"
        
        <div style="text-align: center; margin: 2rem 0; padding: 1rem; background: #e8f5e8; border-radius: 4px;">
            <p><strong>How to use:</strong></p>
            <p>1. Download JSON files above</p>
            <p>2. Go to <a href="https://auspice.us" target="_blank">auspice.us</a></p>
            <p>3. Drag and drop JSON files to view</p>
        </div>
        
        <div style="text-align: center; margin-top: 2rem; padding-top: 1rem; border-top: 1px solid #ddd; color: #666;">
            <p>Built with Nextstrain | <a href="https://github.com/GKamundia/Phylos">Source Code</a></p>
        </div>
    </div>
</body>
</html>
"@

Add-Content -Path $htmlFile -Value $footer -Encoding utf8

Write-Host "✅ Deployment created in .\deploy\" -ForegroundColor Green
Write-Host "📂 Files created:" -ForegroundColor Cyan
Get-ChildItem "deploy" | ForEach-Object { Write-Host "  - $($_.Name) ($([math]::Round($_.Length/1KB, 1)) KB)" }

Write-Host ""
Write-Host "🌐 Next steps:" -ForegroundColor Yellow
Write-Host "1. Upload deploy\ folder contents to GitHub Pages" -ForegroundColor White
Write-Host "2. Or drag deploy\ folder to netlify.com" -ForegroundColor White
Write-Host "3. Or host on any web server" -ForegroundColor White
