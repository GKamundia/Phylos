# Simple Static Deployment Script (No Auspice Installation Required)
# This creates a simple deployment that works with any web server

Write-Host "🚀 Creating simple static deployment..." -ForegroundColor Green

# Check if JSON files exist
if (-not (Test-Path "auspice") -or -not (Get-ChildItem "auspice\*.json" -ErrorAction SilentlyContinue)) {
    Write-Host "❌ No Auspice JSON files found in .\auspice\" -ForegroundColor Red
    exit 1
}

# Create deployment directory
if (Test-Path "deploy") {
    Remove-Item "deploy" -Recurse -Force
}
New-Item -ItemType Directory -Path "deploy" -Force | Out-Null

# Copy JSON files
Write-Host "📁 Copying Auspice JSON files..." -ForegroundColor Cyan
Copy-Item "auspice\*.json" -Destination "deploy\" -Force

# Create index.html that redirects to nextstrain.org
$indexContent = @"
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
            background-color: #f8f9fa;
        }
        .container {
            background: white;
            padding: 2rem;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .dashboard-link {
            display: block;
            padding: 1rem;
            margin: 1rem 0;
            background: #e3f2fd;
            border: 2px solid #2196f3;
            border-radius: 8px;
            text-decoration: none;
            color: #0d47a1;
            transition: all 0.3s;
        }
        .dashboard-link:hover {
            background: #bbdefb;
            transform: translateY(-2px);
        }
        h1 { 
            color: #1565c0; 
            text-align: center;
            border-bottom: 3px solid #2196f3;
            padding-bottom: 1rem;
        }
        .segment-description { 
            font-size: 0.9em; 
            color: #424242; 
            margin-top: 0.5rem; 
            font-style: italic;
        }
        .info-box {
            background: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 4px;
            padding: 1rem;
            margin: 1rem 0;
            color: #856404;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>🦠 Rift Valley Fever Nextstrain Dashboard</h1>
        
        <div class="info-box">
            <strong>📱 How to view:</strong> Click the links below to open each segment in the Nextstrain web application.
        </div>
        
        <p><strong>Genomic surveillance and phylogenetic analysis of Rift Valley Fever virus segments.</strong></p>
        <p>Select a viral segment below to explore the phylogenetic data:</p>
        
        <a href="https://nextstrain.org/fetch/YOUR_GITHUB_PAGES_URL/rvf_L.json" class="dashboard-link" target="_blank">
            <strong>🧬 RVF Segment L (Large)</strong>
            <div class="segment-description">RNA-dependent RNA polymerase gene - The largest genomic segment</div>
        </a>
        
        <a href="https://nextstrain.org/fetch/YOUR_GITHUB_PAGES_URL/rvf_M.json" class="dashboard-link" target="_blank">
            <strong>🧬 RVF Segment M (Medium)</strong>
            <div class="segment-description">Glycoprotein genes (Gn, Gc) and NSm protein - Medium genomic segment</div>
        </a>
        
        <a href="https://nextstrain.org/fetch/YOUR_GITHUB_PAGES_URL/rvf_S.json" class="dashboard-link" target="_blank">
            <strong>🧬 RVF Segment S (Small)</strong>
            <div class="segment-description">Nucleocapsid protein and NSs protein genes - Smallest genomic segment</div>
        </a>
        
        <div style="text-align: center; margin: 2rem 0; padding: 1rem; background: #e8f5e8; border-radius: 4px;">
            <p><strong>Alternative:</strong> Download JSON files and drag them onto <a href="https://auspice.us" target="_blank">auspice.us</a></p>
            <p style="margin-top: 1rem;">
                <a href="rvf_L.json" download>📥 Download L segment</a> | 
                <a href="rvf_M.json" download>📥 Download M segment</a> | 
                <a href="rvf_S.json" download>📥 Download S segment</a>
            </p>
        </div>
        
        <div style="text-align: center; margin-top: 2rem; padding-top: 1rem; border-top: 1px solid #ddd; color: #666;">
            <p><small>Built with <strong>Nextstrain</strong> | 
            <a href="https://github.com/GKamundia/Phylos" target="_blank">View Source Code</a></small></p>
        </div>
    </div>
</body>
</html>
"@

$indexContent | Out-File -FilePath "deploy\index.html" -Encoding utf8

Write-Host "✅ Simple deployment ready in .\deploy\" -ForegroundColor Green
Write-Host "📂 Contents:" -ForegroundColor Cyan
Get-ChildItem "deploy" | Select-Object Name, Length

Write-Host "`n🌐 Deployment options:" -ForegroundColor Yellow
Write-Host "1. Upload 'deploy' folder to any web hosting service" -ForegroundColor White
Write-Host "2. Use GitHub Pages: copy contents to docs/ folder" -ForegroundColor White
Write-Host "3. Use Netlify: drag and drop the 'deploy' folder" -ForegroundColor White

Write-Host "`n💡 Note: Links will use nextstrain.org/fetch to display your data" -ForegroundColor Cyan
Write-Host "Replace 'YOUR_GITHUB_PAGES_URL' in index.html with your actual URL" -ForegroundColor Yellow
