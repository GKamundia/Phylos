@echo off
echo ========================================
echo   RVF Nextstrain Dashboard Deployment
echo ========================================
echo.

REM Check if auspice directory exists
if not exist "auspice" (
    echo ERROR: auspice folder not found!
    echo Please make sure you're in the correct directory.
    pause
    exit /b 1
)

REM Check if JSON files exist
if not exist "auspice\rvf_L.json" (
    echo ERROR: rvf_L.json not found in auspice folder!
    pause
    exit /b 1
)

if not exist "auspice\rvf_M.json" (
    echo ERROR: rvf_M.json not found in auspice folder!
    pause
    exit /b 1
)

if not exist "auspice\rvf_S.json" (
    echo ERROR: rvf_S.json not found in auspice folder!
    pause
    exit /b 1
)

echo ✓ Found all required JSON files
echo.

REM Create simple deployment folder
if exist "simple-deploy" rmdir /s /q "simple-deploy"
mkdir "simple-deploy"
mkdir "simple-deploy\auspice"

REM Copy files
copy "index.html" "simple-deploy\"
copy "auspice\*.json" "simple-deploy\auspice\"

echo ✓ Created deployment folder with all files
echo.
echo Your dashboard is ready! You can now:
echo.
echo 1. GITHUB PAGES:
echo    - Commit and push these files to your GitHub repo
echo    - Enable GitHub Pages in repo settings
echo    - Your dashboard will be live at: https://gkamundia.github.io/Phylos/
echo.
echo 2. NETLIFY DROP:
echo    - Go to: https://app.netlify.com/drop
echo    - Drag and drop the 'simple-deploy' folder
echo    - Get instant live URL
echo.
echo 3. TEST LOCALLY:
echo    - Open 'simple-deploy\index.html' in your browser
echo.

set /p choice="Would you like to open the local version now? (y/n): "
if /i "%choice%"=="y" (
    start "" "simple-deploy\index.html"
)

echo.
echo Deployment files ready in: simple-deploy\
pause
