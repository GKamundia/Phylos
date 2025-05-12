@echo off
:: filepath: c:\Users\Anarchy\Documents\Data_Science\NextStrain\rvf-nextstrain\scripts\run_scheduled_build.bat
:: Script for Windows Task Scheduler to run Nextstrain builds
:: Usage: run_scheduled_build.bat [data_only] [pathogen]
::   data_only  - If set to "true", only runs data acquisition
::   pathogen   - The pathogen to build (defaults to active_pathogen in master_config.yaml)

setlocal enabledelayedexpansion

:: Set working directory to the project root
cd /d "%~dp0.."

:: Parse arguments
set DATA_ONLY=false
set PATHOGEN=

if "%1"=="true" set DATA_ONLY=true
if not "%2"=="" set PATHOGEN=%2

:: Timestamp for logs
for /f "tokens=2-4 delims=/ " %%a in ('date /t') do (
  set DATE=%%c-%%a-%%b
)
for /f "tokens=1-3 delims=: " %%a in ('time /t') do (
  set TIME=%%a-%%b-%%c
)
set TIMESTAMP=%DATE%_%TIME%
set LOG_FILE=logs\scheduled_build_%TIMESTAMP%.log

:: Create logs directory if it doesn't exist
if not exist logs mkdir logs

:: Log start of run
echo Starting Nextstrain build at %date% %time% > %LOG_FILE%

:: Activate the Conda environment if it exists
if exist "C:\Users\%USERNAME%\Miniconda3\Scripts\activate.bat" (
    echo Activating Conda environment... >> %LOG_FILE%
    call C:\Users\%USERNAME%\Miniconda3\Scripts\activate.bat nextstrain
) else (
    echo Warning: Conda environment activation skipped. Make sure nextstrain-cli is available. >> %LOG_FILE%
)

:: Run data acquisition step
echo Running data acquisition step... >> %LOG_FILE%
nextstrain build . --configfile config/master_config.yaml --until prepare_metadata >> %LOG_FILE% 2>&1
if %ERRORLEVEL% neq 0 (
    echo ERROR: Data acquisition failed with code %ERRORLEVEL% >> %LOG_FILE%
    goto :exit_script
)

:: If data_only is true, skip the rest
if "%DATA_ONLY%"=="true" (
    echo Data acquisition completed. Skipping analysis as requested. >> %LOG_FILE%
    goto :exit_script
)

:: Run the full build
echo Running full Nextstrain build... >> %LOG_FILE%
if "%PATHOGEN%"=="" (
    nextstrain build . --configfile config/master_config.yaml >> %LOG_FILE% 2>&1
) else (
    :: Set active_pathogen temporarily
    python scripts/set_active_pathogen.py %PATHOGEN% >> %LOG_FILE% 2>&1
    nextstrain build . --configfile config/master_config.yaml >> %LOG_FILE% 2>&1
    :: Restore original value from backup
    python scripts/restore_config.py >> %LOG_FILE% 2>&1
)

if %ERRORLEVEL% neq 0 (
    echo ERROR: Nextstrain build failed with code %ERRORLEVEL% >> %LOG_FILE%
) else (
    echo Nextstrain build completed successfully at %date% %time% >> %LOG_FILE%
    
    :: Generate build summary
    python scripts/generate_build_summary.py >> %LOG_FILE% 2>&1

    :: Create a backup after successful run
    echo Creating backup of critical data... >> %LOG_FILE%
    python scripts\backup_data.py backup --name "scheduled_win" >> %LOG_FILE% 2>&1
    if %ERRORLEVEL% neq 0 (
        echo WARNING: Backup creation failed with code %ERRORLEVEL% >> %LOG_FILE%
    ) else (
        echo Backup created successfully >> %LOG_FILE%
    )
)

:exit_script
echo Build process finished at %date% %time% >> %LOG_FILE%

:: Send email notification if configured
if exist scripts\send_notification.py (
    echo Sending notification... >> %LOG_FILE%
    python scripts\send_notification.py --log-file %LOG_FILE% >> %LOG_FILE% 2>&1
)

exit /b %ERRORLEVEL%