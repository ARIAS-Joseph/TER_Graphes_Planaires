@echo off
cd /d %~dp0

:: Essayer python dans le PATH d'abord
python --version >nul 2>&1
if %errorlevel% == 0 (
    python server.py
    pause
    exit /b
)

:: Essayer python3
python3 --version >nul 2>&1
if %errorlevel% == 0 (
    python3 server.py
    pause
    exit /b
)

:: Chercher un venv local dans le dossier courant
if exist ".venv\Scripts\python.exe" (
    .venv\Scripts\python.exe server.py
    pause
    exit /b
)

:: Python introuvable
echo.
echo ERREUR : Python n'est pas installe ou pas dans le PATH.
echo Installez Python depuis https://www.python.org/downloads/
echo et cochez "Add Python to PATH" lors de l'installation.
echo.
pause