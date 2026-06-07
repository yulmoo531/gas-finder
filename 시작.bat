@echo off
chcp 65001 > nul
title Gas Finder

echo.
echo [1/2] Installing packages...
python -m pip install flask requests --quiet --quiet

echo [2/2] Starting server...
echo.
start "" http://localhost:5000
python app.py

pause
