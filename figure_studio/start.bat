@echo off
rem Figure Studio 启动器（用 %~dp0 取脚本自身目录，避免在文件里写中文路径）
cd /d "%~dp0"
echo Starting Figure Studio on http://127.0.0.1:8770 ...
start "" http://127.0.0.1:8770
"C:\Python312\python.exe" app.py --port 8770
pause
