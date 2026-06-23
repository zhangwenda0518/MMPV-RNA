@echo off
cd /d "%~dp0"
echo Running Metadata GUI self-test...
echo.
python test_all.py > test_output.log 2>&1
type test_output.log
echo.
echo Test output also saved to test_output.log
pause
