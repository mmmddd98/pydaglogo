@ECHO OFF

set SOURCEDIR=docs
set BUILDDIR=docs\_build

if "%SPHINXBUILD%"=="" (
    where sphinx-build >NUL 2>NUL
    if errorlevel 1 (
        set SPHINXBUILD=python -m sphinx
    ) else (
        set SPHINXBUILD=sphinx-build
    )
)

if "%1"=="" goto help

%SPHINXBUILD% -b %1 %SOURCEDIR% %BUILDDIR%\%1
goto end

:help
echo.
echo Usage:
echo   make.bat html
echo.
echo This runs:
echo   sphinx-build -b html docs docs\_build\html
echo.

:end
