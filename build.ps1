<#
    HydroChrono Build Script
    
    Configures and builds HydroChrono. Dependency paths (Eigen, HDF5, Irrlicht)
    are automatically inherited from your Chrono build - just point to ChronoDir.
    
    Usage:
        .\build.ps1                     # Default build
        .\build.ps1 -Clean -Verbose     # Clean build with details
        .\build.ps1 -Package            # Build and create ZIP
    
    Configuration:
        Create build-config.json with: { "ChronoDir": "C:/path/to/chrono/build/cmake" }
#>

param(
    [ValidateSet("Release", "Debug", "RelWithDebInfo", "MinSizeRel")]
    [string]$BuildType = "Release",
    
    [switch]$Clean,
    [switch]$Verbose,
    [switch]$Package,
    [switch]$NoIrrlicht,    # Disable Irrlicht (auto-enabled if Chrono has it)
    [switch]$NoVSG,         # Disable VSG (auto-enabled if Chrono has it)
    [switch]$MoorDyn,       # Enable MoorDyn mooring coupling
    [switch]$NoDemos,
    [switch]$NoTests,
    [switch]$Help,
    
    [string]$ConfigPath = ".\build-config.json"
)

# Note: We don't use $ErrorActionPreference = "Stop" because it causes
# native command stderr output (like CMake warnings) to throw exceptions.
$SCRIPT_VERSION = "3.0"

# =============================================================================
# Output Helpers
# =============================================================================

function Write-Step { param([string]$msg) Write-Host "`n>> $msg" -ForegroundColor Cyan }
function Write-OK { param([string]$msg) Write-Host "   [OK] $msg" -ForegroundColor Green }
function Write-Warn { param([string]$msg) Write-Host "   [WARN] $msg" -ForegroundColor Yellow }
function Write-Fail { param([string]$msg) Write-Host "   [FAIL] $msg" -ForegroundColor Red }
function Write-Detail { param([string]$msg) Write-Host "   $msg" -ForegroundColor Gray }

# =============================================================================
# Help
# =============================================================================

if ($Help) {
    Write-Host "`nHydroChrono Build System v$SCRIPT_VERSION" -ForegroundColor Cyan
    Write-Host "======================================" -ForegroundColor Cyan
    Write-Host "`nBuilds HydroChrono. Dependencies are auto-detected from Chrono.`n"
    
    Write-Host "USAGE:" -ForegroundColor Yellow
    Write-Host "  .\build.ps1 [options]`n"
    
    Write-Host "OPTIONS:" -ForegroundColor Yellow
    Write-Host "  -BuildType <type>  Release|Debug|RelWithDebInfo|MinSizeRel (default: Release)"
    Write-Host "  -Clean             Remove build directory before building"
    Write-Host "  -Verbose           Show detailed build output"
    Write-Host "  -Package           Create distributable ZIP after building"
    Write-Host "  -NoIrrlicht        Disable Irrlicht (enabled by default if Chrono has it)"
    Write-Host "  -NoVSG             Disable VSG (enabled by default if Chrono has it)"
    Write-Host "  -MoorDyn           Enable MoorDyn mooring coupling"
    Write-Host "                     Requires: git submodule update --init extern/MoorDyn"
    Write-Host "  -NoDemos           Skip demo executables"
    Write-Host "  -NoTests           Skip test targets"
    Write-Host "  -ConfigPath <path> Custom config file (default: build-config.json)`n"
    
    Write-Host "EXAMPLES:" -ForegroundColor Yellow
    Write-Host "  .\build.ps1                    # Standard build"
    Write-Host "  .\build.ps1 -Clean -Verbose    # Clean build with details"
    Write-Host "  .\build.ps1 -Package           # Build and create ZIP`n"
    
    Write-Host "LONG TESTS:" -ForegroundColor Yellow
    Write-Host "  Set `$env:HYDROCHRONO_LONG_TESTS='1' before running ctest to use extended"
    Write-Host "  simulation durations for publication-quality regression reports.`n"
    
    Write-Host "CONFIG FILE:" -ForegroundColor Yellow
    Write-Host '  { "ChronoDir": "C:/path/to/chrono/build/cmake" }' -ForegroundColor Gray
    Write-Host ""
    Write-Host "DOCS:" -ForegroundColor Yellow
    Write-Host "  Full setup guide: docs/_main_pages/developer_docs/build_instructions.md`n"
    exit 0
}

# =============================================================================
# Load Configuration
# =============================================================================

Write-Host "`nHydroChrono Build v$SCRIPT_VERSION" -ForegroundColor Cyan
Write-Host "========================" -ForegroundColor Cyan

Write-Step "Loading configuration"

if (-not (Test-Path $ConfigPath)) {
    Write-Fail "Config file not found: $ConfigPath"
    Write-Host "`nCreate build-config.json:" -ForegroundColor Yellow
    Write-Host '  { "ChronoDir": "C:/path/to/chrono/build/cmake" }' -ForegroundColor Gray
    exit 1
}

try {
    $config = Get-Content -Path $ConfigPath -Raw | ConvertFrom-Json
    $ChronoDir = $config.ChronoDir
    $PythonRoot = $config.PythonRoot
} catch {
    Write-Fail "Failed to parse config: $($_.Exception.Message)"
    exit 1
}

if (-not $ChronoDir) {
    Write-Fail "ChronoDir not set in config file"
    exit 1
}

if (-not (Test-Path $ChronoDir)) {
    Write-Fail "Chrono not found: $ChronoDir"
    exit 1
}

Write-OK "Chrono: $ChronoDir"

# =============================================================================
# Check Chrono Configuration & Auto-Detect Features
# =============================================================================

Write-Step "Detecting Chrono features"

$chronoConfig = Join-Path $ChronoDir "chrono-config.cmake"
if (-not (Test-Path $chronoConfig)) {
    Write-Fail "chrono-config.cmake not found"
    exit 1
}

$chronoContent = Get-Content $chronoConfig -Raw

# Auto-detect available modules
$hasIrrlicht = $chronoContent -match 'Chrono_IRRLICHT_AVAILABLE\s+ON'
$hasVSG = $chronoContent -match 'Chrono_VSG_AVAILABLE\s+ON'
$hasHDF5 = $chronoContent -match 'CHRONO_HDF5_AVAILABLE\s+TRUE'

# Enable features by default if Chrono has them (unless user disabled)
$useIrrlicht = $hasIrrlicht -and (-not $NoIrrlicht)
$useVSG = $hasVSG -and (-not $NoVSG)

Write-Detail "Irrlicht: $(if($hasIrrlicht){if($useIrrlicht){'ON'}else{'OFF (disabled)'}}else{'not available'})"
Write-Detail "VSG: $(if($hasVSG){if($useVSG){'ON'}else{'OFF (disabled)'}}else{'not available'})"
Write-Detail "HDF5: $(if($hasHDF5){'available'}else{'not available'})"

if (-not $hasHDF5) {
    Write-Warn "Chrono built without HDF5 - build may fail"
}

# Eigen3: extract include path from Chrono config so FindEigen3.cmake can
# bypass Eigen3ConfigVersion.cmake (Eigen 5.x rejects the "3.3" request).
$eigen3Include = $null
if ($chronoContent -match 'Eigen3_DIR\s+"([^"]+)"' -and $Matches[1] -notmatch 'NOTFOUND') {
    $eigen3Root = Split-Path (Split-Path $Matches[1])
    $candidate = Join-Path $eigen3Root "include/eigen3"
    if (Test-Path (Join-Path $candidate "Eigen")) {
        $eigen3Include = $candidate
        Write-OK "Eigen3: $eigen3Include"
    } else {
        Write-Warn "Eigen3 headers not found at $candidate"
    }
}
if (-not $eigen3Include -and $chronoContent -match 'EIGEN3_INCLUDE_DIR\s+"([^"]+)"' -and $Matches[1] -notmatch 'NOTFOUND') {
    $candidate = $Matches[1]
    if (Test-Path (Join-Path $candidate "Eigen")) {
        $eigen3Include = $candidate
        Write-OK "Eigen3: $eigen3Include (via EIGEN3_INCLUDE_DIR)"
    } else {
        Write-Warn "Eigen3 headers not found at $candidate"
    }
}

# =============================================================================
# Check CMake
# =============================================================================

Write-Step "Checking tools"

try {
    $cmakeVer = (cmake --version | Select-Object -First 1) -replace 'cmake version ',''
    Write-OK "CMake $cmakeVer"
} catch {
    Write-Fail "CMake not found"
    exit 1
}

# MoorDyn submodule check
if ($MoorDyn -and -not (Test-Path ".\extern\MoorDyn\CMakeLists.txt")) {
    Write-Fail "MoorDyn submodule not found at extern/MoorDyn"
    Write-Host "   Run: git submodule update --init extern/MoorDyn" -ForegroundColor Yellow
    exit 1
}

# =============================================================================
# Clean
# =============================================================================

if ($Clean -and (Test-Path ".\build")) {
    Write-Step "Cleaning"
    Remove-Item -Recurse -Force .\build
    Write-OK "Build directory removed"
}

# =============================================================================
# Configure
# =============================================================================

Write-Step "Configuring"

if (-not (Test-Path ".\build")) {
    New-Item -ItemType Directory -Path build | Out-Null
}

$cmakeArgs = @(
    "-S", ".",
    "-B", "build",
    "-G", "Visual Studio 17 2022",
    "-A", "x64",
    "-DChrono_DIR=`"$($ChronoDir -replace '\\','/')`"",
    "-DHYDROCHRONO_ENABLE_YAML_RUNNER=ON",
    "-DHYDROCHRONO_ENABLE_TESTS=$(if($NoTests){'OFF'}else{'ON'})",
    "-DHYDROCHRONO_ENABLE_DEMOS=$(if($NoDemos){'OFF'}else{'ON'})",
    "-DHYDROCHRONO_ENABLE_IRRLICHT=$(if($useIrrlicht){'ON'}else{'OFF'})",
    "-DHYDROCHRONO_ENABLE_VSG=$(if($useVSG){'ON'}else{'OFF'})",
    "-DCMAKE_BUILD_TYPE=$BuildType",
    "-DHYDROCHRONO_ENABLE_MOORDYN=$(if($MoorDyn){'ON'}else{'OFF'})"
)

# Add Eigen3 include path if extracted from Chrono config
if ($eigen3Include) {
    $cmakeArgs += "-DEIGEN3_INCLUDE_DIR=`"$($eigen3Include -replace '\\','/')`""
}

# Add Python if specified
if ($PythonRoot -and (Test-Path $PythonRoot)) {
    $cmakeArgs += "-DPython3_ROOT_DIR=`"$($PythonRoot -replace '\\','/')`""
}

# Get Irrlicht path from Chrono if needed
if ($useIrrlicht) {
    if ($chronoContent -match 'Irrlicht_INCLUDE_DIR\s+"([^"]+)"') {
        $cmakeArgs += "-DIrrlicht_INCLUDE_DIR=`"$($Matches[1])`""
    }
    if ($chronoContent -match 'set\(Irrlicht_LIBRARY\s+"([^"]+)"') {
        $irrlichtRoot = Split-Path (Split-Path (Split-Path $Matches[1])) 
        $cmakeArgs += "-DIrrlicht_ROOT=`"$($irrlichtRoot -replace '\\','/')`""
    }
}

Write-Detail "Build: $BuildType | Tests: $(if($NoTests){'OFF'}else{'ON'}) | Demos: $(if($NoDemos){'OFF'}else{'ON'}) | MoorDyn: $(if($MoorDyn){'ON'}else{'OFF'})"

if ($Verbose) {
    cmake @cmakeArgs
} else {
    $null = cmake @cmakeArgs 2>&1
}

if ($LASTEXITCODE -ne 0) {
    Write-Fail "CMake failed"
    Write-Host "   Run with -Verbose for details" -ForegroundColor Yellow
    exit 1
}

Write-OK "Configuration complete"

# =============================================================================
# Build
# =============================================================================

Write-Step "Building"

$t0 = Get-Date

if ($Verbose) {
    cmake --build build --config $BuildType
} else {
    $null = cmake --build build --config $BuildType 2>&1
}

$buildSecs = [math]::Round(((Get-Date) - $t0).TotalSeconds, 1)

if ($LASTEXITCODE -ne 0) {
    Write-Fail "Build failed"
    Write-Host "   Run with -Verbose for details" -ForegroundColor Yellow
    exit 1
}

Write-OK "Built in $buildSecs seconds"

# =============================================================================
# Copy Third-Party DLLs (not handled by Chrono's DLL copy command)
# =============================================================================

$binPath = ".\build\bin\$BuildType"

if ($useIrrlicht) {
    # Get Irrlicht DLL path from Chrono config
    if ($chronoContent -match 'Irrlicht_DLL\s+"([^"]+)"') {
        $irrlichtDll = $Matches[1]
        if (Test-Path $irrlichtDll) {
            $destDll = Join-Path $binPath "Irrlicht.dll"
            if (-not (Test-Path $destDll)) {
                Copy-Item $irrlichtDll $destDll -Force
                Write-OK "Copied Irrlicht.dll"
            }
        }
    }
}

if ($useVSG) {
    # Get VSG DLL directory from Chrono config (set by find_package(vsg))
    # VSG_DLL_DIR is set by Chrono when it finds VSG
    $vsgDllDir = $null
    
    if ($chronoContent -match 'VSG_DLL_DIR\s+([^\s\)]+)') {
        $vsgDllDir = $Matches[1].Trim('"')
    }
    
    # Fallback: derive from vsg_DIR (lib/cmake/vsg -> bin)
    if (-not $vsgDllDir -or -not (Test-Path $vsgDllDir)) {
        if ($chronoContent -match 'vsg_DIR\s+"([^"]+)"') {
            $vsgCmakeDir = $Matches[1]
            # Go up from lib/cmake/vsg to root, then into bin
            $vsgRoot = Split-Path (Split-Path (Split-Path $vsgCmakeDir))
            $vsgDllDir = Join-Path $vsgRoot "bin"
        }
    }
    
    if ($vsgDllDir -and (Test-Path $vsgDllDir)) {
        $vsgDlls = Get-ChildItem -Path $vsgDllDir -Filter "*.dll" -ErrorAction SilentlyContinue
        $copiedCount = 0
        foreach ($dll in $vsgDlls) {
            $destDll = Join-Path $binPath $dll.Name
            if (-not (Test-Path $destDll)) {
                Copy-Item $dll.FullName $destDll -Force
                $copiedCount++
            }
        }
        if ($copiedCount -gt 0) {
            Write-OK "Copied $copiedCount VSG DLLs from $vsgDllDir"
        }
    } else {
        Write-Warn "VSG enabled but could not find VSG DLL directory"
    }
}

# Copy HDF5 and its dependencies (szip, zlib) from vcpkg
if ($chronoContent -match 'HDF5_DIR\s+"([^"]+)"') {
    $hdf5Root = Split-Path (Split-Path $Matches[1])
    $hdf5BinDir = Join-Path $hdf5Root "bin"
    if (Test-Path $hdf5BinDir) {
        $hdf5Dlls = Get-ChildItem -Path $hdf5BinDir -Include "hdf5*.dll","szip*.dll","zlib*.dll" -Recurse -ErrorAction SilentlyContinue
        $copiedCount = 0
        foreach ($dll in $hdf5Dlls) {
            $destDll = Join-Path $binPath $dll.Name
            if (-not (Test-Path $destDll)) {
                Copy-Item $dll.FullName $destDll -Force
                $copiedCount++
            }
        }
        if ($copiedCount -gt 0) {
            Write-OK "Copied $copiedCount HDF5/compression DLLs from $hdf5BinDir"
        }
    }
}

# Verify MoorDyn DLL is in the output directory
if ($MoorDyn) {
    $moordynDll = Join-Path $binPath "moordyn.dll"
    if (Test-Path $moordynDll) {
        Write-OK "MoorDyn DLL: moordyn.dll"
    } else {
        # CMake places moordyn.dll in bin/<config> alongside other targets,
        # but fall back to searching the build tree if not found.
        $found = Get-ChildItem -Path ".\build" -Filter "moordyn.dll" -Recurse -ErrorAction SilentlyContinue |
            Select-Object -First 1
        if ($found) {
            Copy-Item $found.FullName $moordynDll -Force
            Write-OK "MoorDyn DLL: copied from $($found.Directory.Name)"
        } else {
            Write-Warn "MoorDyn enabled but moordyn.dll not found in build tree"
        }
    }
}

# =============================================================================
# Verify
# =============================================================================
$exe = Join-Path $binPath "run_hydrochrono.exe"

if (Test-Path $exe) {
    $size = [math]::Round((Get-Item $exe).Length / 1MB, 1)
    Write-OK "run_hydrochrono.exe ($size MB)"
}

# =============================================================================
# Package
# =============================================================================

if ($Package) {
    Write-Step "Packaging"
    
    $prefix = ".\build\install"
    
    if ($Verbose) {
        cmake --install build --config $BuildType --prefix $prefix
    } else {
        $null = cmake --install build --config $BuildType --prefix $prefix 2>&1
    }
    
    if ($LASTEXITCODE -ne 0) {
        Write-Fail "Install failed"
        exit 1
    }
    
    Push-Location build
    if ($Verbose) { cpack -C $BuildType } else { $null = cpack -C $BuildType 2>&1 }
    $packResult = $LASTEXITCODE
    Pop-Location
    
    if ($packResult -ne 0) {
        Write-Fail "CPack failed"
        exit 1
    }
    
    $pkg = Get-ChildItem ".\build\HydroChrono-*.zip" -ErrorAction SilentlyContinue | Select-Object -First 1
    if ($pkg) {
        Write-OK "$($pkg.Name) ($([math]::Round($pkg.Length/1MB, 1)) MB)"
    }
}

# =============================================================================
# Done
# =============================================================================

Write-Host "`n========================================" -ForegroundColor Green
Write-Host "BUILD SUCCESSFUL" -ForegroundColor Green
Write-Host "========================================`n" -ForegroundColor Green

Write-Host "Output: $binPath" -ForegroundColor Cyan
Write-Host ""
Write-Host "Tests:  ctest -C $BuildType -L regression --test-dir build" -ForegroundColor Gray
Write-Host "        ctest -C $BuildType -L unit       --test-dir build" -ForegroundColor Gray
Write-Host "        Add -V for verbose output, --output-on-failure for failures only" -ForegroundColor DarkGray
Write-Host ""
Write-Host "Long:   `$env:HYDROCHRONO_LONG_TESTS='1'" -ForegroundColor Gray
Write-Host "        ctest -C $BuildType -L regression --test-dir build" -ForegroundColor Gray
Write-Host "        Runs with extended simulation durations" -ForegroundColor DarkGray
Write-Host ""
Write-Host "Report: python tests/regression/utilities/generate_report.py --build-dir build --pdf" -ForegroundColor Gray
Write-Host "        Generates regression test report (markdown + PDF) in build/bin/report/" -ForegroundColor DarkGray
Write-Host "        Requires: pip install pypandoc  (or pandoc on PATH)" -ForegroundColor DarkGray
Write-Host ""
