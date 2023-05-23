$ROOT_DIR = $ARGS[0]
$BUILD_DIR = $ARGS[1]

Write-Output $ROOT_DIR
Write-Output $BUILD_DIR

Write-Output "================================== Sphere Decay =================================="
Set-Location $BUILD_DIR\demos
Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_sphere_decay.exe $ROOT_DIR\demos --nogui"
Copy-Item "$ROOT_DIR\demos\sphere\postprocessing\plot_sphere_decay.plt" -Destination "$BUILD_DIR\demos\results\"
Copy-Item "$ROOT_DIR\demos\sphere\postprocessing\sphere_decay_comparison.txt" -Destination "$BUILD_DIR\demos\results\"
Set-Location $BUILD_DIR\demos\results
gnuplot plot_sphere_decay.plt

Write-Output "================================== RM3 Decay =================================="
Set-Location $BUILD_DIR\demos
Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_rm3_decay.exe $ROOT_DIR\demos --nogui"
Copy-Item "$ROOT_DIR\demos\rm3\postprocessing\plot_rm3_decay.plt" "$BUILD_DIR\demos\results\"
Copy-Item "$ROOT_DIR\demos\rm3\postprocessing\rm3_WECSim_decay.txt" -Destination "$BUILD_DIR\demos\results\"
Set-Location $BUILD_DIR\demos\results
gnuplot plot_rm3_decay.plt

Write-Output "================================== OSWEC Decay =================================="
Set-Location $BUILD_DIR\demos
Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_oswec_decay.exe $ROOT_DIR\demos --nogui"
Copy-Item "$ROOT_DIR\demos\oswec\postprocessing\plot_oswec_decay.plt" "$BUILD_DIR\demos\results\"
Copy-Item "$ROOT_DIR\demos\oswec\postprocessing\wecsim_oswec_decay.txt" -Destination "$BUILD_DIR\demos\results\"
Set-Location $BUILD_DIR\demos\results
gnuplot plot_oswec_decay.plt

Write-Output "================================== F3OF DT1 =================================="
Set-Location $BUILD_DIR\demos
Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_F3OF_DT1.exe $ROOT_DIR\demos --nogui"
Copy-Item "$ROOT_DIR\demos\f3of\postprocessing\F3OF_DT1_SURGE.plt" "$BUILD_DIR\demos\results\"
Copy-Item -Path "$ROOT_DIR\demos\f3of\reference_data\*" -Destination "$BUILD_DIR\demos\results\reference_data" -Recurse
Set-Location $BUILD_DIR\demos\results
gnuplot F3OF_DT1_SURGE.plt

Write-Output "================================== F3OF DT2 =================================="
Set-Location $BUILD_DIR\demos
Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_F3OF_DT2.exe $ROOT_DIR\demos --nogui"
Copy-Item "$ROOT_DIR\demos\f3of\postprocessing\F3OF_DT2_PITCH.plt" "$BUILD_DIR\demos\results\"
Set-Location $BUILD_DIR\demos\results
gnuplot F3OF_DT2_PITCH.plt

Write-Output "================================== F3OF DT3 =================================="
Set-Location $BUILD_DIR\demos
Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_F3OF_DT3.exe $ROOT_DIR\demos --nogui"
Copy-Item "$ROOT_DIR\demos\f3of\postprocessing\F3OF_DT3_PITCH.plt" "$BUILD_DIR\demos\results\"
Set-Location $BUILD_DIR\demos\results
gnuplot F3OF_DT3_PITCH.plt

Write-Output "================================== Sphere Reg Waves =================================="
Set-Location $BUILD_DIR\demos
Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_sphere_reg_waves.exe $ROOT_DIR\demos --nogui"
Write-Output "Sphere reg waves completed, no plot files yet"
# Copy-Item "$ROOT_DIR\demos\f3of\postprocessing\F3OF_DT3_PITCH.plt" "$BUILD_DIR\demos\results\"
# Copy-Item "$ROOT_DIR\demos\f3of\reference_data" -Destination "$BUILD_DIR\demos\results\"
# Set-Location $BUILD_DIR\demos\results
# gnuplot F3OF_DT3_PITCH.plt

Write-Output "================================== RM3 Reg Waves =================================="
Set-Location $BUILD_DIR\demos
Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_rm3_reg_waves.exe $ROOT_DIR\demos --nogui"
Write-Output "RM3 reg waves completed, no plot files yet"
# Copy-Item "$ROOT_DIR\demos\f3of\postprocessing\F3OF_DT3_PITCH.plt" "$BUILD_DIR\demos\results\"
# Copy-Item "$ROOT_DIR\demos\f3of\reference_data" -Destination "$BUILD_DIR\demos\results\"
# Set-Location $BUILD_DIR\demos\results
# gnuplot F3OF_DT3_PITCH.plt

Write-Output "================================== Sphere Irreg Waves =================================="
Set-Location $BUILD_DIR\demos
Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_sphere_irreg_waves.exe $ROOT_DIR\demos --nogui"
Write-Output "Sphere irreg waves completed, no plot files yet"
# Copy-Item "$ROOT_DIR\demos\f3of\postprocessing\F3OF_DT3_PITCH.plt" "$BUILD_DIR\demos\results\"
# Copy-Item "$ROOT_DIR\demos\f3of\reference_data" -Destination "$BUILD_DIR\demos\results\"
# Set-Location $BUILD_DIR\demos\results
# gnuplot F3OF_DT3_PITCH.plt
# Invoke-Expression -Command "$BUILD_DIR\demos\RelWithDebInfo\demo_sphere_irreg_waves.exe $ROOT_DIR\demos --nogui"

Set-Location $ROOT_DIR\demos