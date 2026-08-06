$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true

$imaging_roots=@("I:\imaging\analysis", "K:\imaging\analysis")
# $imaging_roots=@("E:\imaging\analysis", "F:\imaging\analysis", "H:\imaging\analysis")
foreach ($imaging_root in $imaging_roots) {
    $channels=@("561", "640")
    foreach ($channel in $channels) {
        $files = Get-ChildItem "$imaging_root\**$channel`*cell_classification_model*.yml" -File -Recurse

        foreach ($yml in $files) {
            $root=(Get-Item "$yml").DirectoryName
            $tiff_root=Split-Path -Path "$yml\..\.." -Parent -Resolve
            $name=Split-Path -Path "$tiff_root" -Leaf

            $output="$root\$name`_$channel`_cells_located_model_V2.yml"
            $regions_ccf="$root\..\registration_v2\upsampled*ccfv3augmented_mouse_25um_annotation.zarr.zip"
            $regions_kim="$root\..\registration_v2\upsampled*kim_mouse_25um_annotation.zarr.zip"
            $regions_kim_iso="$root\..\registration_v2\upsampled*kim_mouse_isotropic_20um_annotation.zarr.zip"
            $regions_allen="$root\..\registration_v2\upsampled*allen_mouse_25um_annotation.zarr.zip"

            if (-not (Test-Path -Path "$output" -PathType Leaf) -and (Test-Path -Path "$regions_ccf" -PathType Leaf) -and (Test-Path -Path "$regions_kim" -PathType Leaf) -and (Test-Path -Path "$regions_kim_iso" -PathType Leaf) -and (Test-Path -Path "$regions_allen" -PathType Leaf)) {
                $regions_ccf=(Get-ChildItem -Path "$regions_ccf" | Select-Object -First 1).FullName
                $regions_kim=(Get-ChildItem -Path "$regions_kim" | Select-Object -First 1).FullName
                $regions_kim_iso=(Get-ChildItem -Path "$regions_kim_iso" | Select-Object -First 1).FullName
                $regions_allen=(Get-ChildItem -Path "$regions_allen" | Select-Object -First 1).FullName

                echo "Localizing $yml -> $output"
                cell_3d_locate_cells -c "$yml" -r "$regions_ccf" "$regions_kim" "$regions_kim_iso" "$regions_allen" --vaa3d-atlas-path "D:\ccfv3augmented_mouse_25um_v1_0_CPL_v1_vaa3d_atlas.csv" "D:\kim_mouse_25um_v1_1_CPL_v2_vaa3d_atlas.csv" "D:\kim_mouse_isotropic_20um_v1_0_CPL_v1_vaa3d_atlas.csv" "D:\allen_mouse_25um_v1_2_CPL_v1_vaa3d_atlas.csv" --atlas-name "ccfv3_aug_25_v1_0_cpl_v1" "kim_25_v1_1_cpl_v2" "kim_iso_20_v1_0_cpl_v1" "allen_25_v1_2_cpl_v1" -o "$output" --workers 6 --exclude-non-cells
            }
        }
    }
}
