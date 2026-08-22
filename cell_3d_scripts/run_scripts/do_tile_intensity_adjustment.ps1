$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true

# $imaging_roots=@("E:\imaging\analysis", "F:\imaging\analysis", "H:\imaging\analysis", "I:\imaging\analysis", "K:\imaging\analysis", "L:\imaging\analysis", "M:\imaging\analysis", "N:\imaging\analysis")
$imaging_roots=@("K:\imaging\analysis", "L:\imaging\analysis", "M:\imaging\analysis", "N:\imaging\analysis")
foreach ($imaging_root in $imaging_roots) {

    $channels=@("561", "640")
    foreach ($channel in $channels) {
        $files = Get-ChildItem "$imaging_root\**$channel`*.tif" -File -Recurse

        foreach ($ch_tiff in $files) {
            $root=Split-Path -Path "$ch_tiff" -Parent
            $name=Split-Path -Path "$root" -Leaf

            $tiff488="$root\$name`_BS_488.tif"
            $cells_path="$root\results\cellfinder_count_v2\$name`_$channel`_cells_data_located_model_V2.yml"
            $projection_path="$root\results\cellfinder_count_v2\$name`_488_z_mean.tif"
            $adjustment_path="$root\results\cellfinder_count_v2\$name`_488_tile_intensity_adjustment.tif"

            if (Test-Path -Path "$adjustment_path") {
                $splat_args=@()
            } else {
                $splat_args=@("-i",  "$tiff488")
            }

            echo "Setting tile intensity adjustment in $cells_path"
            cell_3d_tile_intensity_adjustment --skip-processed-cells @splat_args -ic "$cells_path" -oc "$cells_path" -p "$projection_path" -a "$adjustment_path" --num-tiles 4 6 --sampling-step 10 --n-patch-pixels 200
        }
    }
}
