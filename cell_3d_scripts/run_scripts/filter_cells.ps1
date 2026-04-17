$ErrorActionPreference = "Stop"

$imaging_roots=@("E:\imaging\analysis", "G:\imaging\analysis", "I:\imaging\analysis", "M:\imaging\analysis")
$atlases=@("kim_25um_v1_1")#, "allen_25um_v1_1")

foreach ($imaging_root in $imaging_roots) {

    foreach ($atlas in $atlases) {
        $channels=@("561", "640")
        foreach ($channel in $channels) {
            $files = Get-ChildItem "$imaging_root\**$channel`_cell_data_model*$atlas`.yml" -File -Recurse -Exclude @("good*", "bad*")

            foreach ($yml in $files) {
                $filters=@("peak")#, "fixed", "peak_int_xy", "peak_int_xyz", "mean", "p5", "p25", "p50", "p75")
                foreach ($filter in $filters) {
                    $root=Split-Path -Path "$yml" -Parent
                    $name=Split-Path -Path "$yml" -Leaf

                    $output="$root\filtered_cells\$filter\$channel"
                    $good="$output\good_$name"
                    $bad="$output\bad_$name"

                    if ($filter -eq "peak_int_xy") {
                        $splat_args=@("-cf", "r_xy>peak")
                        $intensity="peak"
                    } elseif ($filter -eq "peak_int_xyz") {
                        $splat_args=@("-cf", "r_xy>peak", "-cf", "r_z>peak")
                        $intensity="peak"
                    } elseif ($filter -eq "fixed") {
                        $splat_args=@("-cf", "center_intensity>=5605")
                        $intensity="peak"
                    } elseif ($filter -eq "peak_r") {
                        $splat_args=@("-cf", "intensity_ratio>=peak")
                        $intensity="peak"
                    } else {
                        $splat_args=@()
                        $intensity="$filter"
                    }

                    if (-not (Test-Path -Path "$good" -PathType Leaf)) {
                        echo "Filtering $yml -> $output"
                        cell_3d_filter_cells -c "$yml" -o "$good" --output-removed-cells-path "$bad" -plots "$output" -cf "r_xy>0" -cf "r_z>0" -cf "r_z_max_std<=4" -cf "r_xy_max_std<=2" -cf "center_intensity>=$intensity" @splat_args
                    }
                }
            }
        }
    }
}
