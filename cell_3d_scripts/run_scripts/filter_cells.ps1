$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true

$imaging_roots=@("E:\imaging\analysis", "G:\imaging\analysis", "I:\imaging\analysis", "M:\imaging\analysis")
foreach ($imaging_root in $imaging_roots) {

    $channels=@("561", "640")
    foreach ($channel in $channels) {
        $files = Get-ChildItem "$imaging_root\**$channel`_cells_data_located_model_V2.yml" -File -Recurse -Exclude @("good*", "bad*", "*light*")

        foreach ($yml in $files) {
            $root=Split-Path -Path "$yml" -Parent
            $name=Split-Path -Path "$yml" -Leaf

            $output="$root\filtered_cells\$channel"
            $good="good_$name"
            $bad="bad_$name"

            if (-not (Test-Path -Path "$good" -PathType Leaf)) {
                echo "Filtering $yml -> $output"
                cell_3d_filter_cells -c "$yml" -rp "$output" --output-cells-name "$good" --output-removed-cells-name "$bad" --output-plots-name "plots" `
                apply-filter --subdir "g0sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian0" `
                apply-filter --subdir "g1sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian1" `
                apply-filter --subdir "g2sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian2" `
                apply-filter --subdir "g3sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian3" `
                apply-filter --subdir "g2sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_mean2" `
                apply-filter --subdir "g3sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_mean3" `
                apply-filter --subdir "g2sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_p2,50"
            }
        }
    }
}
