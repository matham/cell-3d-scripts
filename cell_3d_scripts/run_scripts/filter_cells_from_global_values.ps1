$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true

$imaging_roots=@("N:\imaging\analysis")#, "F:\imaging\analysis", "H:\imaging\analysis", "I:\imaging\analysis", "K:\imaging\analysis", "L:\imaging\analysis", "M:\imaging\analysis", "N:\imaging\analysis")
foreach ($imaging_root in $imaging_roots) {

    $channels=@("561", "640")
    foreach ($channel in $channels) {
        $files = Get-ChildItem "$imaging_root\*MF1*$channel`_cells_data_located_model_V2.yml" -File -Recurse -Exclude @("good*", "bad*", "*light*") | Where-Object { $_.Directory -like "*2026*"}

        foreach ($yml in $files) {
            $root=Split-Path -Path "$yml" -Parent
            $name=Split-Path -Path "$yml" -Leaf

            $output="$root\filtered_cells\$channel"
            $good="good_$name"

            if (-not (Test-Path -Path "$good" -PathType Leaf)) {
                echo "Filtering $yml -> $output"
                cell_3d_filter_cells -c "$yml" -rp "$output" --output-cells-name "$good" --output-plots-name "plots" `
                apply-filter --subdir "intensity_g0sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=3.8413141149100887" `
                apply-filter --subdir "intensity_g2sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=4.272289056405588" `
                apply-filter --subdir "intensity_g3sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=4.487776527153338" `
                apply-filter --subdir "intensity_g2sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=4.4924685067627"
            }
        }
    }
}
