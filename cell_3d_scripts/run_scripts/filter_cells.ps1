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

            if (-not (Test-Path -Path "$good" -PathType Leaf)) {
                echo "Filtering $yml -> $output"
                cell_3d_filter_cells -c "$yml" -rp "$output" --output-cells-name "$good" --output-plots-name "plots" `
                apply-filter --subdir "intensity_g0sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussianS0" `
                apply-filter --subdir "intensity_g2sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussianS2" `
                apply-filter --subdir "intensity_g3sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussianS3" `
                apply-filter --subdir "intensity_g2sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_mean2" `
                apply-filter --subdir "intensity_g3sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_mean3" `
                apply-filter --subdir "intensity_g2sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_p2,50" `
                apply-filter --subdir "intensity_g3sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_p3,50" `
                apply-filter --subdir "intensity_biGaussianTrough" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=biGaussianTrough" `
                apply-filter --subdir "intensity_biGaussianS0" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=biGaussianS0" `
                apply-filter --subdir "intensity_biGaussianS-2" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=biGaussianS-2" `
                `
                apply-filter --subdir "intensity_rel_g0sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussianS0" `
                apply-filter --subdir "intensity_rel_g2sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussianS2" `
                apply-filter --subdir "intensity_rel_g3sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussianS3" `
                apply-filter --subdir "intensity_rel_g2sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussian_mean2" `
                apply-filter --subdir "intensity_rel_g3sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussian_mean3" `
                apply-filter --subdir "intensity_rel_g2sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussian_p2,50" `
                apply-filter --subdir "intensity_rel_g3sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussian_p3,50" `
                apply-filter --subdir "intensity_rel_biGaussianTrough" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=biGaussianTrough" `
                apply-filter --subdir "intensity_rel_biGaussianS0" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=biGaussianS0" `
                apply-filter --subdir "intensity_rel_biGaussianS-2" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=biGaussianS-2" `
                `
                apply-filter --subdir "intensity_total_g0sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussianS0" `
                apply-filter --subdir "intensity_total_g2sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussianS2" `
                apply-filter --subdir "intensity_total_g3sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussianS3" `
                apply-filter --subdir "intensity_total_g2sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussian_mean2" `
                apply-filter --subdir "intensity_total_g3sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussian_mean3" `
                apply-filter --subdir "intensity_total_g2sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussian_p2,50" `
                apply-filter --subdir "intensity_total_g3sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussian_p3,50" `
                apply-filter --subdir "intensity_total_biGaussianTrough" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=biGaussianTrough" `
                apply-filter --subdir "intensity_total_biGaussianS0" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=biGaussianS0" `
                apply-filter --subdir "intensity_total_biGaussianS-2" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=biGaussianS-2"
            }
        }
    }
}
