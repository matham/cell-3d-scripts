$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true


cell_3d_filter_cells --cells-path-glob `
-c "E:\imaging\analysis\2025*\**\MF1*640*_cells_data_located_model_V2.yml" `
-c "F:\imaging\analysis\2025*\**\MF1*640*_cells_data_located_model_V2.yml" `
-c "H:\imaging\analysis\2025*\**\MF1*640*_cells_data_located_model_V2.yml" `
-c "I:\imaging\analysis\2025*\**\MF1*640*_cells_data_located_model_V2.yml" `
-c "K:\imaging\analysis\2025*\**\MF1*640*_cells_data_located_model_V2.yml" `
-c "L:\imaging\analysis\2025*\**\MF1*640*_cells_data_located_model_V2.yml" `
-c "M:\imaging\analysis\2025*\**\MF1*640*_cells_data_located_model_V2.yml" `
-c "N:\imaging\analysis\2025*\**\MF1*640*_cells_data_located_model_V2.yml" `
-rp "D:\cell_metadata_plots\MF1_640_yidan" --output-plots-name "plots" --output-raw-name "raw_data.npz" `
apply-filter --subdir "intensity_g0sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussianS0" `
apply-filter --subdir "intensity_g1sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussianS1" `
apply-filter --subdir "intensity_g2sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussianS2" `
apply-filter --subdir "intensity_g3sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussianS3" `
apply-filter --subdir "intensity_g2sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_mean2" `
apply-filter --subdir "intensity_g3sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_mean3" `
apply-filter --subdir "intensity_g2sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_p2,50" `
apply-filter --subdir "intensity_g3sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity>=gaussian_p3,50" `
`
apply-filter --subdir "intensity_rel_g0sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussianS0" `
apply-filter --subdir "intensity_rel_g1sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussianS1" `
apply-filter --subdir "intensity_rel_g2sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussianS2" `
apply-filter --subdir "intensity_rel_g3sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussianS3" `
apply-filter --subdir "intensity_rel_g2sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussian_mean2" `
apply-filter --subdir "intensity_rel_g3sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussian_mean3" `
apply-filter --subdir "intensity_rel_g2sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussian_p2,50" `
apply-filter --subdir "intensity_rel_g3sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_mean_intensity_rel>=gaussian_p3,50" `
`
apply-filter --subdir "intensity_total_g0sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussianS0" `
apply-filter --subdir "intensity_total_g1sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussianS1" `
apply-filter --subdir "intensity_total_g2sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussianS2" `
apply-filter --subdir "intensity_total_g3sigma" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussianS3" `
apply-filter --subdir "intensity_total_g2sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussian_mean2" `
apply-filter --subdir "intensity_total_g3sigma_mean" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussian_mean3" `
apply-filter --subdir "intensity_total_g2sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussian_p2,50" `
apply-filter --subdir "intensity_total_g3sigma_med" -cf "log10:paor_d1_um5>1" -cf "log10:paor_d2_um5>1" -cf "log10:paor_d3_um5>1" -cf "log10:paor_volume_um3>1" -cf "log10:paor_intensity_total>=gaussian_p3,50"
