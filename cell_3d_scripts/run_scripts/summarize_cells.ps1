$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true

$output_vaa3d = $true
$output_analysis_grouped = $false
$output_hypothesis_grouped_short = $false

$atlases=@(
    @("allen_mouse_25um", "allen_mouse_25um_v1_2_CPL_v1", "allen_25_v1_2_cpl_v1", "analysis_groups_v02", "hypothesis_groups_a_v01"),
    @("kim_mouse_25um", "kim_mouse_25um_v1_1_CPL_v2", "kim_25_v1_1_cpl_v2", "analysis_groups_v02", "hypothesis_groups_a_v01"),
    @("kim_mouse_isotropic_20um", "kim_mouse_isotropic_20um_v1_0_CPL_v1", "kim_iso_20_v1_0_cpl_v1", "analysis_groups_v02", "hypothesis_groups_a_v01"),
    @("ccfv3augmented_mouse_25um", "ccfv3augmented_mouse_25um_v1_0_CPL_v1", "ccfv3_aug_25_v1_0_cpl_v1", "analysis_groups_v02", "hypothesis_groups_a_v01")
)

$imaging_roots=@("N:\imaging\analysis")#, "F:\imaging\analysis", "H:\imaging\analysis", "I:\imaging\analysis", "K:\imaging\analysis", "L:\imaging\analysis", "M:\imaging\analysis", "N:\imaging\analysis")
foreach ($imaging_root in $imaging_roots) {
    foreach ($atlas_pair in $atlases) {
        $atlas_registration=$atlas_pair[0]
        $atlas=$atlas_pair[1]
        $atlas_key=$atlas_pair[2]
        $analysis_groups=$atlas_pair[3]
        $hypothesis_groups=$atlas_pair[4]

        if ($output_analysis_grouped) {
            $splat_args=@("--merged-atlas-path", "D:\$atlas`_$analysis_groups.csv")
            $output_type="$analysis_groups"
        } elseif ($output_hypothesis_grouped_short) {
            $splat_args=@("--merged-atlas-path", "D:\$atlas`_$hypothesis_groups.csv")
            $output_type="$hypothesis_groups"
        } elseif ($output_vaa3d) {
            $splat_args=@("--output-vaa3d-format")
            $output_type="vaa3d"
        } else {
            $splat_args=@()
            $output_type="compact"
        }

        $channels=@("561", "640")
        foreach ($channel in $channels) {
            $files = Get-ChildItem "$imaging_root\**$channel`_cells_data_located_model_V2.yml" -File -Recurse -Exclude @(,"bad*") | Where-Object { $_.name -like "*MNP*"}

            foreach ($yml in $files) {
                $root=(Get-Item "$yml").DirectoryName
                $name=(Get-Item "$yml").Name
                $basename=(Get-Item "$yml").Basename

                if ($name -like "*good*") {
                    $threshold=(Get-Item "$yml").Directory.Name
                    $basename=$basename.substring(5)

                    $regions_vol_pat="$root\..\..\..\..\registration_v2\region_volumes*$atlas_registration`*.csv"
                    $output="$root\..\..\..\summaries\$threshold\$output_type\$atlas_key`_$threshold`_$basename`.csv"
                } else {
                    $regions_vol_pat="$root\..\registration_v2\region_volumes*$atlas_registration`*.csv"
                    $output="$root\summaries\all\$output_type\$atlas_key`_all_$basename`.csv"
                }

              if (-not (Test-Path -Path "$output" -PathType Leaf) -and (Test-Path -Path "$regions_vol_pat" -PathType Leaf)) {
                $regions_vol_path=(Get-ChildItem -Path "$regions_vol_pat" | Select-Object -First 1).FullName

                echo "Summarizing $yml -> $output"
                cell_3d_summarize_regions -c "$yml" --atlas-name "$atlas_key" --vaa3d-atlas-path "D:\$atlas`_vaa3d_atlas.csv" --regions-volume-path "$regions_vol_path" -o "$output" @splat_args
              }
            }
        }
    }
}
