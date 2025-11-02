param ([switch] $output_vaa3d=$false, [switch] $output_grouped=$false, [switch] $output_grouped_short=$false, $threshold="no")
$ErrorActionPreference = "Stop"

$merged_atlas_path="D:\Kim mouse 25um_v1.1 Analysis Groups v02.csv"
$merged_short_atlas_path="D:\Kim mouse 25um_v1.1 Hypothesis Group A v01 List.csv"
$atlas="kim_25um_v1_1"
$atlas_reg="kim_mouse_25um"

if ($output_grouped) {
	$splat_args=@("--merged-atlas-path", "$merged_atlas_path")
	$output_prefix="region_summary_analysis_groups_v02"
} elseif ($output_grouped_short) {
	$splat_args=@("--merged-atlas-path", "$merged_short_atlas_path")
	$output_prefix="region_summary_hypothesis_groups_a_v01"
} elseif ($output_vaa3d) {
	$splat_args=@("--output-vaa3d-format")
	$output_prefix="region_summary_vaa3d"
} else {
	$splat_args=@()
	$output_prefix="region_summary_bare"
}

if ($threshold -eq "no") {
	$splat_threshold_args=@()
	$threshold_label="no_threshold"
} else {
	$splat_threshold_args=@("-cf", "center_intensity>=$threshold")
	$threshold_label="$threshold`_threshold"
}

$imaging_roots=@("E:\imaging\fused", "G:\imaging\fused")
foreach ($imaging_root in $imaging_roots) {

	$channels=@("561", "640")
	foreach ($channel in $channels) {
		$files = Get-ChildItem "$imaging_root\**$channel`*_cell_data_model*$atlas`.yml" -File -Recurse -Exclude @("good*", "bad*")

		foreach ($yml in $files) {
			$root=(Get-Item "$yml").DirectoryName
			$name=(Get-Item "$yml").Name
			$basename=(Get-Item "$yml").Basename

			$regions_vol_pat="$root\..\..\registration\region_volumes*$atlas_reg`*.csv"
			$output="$root\summaries\$output_prefix\$output_prefix`_$threshold_label`_$basename`.csv"

		  if (-not (Test-Path -Path "$output" -PathType Leaf) -and (Test-Path -Path "$regions_vol_pat" -PathType Leaf)) {
			$regions_vol_path=(Get-ChildItem -Path "$regions_vol_pat" | Select-Object -First 1).FullName

			echo "Summarizing $yml -> $output"
			cell_3d_summarize_regions -c "$yml" --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" --regions-volume-path "$regions_vol_path" -o "$output" -cf "r_xy>0" -cf "r_z>0" -cf "r_z_max_std<=4" -cf "r_xy_max_std<=2" @splat_threshold_args @splat_args
		  }
		}
	}
}
