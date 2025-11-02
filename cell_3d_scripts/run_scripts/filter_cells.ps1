$ErrorActionPreference = "Stop"

$atlas="kim_25um_v1_1"

$imaging_roots=@("E:\imaging\fused", "G:\imaging\fused")
foreach ($imaging_root in $imaging_roots) {

	$channels=@("561", "640")
	foreach ($channel in $channels) {
		$files = Get-ChildItem "$imaging_root\**$channel`_cell_data_model*$atlas`.yml" -File -Recurse -Exclude @("good*", "bad*")

		foreach ($yml in $files) {
			$intensities=@("peak", "mean", "p25", "p50", "p75")
			foreach ($intensity in $intensities) {
				$root=Split-Path -Path "$yml" -Parent
				$name=Split-Path -Path "$yml" -Leaf

				$output="$root\filtered_cells\$intensity\$channel"
				$good="$output\good_$name"
				$bad="$output\bad_$name"

				if (-not (Test-Path -Path "$good" -PathType Leaf)) {
					echo "Filtering $yml -> $output"
					cell_3d_filter_cells -c "$yml" -o "$good" --output-removed-cells-path "$bad" -plots "$output" -cf "r_xy>0" -cf "r_z>0" -cf "r_z_max_std<=4" -cf "r_xy_max_std<=2" -cf "center_intensity>=$intensity"
				}
			}
		}
	}
}
