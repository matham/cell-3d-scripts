$ErrorActionPreference = "Stop"

$imaging_roots=@("E:\imaging\fused", "G:\imaging\fused")
foreach ($imaging_root in $imaging_roots) {

	$channels=@("561", "640")
	foreach ($channel in $channels) {
		$files = Get-ChildItem "$imaging_root\**$channel`*.tif" -File -Recurse

		$tiled_threshold="1"
		if ("$channel" -eq "561") {
			$tiled_threshold="1.25"
		}

		foreach ($ch_tiff in $files) {
			$root=Split-Path -Path "$ch_tiff" -Parent
			$name=Split-Path -Path "$root" -Leaf

			$tiff488="$root\$name`_BS_488.tif"
			$cells_dir="$root\cells_$channel"
			$res_root="$root\results\cellfinder_count"
			$xml="$res_root\$name`_$channel`_cell_classification_model_V1.xml"

			if (-not (Test-Path -Path "$res_root")) {
				mkdir "$res_root"
			}
			if (-not (Test-Path -Path "$cells_dir")) {
				mkdir "$cells_dir"
			}
			if (-not (Test-Path -Path "$xml")) {
				echo "Detecting cells in $ch_tiff"
				brainmapper -s "$ch_tiff" -b "$tiff488" -o "$cells_dir" -v 4 2.03 2.03 --no-register --no-analyse --no-figures --max-cluster-size 10000 --soma-diameter 8 --ball-xy-size 8 --ball-z-size 8 --ball-overlap-fraction 0.65 --log-sigma-size 0.2 --threshold 0.5 --tiled-threshold "$tiled_threshold" --tiled-threshold-tile-size 5 --soma-spread-factor 4 --detection-batch-size 1 --trained-model "D:\models\model_20250925_V1.keras" --classification-batch-size 128 --pin-memory --orientation psl --norm-channels --norm-sampling 32 --classification-max-workers 6 --detect-coi
				cp "$cells_dir\points\cell_classification.xml" "$xml"
				cp "$cells_dir\*.log" "$res_root"
			}
		}
	}
}
