$ErrorActionPreference = "Stop"

$imaging_roots=@("E:\imaging\fused", "G:\imaging\fused")
foreach ($imaging_root in $imaging_roots) {

	$channels=@("561", "640")
	foreach ($channel in $channels) {
		$files = Get-ChildItem "$imaging_root\**$channel`*cell_classification_model*.xml" -File -Recurse

		foreach ($xml in $files) {
		  $root=Split-Path -Path "$xml\..\.." -Parent -Resolve
		  $name=Split-Path -Path "$root" -Leaf

		  $tiff="$root\$name`_BS_$channel`.tif"
		  $res_root="$root\results\cellfinder_count"
		  $yml="$res_root\$name`_$channel`_cell_data_model_V1.yml"

		  if ((Test-Path -Path "$xml" -PathType Leaf) -and -not (Test-Path -Path "$yml" -PathType Leaf)) {
			echo "Analyzing $xml -> $yml"
			cell_meta_3d -s "$tiff" -c "$xml" -o "$yml" --voxel-size 4 2.03 2.03 --batch-size 132 --max-workers 6 --cube-size 100 50 50 --initial-center-search-radius 4 2 2 --initial-center-search-volume 4 2 2 --lateral-intensity-algorithm area_margin --lateral-max-radius 16 --lateral-decay-length 10 --lateral-decay-algorithm gaussian --axial-intensity-algorithm center_line --axial-max-radius 32 --axial-decay-length 32 --axial-decay-algorithm gaussian
		  }
		}
	}
}
