$ErrorActionPreference = "Stop"

$atlas="kim_25um_v1_1"
$atlas_reg="kim_mouse_25um"

$imaging_roots=@("E:\imaging\fused", "G:\imaging\fused")
foreach ($imaging_root in $imaging_roots) {

	$channels=@("561", "640")
	foreach ($channel in $channels) {
		$files = Get-ChildItem "$imaging_root\**$channel`_cell_data_model*.yml" -File -Recurse -Exclude @("good*", "bad*")

		foreach ($yml in $files) {
		  $root=(Get-Item "$yml").DirectoryName
		  $name=(Get-Item "$yml").Name
		  $basename=(Get-Item "$yml").Basename

		  $output="$root\$atlas\$basename`_$atlas`.yml"
		  $regions_path="$root\..\registration\upsampled*$atlas_reg`*annotation.zarr.zip"

		  if (-not (Test-Path -Path "$output" -PathType Leaf) -and (Test-Path -Path "$regions_path" -PathType Leaf)) {
			$regions=(Get-ChildItem -Path "$regions_path" | Select-Object -First 1).FullName

			echo "Localizing $yml -> $output"
			cell_3d_locate_cells -c "$yml" -r "$regions" --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" -o "$output" --workers 6
		  }
		}
	}
}
