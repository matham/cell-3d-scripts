$ErrorActionPreference = "Stop"

$imaging_roots=@("E:\imaging\analysis", "G:\imaging\analysis", "I:\imaging\analysis", "M:\imaging\analysis")
$atlases=@(,@("kim_mouse_25um", "kim_25um_v1_1"))#, @("allen_mouse_25um", "allen_25um_v1_1"))

foreach ($imaging_root in $imaging_roots) {
    foreach ($atlas_pair in $atlases) {
        $atlas_reg=$atlas_pair[0]
        $atlas=$atlas_pair[1]

        $channels=@("561", "640")
        foreach ($channel in $channels) {
            $files = Get-ChildItem "$imaging_root\**$channel`_cell_data_model*.yml" -File -Recurse -Exclude @("good*", "bad*")

            foreach ($yml in $files) {
              $root=(Get-Item "$yml").DirectoryName
              $basename=(Get-Item "$yml").Basename

              $output="$root\$atlas\$basename`_$atlas`.yml"
              $regions_path="$root\..\registration\upsampled*$atlas_reg`*annotation.zarr.zip"

              if (-not (Test-Path -Path "$output" -PathType Leaf) -and (Test-Path -Path "$regions_path" -PathType Leaf)) {
                $regions=(Get-ChildItem -Path "$regions_path" | Select-Object -First 1).FullName

                echo "Localizing $yml -> $output"
                cell_3d_locate_cells -c "$yml" -r "$regions" --vaa3d-atlas-path "D:\$atlas`, CPL v1 vaa3d atlas.csv" -o "$output" --workers 6 --exclude-non-cells
              }
            }
        }
    }
}

foreach ($imaging_root in $imaging_roots) {
    foreach ($atlas_pair in $atlases) {
        $atlas_reg=$atlas_pair[0]
        $atlas=$atlas_pair[1]

        $channels=@("561", "640")
        foreach ($channel in $channels) {
            $files = Get-ChildItem "$imaging_root\**$channel`*cell_classification_model*.xml" -File -Recurse

            foreach ($xml in $files) {
              $root=(Get-Item "$xml").DirectoryName
              $tiff_root=Split-Path -Path "$xml\..\.." -Parent -Resolve
              $name=Split-Path -Path "$tiff_root" -Leaf

              $output="$root\$atlas\$name`_$channel`_non_cells_model_V1_$atlas`.yml"
              $regions_path="$root\..\registration\upsampled*$atlas_reg`*annotation.zarr.zip"

              if (-not (Test-Path -Path "$output" -PathType Leaf) -and (Test-Path -Path "$regions_path" -PathType Leaf)) {
                $regions=(Get-ChildItem -Path "$regions_path" | Select-Object -First 1).FullName

                echo "Localizing $xml -> $output"
                cell_3d_locate_cells -c "$xml" -r "$regions" --vaa3d-atlas-path "D:\$atlas`, CPL v1 vaa3d atlas.csv" -o "$output" --workers 6 --exclude-cells
              }
            }
        }
    }
}
