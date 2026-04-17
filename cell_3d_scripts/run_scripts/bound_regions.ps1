$ErrorActionPreference = "Stop"

$atlases=@("kim_mouse_25um", "allen_mouse_25um")
$imaging_roots=@("E:\imaging\analysis", "G:\imaging\analysis", "I:\imaging\analysis", "M:\imaging\analysis")

foreach ($imaging_root in $imaging_roots) {

    foreach ($atlas in $atlases) {
        $zarrs = Get-ChildItem "$imaging_root\**upsampled_*_488_$atlas`_annotation.zarr.zip" -File -Recurse

        foreach ($zarr in $zarrs) {
          $root=(Get-Item "$zarr").DirectoryName
          $basename=(Get-Item "$zarr").Basename
          # remove .zarr
          $basename=$basename.substring(0,$basename.length-5)

          $output="$root\$basename`_regions_cuboid`.yml"

          if (-not (Test-Path -Path "$output" -PathType Leaf)) {

            echo "Bounding $zarr -> $output"
            python "C:\Users\CPLab\PycharmProjects\cell_3d_scripts\cell_3d_scripts\bound_regions.py" -r "$zarr" -b "$output" --workers 12
          }
        }
    }
}
