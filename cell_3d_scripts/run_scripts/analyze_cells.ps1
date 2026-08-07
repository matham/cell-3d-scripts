$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true

$imaging_roots=@("E:\imaging\analysis", "F:\imaging\analysis", "H:\imaging\analysis", "I:\imaging\analysis", "K:\imaging\analysis")
foreach ($imaging_root in $imaging_roots) {

    $channels=@("561", "640")
    foreach ($channel in $channels) {
        $files = Get-ChildItem "$imaging_root\**$channel`*cells_located_model*.yml" -File -Recurse

        foreach ($yml_in in $files) {
          $root=Split-Path -Path "$yml_in\..\.." -Parent -Resolve
          $name=Split-Path -Path "$root" -Leaf

          $tiff="$root\$name`_BS_$channel`.tif"
          $res_root="$root\results\cellfinder_count_v2"
          $yml_out="$res_root\$name`_$channel`_cells_data_located_model_V2.yml"
          $seg_out="$res_root\$name`_$channel`_cells_segmentation_model_V2.h5"

          if ((Test-Path -Path "$yml_in" -PathType Leaf) -and -not (Test-Path -Path "$yml_out" -PathType Leaf)) {
            echo "Analyzing $yml_in -> $yml_out"
            cell_meta_3d -s "$tiff" -c "$yml_in" -o "$yml_out" --voxel-size 4 2.03 2.03 --batch-size 128 --max-workers 12 --cube-size 108 92 92 --initial-center-search-radius 20 20 20 --lateral-intensity-algorithm area_margin --lateral-max-radius 16 --lateral-decay-length 10 --lateral-decay-algorithm gaussian --axial-intensity-algorithm center_line --axial-max-radius 32 --axial-decay-length 32 --axial-decay-algorithm gaussian --segmentation-path "$seg_out" --axial-decay-fraction 0.6666 --lateral-decay-fraction 0.6666 --seg-super-voxel 4 2 2 --seg-decay-fraction 0.6666 --seg-padding-factor 2.0
          }
        }
    }
}
