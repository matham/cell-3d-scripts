$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true

$imaging_roots=@("I:\imaging\analysis", "K:\imaging\analysis")
# $imaging_roots=@("E:\imaging\analysis", "F:\imaging\analysis", "H:\imaging\analysis")
foreach ($imaging_root in $imaging_roots) {

    $channels=@("561", "640")
    foreach ($channel in $channels) {
        $files = Get-ChildItem "$imaging_root\**$channel`*.tif" -File -Recurse

        $tiled_threshold="1.75"
        $threshold="1"

        foreach ($ch_tiff in $files) {
            $root=Split-Path -Path "$ch_tiff" -Parent
            $name=Split-Path -Path "$root" -Leaf

            $tiff488="$root\$name`_BS_488.tif"
            $cells_dir="$root\cells_$channel"
            $res_root="$root\results\cellfinder_count_v2"
            $yml="$res_root\$name`_$channel`_cell_classification_model_V2.yml"

            if (-not (Test-Path -Path "$res_root")) {
                mkdir "$res_root"
            }
            if (-not (Test-Path -Path "$cells_dir")) {
                mkdir "$cells_dir"
            }
            if (-not (Test-Path -Path "$yml")) {
                echo "Detecting cells in $ch_tiff"
                brainmapper -s "$ch_tiff" -b "$tiff488" -o "$cells_dir" -v 4 2.03 2.03 --no-register --no-analyse --no-figures --max-cluster-size 100000 --soma-diameter 8 --ball-xy-size 6 --ball-z-size 8 --ball-overlap-fraction 0.65 --log-sigma-size 0.2 --threshold "$threshold" --tiled-threshold "$tiled_threshold" --tiled-threshold-tile-size 3 --soma-spread-factor 5 --detection-batch-size 1 --trained-model "D:\models\model_20260717_V2.keras" --classification-batch-size 128 --pin-memory --orientation psl --norm-channels --norm-sampling 32 --classification-max-workers 6 --detect-coi
                cp "$cells_dir\points\cell_classification.yml" "$yml"
                cp "$cells_dir\*.log" "$res_root"
            }
        }
    }
}
