$thresholds=@("no", "peak", "mean", "p25", "p50", "p75")
foreach ($threshold in $thresholds) {
	$fmts=@("bare", "analysis_groups_v02", "hypothesis_groups_a_v01")
	foreach ($fmt in $fmts) {
		cell_3d_collate_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" -c 561 -c 640 -k "Image_name" -cols 7 -o "D:\collated_down\$threshold`_threshold\collated_down_$fmt`_$threshold`_threshold.csv" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*kim_25um_v1_1.csv"
	}
}
