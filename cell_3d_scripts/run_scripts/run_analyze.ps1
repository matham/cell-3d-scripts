$thresholds=@("no", "peak", "mean", "p25", "p50", "p75")
foreach ($threshold in $thresholds) {
	$fmts=@("vaa3d", "analysis_groups_v02", "hypothesis_groups_a_v01")
	foreach ($fmt in $fmts) {
		if ($fmt -eq "vaa3d") {
			$output_fmt="atlas_leafs"
			$splat_args4=@("--only-leafs")
		} else {
			$output_fmt="$fmt"
			$splat_args4=@()
		}

		$groups=@("core", "main", "all")
		foreach ($group in $groups) {
			$models=@("condition", "condition_genotype_sex")
			foreach ($model in $models) {
				if ($group -eq "core") {
					$splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
				} elseif ($group -eq "main") {
					$splat_args2=@("-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
				} else {
					$splat_args2=@()
				}

				if ($model -eq "condition") {
					$splat_args3=@("-gk", "Condition")
				} else {
					$splat_args3=@("-gk", "Condition", "-gk", "Genotype", "-gk", "Sex")
				}

				$measures=@("cell_density", "cell_count")
				foreach ($measure in $measures) {
					if ($fmt -eq "vaa3d") {
						$splat_args1=@("-oam", "D:\$measure`_collated_simple\$model\$group\metadata_$model`-$group`.csv")
					} else {
						$splat_args1=@()
					}
					$measure_key="Cell density (mm3)"
					if ($measure -eq "cell_count") {
						$measure_key="Total cells"
					}


					cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" -c 561 -c 640 -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*kim_25um_v1_1.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "$measure_key" -ir 8 --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" -oas "D:\$measure`_collated_simple\$model\$group\$output_fmt\$measure`-$output_fmt`-$threshold`_threshold-{channel}-$model`-$group`.csv"  @splat_args1 @splat_args2 @splat_args3 @splat_args4

					cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" -c 561 -c 640 -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*kim_25um_v1_1.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "$measure_key" -ir 8 --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" -oa "D:\$measure`_analyzed\$model\$group\$fmt\deseq2\$measure`_deseq2_$fmt`-$threshold`_threshold-{channel}-$model`-$group`.csv" @splat_args2 @splat_args3 @splat_args4 --deseq2

					cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" -c 561 -c 640 -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*kim_25um_v1_1.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "$measure_key" -ir 8 --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" -oa "D:\$measure`_analyzed\$model\$group\$fmt\ttest\$measure`_ttest_$fmt`-$threshold`_threshold-{channel}-$model`-$group`.csv" @splat_args2 @splat_args3 @splat_args4 --t-test --t-test-false-discovery
				}

				cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" -c 561 -c 640 -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*kim_25um_v1_1.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "Cell density (mm3)" -ir 8 --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" -or "D:\collated_wide\$model\$group\$fmt\wide_raw_$fmt`-$threshold`_threshold-{channel}-$model`-$group`.csv" -og "D:\collated_wide\$model\$group\$fmt\wide_averages_$fmt`-$threshold`_threshold-{channel}-$model`-$group`.csv" @splat_args2 @splat_args3
			}
		}
	}
}
