$thresholds=@("no", "peak", "mean", "p5", "p25", "p50", "p75")
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

        $models=@("core_cond", "core_cond_sex_buffer", "core_cond_trap2", "core_cond_wt", "core_cond_wt_sex_buffer", "all_cond", "all_cond_sex_buffer")
        foreach ($model in $models) {
            if ($model -eq "core_cond") {
                $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                $splat_args3=@("-gk", "Condition", "-md", "Condition")
            } elseif ($model -eq "core_cond_sex_buffer") {
                $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-gk", "Buffer", "-md", "     + Condition*Sex + Condition*buffer")
            } elseif ($model -eq "all_cond") {
                $splat_args2=@("-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                $splat_args3=@("-gk", "Condition", "-md", "Condition")
            } elseif ($model -eq "all_cond_sex_buffer") {
                $splat_args2=@("-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-gk", "Buffer", "-md", "Condition + Condition*Sex + Condition*buffer")
            } elseif ($model -eq "core_cond_trap2") {
                $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP_WT")
                $splat_args3=@("-gk", "Condition", "-gk", "Genotype", "-md", "Condition + Condition*Genotype")
            } elseif ($model -eq "core_cond_wt") {
                $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Genotype=FosTRAP2")
                $splat_args3=@("-gk", "Condition", "-gk", "Genotype", "-md", "Condition + Condition*Genotype")
            } elseif ($model -eq "core_cond_wt_sex_buffer") {
                $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Genotype=FosTRAP2")
                $splat_args3=@("-gk", "Condition", "-gk", "Genotype", "-gk", "Sex", "-gk", "Buffer", "-md", "Condition + Condition*Genotype + Condition*Sex + Condition*buffer")
            } else {
                exit
            }

            $measures=@("cell_count", "cell_density", "total_volume", "total_volume_mm3")
            foreach ($measure in $measures) {
                if ($fmt -eq "vaa3d") {
                    $splat_args1=@("-oam", "D:\$measure`_deseq2_simple\$model\metadata_$model`.csv")
                } else {
                    $splat_args1=@()
                }
                $splat_args5=@()
                if ($measure -eq "cell_count") {
                    $measure_key="Total cells"
                    $splat_args5=@("--deseq2-norm", "region_volume")
                } elseif ($measure -eq "cell_density") {
                    $measure_key="Cell density (mm3)"
                } elseif ($measure -eq "total_volume") {
                    $measure_key="Total voxels"
                } elseif ($measure -eq "total_volume_mm3") {
                    $measure_key="Total volume (mm3)"
                }

                cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" -c 561 -c 640 -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*kim_25um_v1_1.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "$measure_key" -ir 8 --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" --output-path-all-stat-deseq2 "D:\$measure`_deseq2_simple\$model\$output_fmt\data-$measure`-$output_fmt`-{channel}-$model`.csv" --output-path-all-norm-deseq2 "D:\$measure`_deseq2_simple\$model\$output_fmt\norm-$measure`-$output_fmt`-{channel}-$model`.csv" @splat_args1 @splat_args2 @splat_args3 @splat_args4 -vk "Total volume (mm3)"

                cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" -c 561 -c 640 -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*kim_25um_v1_1.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "$measure_key" -ir 8 --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" -oa "D:\$measure`_analyzed\$model\$fmt\deseq2\$measure`_deseq2_$fmt`-$threshold`_threshold-{channel}-$model`.csv" @splat_args2 @splat_args3 @splat_args4 --deseq2 @splat_args5 -vk "Total volume (mm3)"

                cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" -c 561 -c 640 -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*kim_25um_v1_1.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "$measure_key" -ir 8 --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" -oa "D:\$measure`_analyzed\$model\$fmt\ttest\$measure`_ttest_$fmt`-$threshold`_threshold-{channel}-$model`.csv" @splat_args2 @splat_args3 @splat_args4 --t-test --t-test-false-discovery -vk "Total volume (mm3)"
            }

            cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" -c 561 -c 640 -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*kim_25um_v1_1.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "Cell density (mm3)" -ir 8 --vaa3d-atlas-path "D:\Kim 25um v1.1, CPL v1 vaa3d atlas.csv" -or "D:\collated_wide\$model\$fmt\wide_raw_$fmt`-$threshold`_threshold-{channel}-$model`.csv" -og "D:\collated_wide\$model\$fmt\wide_averages_$fmt`-$threshold`_threshold-{channel}-$model`.csv" @splat_args2 @splat_args3 -vk "Total volume (mm3)"
        }
    }
}
