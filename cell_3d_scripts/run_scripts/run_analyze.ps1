$atlases=@("kim_25um_v1_1")#, "allen_25um_v1_1")
foreach ($atlas in $atlases) {

    $thresholds=@("fixed")#, "peak_int_xyz", "peak", "peak_int_xy", "no", "mean", "p5", "p25", "p50", "p75")
    foreach ($threshold in $thresholds) {
        $fmts=@("vaa3d")#, "analysis_groups_v02", "hypothesis_groups_a_v01")
        foreach ($fmt in $fmts) {
            if ($fmt -eq "vaa3d") {
                $output_fmt="atlas_leafs"
                $splat_args4=@("--only-leafs")
            } else {
                $output_fmt="$fmt"
                $splat_args4=@()
            }

    #         $models=@("core_cond_wt", "core_cond_wt_sex", "core_cond_wt_sex_buffer", "core_cond_sex", "all_cond_sex", "core_cond", "core_cond_sex_buffer", "core_cond_trap2", "all_cond", "all_cond_sex_buffer")
    #         $models=@("cond_paired_wo12", "cond_paired_wb")
            $models=@("all_sex", "core", "core_sex", "all")#, "core_merged", "core_merged_sex", "all_merged", "all_merged_sex")
            foreach ($model in $models) {
    #             $channel_splat_args=@("-c", "561", "-c", "640")
                $channel_splat_args=@("-c", "561")
                if ($model -eq "core_cond") {
                    $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                    $splat_args3=@("-gk", "Condition", "-md", "Condition")
                } elseif ($model -eq "core_cond_sex_buffer") {
                    $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-gk", "Buffer", "-md", "Condition + Condition*Sex + Condition*buffer")
                } elseif ($model -eq "core_cond_sex") {
                    $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-md", "Condition + Condition*Sex")
                } elseif ($model -eq "all_cond") {
                    $splat_args2=@("-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                    $splat_args3=@("-gk", "Condition", "-md", "Condition")
                } elseif ($model -eq "all_cond_sex_buffer") {
                    $splat_args2=@("-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-gk", "Buffer", "-md", "Condition + Condition*Sex + Condition*buffer")
                } elseif ($model -eq "all_cond_sex") {
                    $splat_args2=@("-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP2", "-eg", "Genotype=FosTRAP_WT")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-md", "Condition + Condition*Sex")
                } elseif ($model -eq "core_cond_trap2") {
                    $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Condition=TMT_IHC", "-eg", "Genotype=FosTRAP_WT")
                    $splat_args3=@("-gk", "Condition", "-gk", "Genotype", "-md", "Condition + Condition*Genotype")
                } elseif ($model -eq "core_cond_wt") {
                    $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Genotype=FosTRAP2")
                    $splat_args3=@("-gk", "Condition", "-md", "Condition")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "core_cond_wt_sex_buffer") {
                    $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Genotype=FosTRAP2")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-gk", "Buffer", "-md", "Condition + Condition*Sex + Condition*buffer")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "core_cond_wt_sex") {
                    $splat_args2=@("-eg", "Condition=2MBA", "-eg", "Condition=IAMM", "-eg", "Genotype=FosTRAP2")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-md", "Condition + Condition*Sex")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "cond_W_O1") {
                    $splat_args2=@("-eg", "Exposure=B", "-eg", "Exposure=W", "-eg", "Exposure=W_O2", "-eg", "Cohort_of_xpt=old")
                    $splat_args3=@("-gk", "Condition", "-md", "Condition")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "cond_W_O1_sex") {
                    $splat_args2=@("-eg", "Exposure=B", "-eg", "Exposure=W", "-eg", "Exposure=W_O2", "-eg", "Cohort_of_xpt=old")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-md", "Condition*Sex")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "cond_W_O2") {
                    $splat_args2=@("-eg", "Exposure=B", "-eg", "Exposure=W", "-eg", "Exposure=W_O1", "-eg", "Cohort_of_xpt=old")
                    $splat_args3=@("-gk", "Condition", "-md", "Condition")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "cond_W_O2_sex") {
                    $splat_args2=@("-eg", "Exposure=B", "-eg", "Exposure=W", "-eg", "Exposure=W_O1", "-eg", "Cohort_of_xpt=old")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-md", "Condition*Sex")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "cond_W") {
                    $splat_args2=@("-eg", "Exposure=B", "-eg", "Exposure=W_O1", "-eg", "Exposure=W_O2", "-eg", "Cohort_of_xpt=old")
                    $splat_args3=@("-gk", "Condition", "-md", "Condition")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "cond_W_sex") {
                    $splat_args2=@("-eg", "Exposure=B", "-eg", "Exposure=W_O1", "-eg", "Exposure=W_O2", "-eg", "Cohort_of_xpt=old")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-md", "Condition*Sex")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "cond_B") {
                    $splat_args2=@("-eg", "Exposure=W", "-eg", "Exposure=W_O1", "-eg", "Exposure=W_O2", "-eg", "Cohort_of_xpt=old")
                    $splat_args3=@("-gk", "Condition", "-md", "Condition")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "cond_B_sex") {
                    $splat_args2=@("-eg", "Exposure=W", "-eg", "Exposure=W_O1", "-eg", "Exposure=W_O2", "-eg", "Cohort_of_xpt=old")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-md", "Condition*Sex")
                    $channel_splat_args=@("-c", "640")
                } elseif ($model -eq "cond_paired_wo12") {
                    $splat_args2=@("-eg", "Exposure=W", "-eg", "Exposure=B", "-eg", "Cohort_of_xpt=old", "-eg", "Condition=Immuno")
                    $splat_args3=@("-gk", "ExpCond", "-md", "ExpCond")
                    $channel_splat_args=@("-c", "561")
                } elseif ($model -eq "cond_paired_wb") {
                    $splat_args2=@("-eg", "Exposure=W_O1", "-eg", "Exposure=W_O2", "-eg", "Cohort_of_xpt=old", "-eg", "Condition=Immuno")
                    $splat_args3=@("-gk", "ExpCond", "-md", "ExpCond")
                    $channel_splat_args=@("-c", "561")
                } elseif ($model -eq "core") {
                    $splat_args2=@("-eg", "Condition=backward", "-eg", "Condition=backward_24", "-eg", "Condition=shock_24", "-eg", "Condition=control_24", "-eg", "Brain_size=Largee")
                    $splat_args3=@("-gk", "Condition", "-md", "Condition")
                } elseif ($model -eq "core_sex") {
                    $splat_args2=@("-eg", "Condition=backward", "-eg", "Condition=backward_24", "-eg", "Condition=shock_24", "-eg", "Condition=control_24", "-eg", "Brain_size=Largee")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-md", "Condition + Condition*Sex")
                } elseif ($model -eq "all") {
                    $splat_args2=@("-eg", "Brain_size=Largee")
                    $splat_args3=@("-gk", "Condition", "-md", "Condition")
                } elseif ($model -eq "all_sex") {
                    $splat_args2=@("-eg", "Brain_size=Largee")
                    $splat_args3=@("-gk", "Condition", "-gk", "Sex", "-md", "Condition + Condition*Sex")
                } elseif ($model -eq "core_merged") {
                    $splat_args2=@("-eg", "ConditionMerged=backward", "-eg", "Brain_size=Largee")
                    $splat_args3=@("-gk", "ConditionMerged", "-md", "ConditionMerged")
                } elseif ($model -eq "core_merged_sex") {
                    $splat_args2=@("-eg", "ConditionMerged=backward", "-eg", "Brain_size=Largee")
                    $splat_args3=@("-gk", "ConditionMerged", "-gk", "Sex", "-md", "ConditionMerged + ConditionMerged*Sex")
                } elseif ($model -eq "all_merged") {
                    $splat_args2=@("-eg", "Brain_size=Largee")
                    $splat_args3=@("-gk", "ConditionMerged", "-md", "ConditionMerged")
                } elseif ($model -eq "all_merged_sex") {
                    $splat_args2=@("-eg", "Brain_size=Largee")
                    $splat_args3=@("-gk", "ConditionMerged", "-gk", "Sex", "-md", "ConditionMerged + ConditionMerged*Sex")
                } else {
                    exit
                }

                $measures=@("cell_count")#, "cell_density", "total_volume", "total_volume_mm3")
                foreach ($measure in $measures) {
                    if ($fmt -eq "vaa3d") {
                        $splat_args1=@("-oam", "D:\$measure`_deseq2_simple_$atlas\$model\metadata_$model`.csv", "--output-all-stat-deseq2-valid-leafs-only")
                    } else {
                        $splat_args1=@("--output-all-stat-deseq2-valid-leafs-only")
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

                    cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Cynthia.csv" -r "E:\imaging\fused" @channel_splat_args -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*$atlas`.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "$measure_key" -ir 8 -er 512 --vaa3d-atlas-path "D:\$atlas`, CPL v1 vaa3d atlas.csv" --output-path-all-stat-deseq2 "D:\$measure`_deseq2_simple_$atlas\$model\$output_fmt\data-$measure`-$output_fmt`-$threshold`_threshold-{channel}-$model`.csv" --output-path-all-norm-deseq2 "D:\$measure`_deseq2_simple_$atlas\$model\$output_fmt\norm-$measure`-$output_fmt`-$threshold`_threshold-{channel}-$model`.csv" @splat_args1 @splat_args2 @splat_args3 @splat_args4 -vk "Total volume (mm3)"

                    cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Cynthia.csv" -r "E:\imaging\fused" @channel_splat_args -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*$atlas`.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "$measure_key" -ir 8 -er 512 --vaa3d-atlas-path "D:\$atlas`, CPL v1 vaa3d atlas.csv" -oa "D:\$measure`_analyzed_$atlas\$model\$fmt\deseq2\$measure`_deseq2_$fmt`-$threshold`_threshold-{channel}-$model`.csv" @splat_args2 @splat_args3 @splat_args4 --deseq2 @splat_args5 -vk "Total volume (mm3)"

                    #cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" @channel_splat_args -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*$atlas`.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "$measure_key" -ir 8 -er 512 --vaa3d-atlas-path "D:\$atlas`, CPL v1 vaa3d atlas.csv" -oa "D:\$measure`_analyzed_$atlas\$model\$fmt\ttest\$measure`_ttest_$fmt`-$threshold`_threshold-{channel}-$model`.csv" @splat_args2 @splat_args3 @splat_args4 --t-test --t-test-false-discovery -vk "Total volume (mm3)"
                }

                #cell_3d_analyze_summaries -m "D:\Mega data for lightsheet mice - Yidan.csv" -r "E:\imaging\fused" @channel_splat_args -sk "Image_name" -g "**/region_summary_$fmt`_$threshold`_threshold*{name}*{channel}*cell_data*$atlas`.csv" -ak "Total voxels" -ak "Total volume (mm3)" -ak "Total cells" -ak "Cell density (mm3)" -dk "Cell density (mm3)" -ir 8 -er 512 --vaa3d-atlas-path "D:\$atlas`, CPL v1 vaa3d atlas.csv" -or "D:\collated_wide_$atlas\$model\$fmt\wide_raw_$fmt`-$threshold`_threshold-{channel}-$model`.csv" -og "D:\collated_wide_$atlas\$model\$fmt\wide_averages_$fmt`-$threshold`_threshold-{channel}-$model`.csv" @splat_args2 @splat_args3 -vk "Total volume (mm3)"
            }
        }
    }
}
