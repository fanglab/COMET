COMET v1.1.0-brain comet dmr output guide
metilene | target 'AD_Microglia' (2 samples) vs background 'MCI_NCI_Microglia' (4 samples)
Parameters and tool version: 02.run_info/run_parameters.txt; inputs: 02.run_info/input_files.tsv

01.dmr/
  AD_Microglia_vs_MCI_NCI_Microglia.metilene.all.tsv
      every region reported by metilene, one row each:
        chrom, start, end            BED coordinates (0-based start, end exclusive)
        n_cpg                        CpGs in the region as counted by the tool; used for
                                     --min-cpg and the depth filter
        n_cpg_any_sample             CpGs in the region covered in at least one sample
        mean_AD_Microglia, mean_MCI_NCI_Microglia       mean methylation (fraction) of each group
        diff_AD_Microglia_minus_MCI_NCI_Microglia     mean_AD_Microglia - mean_MCI_NCI_Microglia
        p_value, q_value, statistic  metilene MWU p-value, BH q-value, 2D-KS p-value
        tool_significant             1 if metilene BH q-value < 0.05 (metilene also requires |diff| >= 0.2)
        min_mean_depth               lowest per-sample mean depth (summed n_total / n_cpg)
        min_individual_diff          smallest gap between any target and any background sample,
                                     in the direction of the group difference (sample level =
                                     summed n_mod / n_total; NA if a sample has no coverage)
        pass_dmr_diff                1 if |diff| >= 0.2 (--dmr-diff)
        pass_individual_diff         1 if min_individual_diff >= 0.2
        pass_min_cpg                 1 if n_cpg >= 5
        pass_min_depth               1 if every sample has mean depth >= 3
        significant                  1 if tool_significant, pass_dmr_diff, pass_min_cpg and
                                     pass_min_depth are 1, and pass_individual_diff
                                     (--dmr-level individual)
        mean_<sample>                methylation level of each input file over the region
                                     (summed n_mod / summed n_total; NA if not covered)
        meth_<sample>, total_<sample>  summed n_mod and n_total of each sample over the region
  AD_Microglia_vs_MCI_NCI_Microglia.metilene.sig.tsv     rows of all.tsv with significant = 1
  AD_Microglia_vs_MCI_NCI_Microglia.metilene.sig.bed     the same regions as BED; name = AD_Microglia:<difference>
  AD_Microglia_vs_MCI_NCI_Microglia.metilene.raw.tsv     unmodified metilene output

02.run_info/
  OUTPUT_README.txt      this guide
  run_parameters.txt     command line, parameters, tool versions
  input_files.tsv        sample name, group and path of every input file
  logs/                  tool logs and summary.txt (region counts)
