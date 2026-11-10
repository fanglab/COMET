COMET v1.1.0 output guide
sample brain_03 | platform pacbio | tissue brain | cell types Exc_neuron,Inh_neuron,Astrocyte,Microglia,Oligodendrocyte,OPC
Files are prefixed by the sample name. Parameters and tool versions: 03.run_info/run_parameters.txt

01.cell_type_assignment/  read-level cell-type assignment
  cell_type_prob/brain_03.cell_type_prob.tsv
      one row per scored read (>= 5 informative CpGs in the reference panel); columns:
      read_id, n_cpg (informative CpGs), label (panel column 1..6, 0 = unassigned),
      cell_type, P_<cell type> (posterior of each cell type), rescued_lineage ('.' if none).
      A read is assigned when its top posterior >= 0.65 and >= 2x the runner-up.
  read_list/cell_type/brain_03.<cell type>.reads.txt
      IDs of the reads assigned to each cell type.
  read_list/brain_lineage_rescued/brain_03.<lineage>.reads.txt
      unassigned reads rescued at lineage level (Neuron=1,2, Macroglia=3,5,6; summed posterior
      >= 0.65 and every member > best outside type + 0.10).
      Neuron = Exc_neuron + Inh_neuron; Macroglia = Astrocyte + Oligodendrocyte + OPC.
  BAM/cell_type/brain_03.<cell type>.bam (.csi)
  BAM/brain_lineage_rescued/brain_03.<lineage>.bam (.csi)
      the reads of each set taken from the input BAM with all tags (MM/ML) kept, indexed;
      usable like a sorted-cell library by any methylation tool.
  BED/bulk/        all reads (bulk reference)
  BED/cell_type/   each cell type
  BED/brain_lineage_rescued/   each rescued lineage
      per-site calls, bgzipped + tabix-indexed, columns: chrom, start, end, n_mod, n_total
      (0-based start; methylation level = n_mod / n_total)
      *.CG.bed.gz    5mC at CpGs from COMET's read-level calls; a call counts as modified at
                     ML >= tau and unmodified at ML <= 255 - tau (tau = pass-number cutoff);
                     both strands are merged onto the C of the CpG

02.cell_type_quantification/  cell-type quantification
  brain_03.composition.tsv
      group, n_reads and three proportions:
        frac_single          among reads assigned to a single cell type
        frac_single_lineage  among single-type + lineage-rescued reads
        frac_all             among all scored reads, including unassigned
  brain_03.gw_methylation.tsv
      genome-wide CpG methylation of each read set from all CpG calls of its reads
      (in or outside the reference panel): group, n_reads, n_meth_calls, n_unmeth_calls, methylation.
      Bulk = every read with CpG calls; Unassigned = scored reads not assigned or rescued.

03.run_info/
  OUTPUT_README.txt     this guide
  run_parameters.txt    command line, parameters, reference panel, tool versions
  brain_03.pass_number_table.tsv
      ML cutoff per CCS pass number: np, tau0, n_meth(tau0), n_unmeth(tau0), tau_np,
      n_meth(tau_np), n_unmeth(tau_np); tau_np = 256 means calls from that pass number are dropped
  logs/                 per-step logs

Files in this run:
  01.cell_type_assignment/BAM/brain_lineage_rescued/brain_03.Macroglia.bam
  01.cell_type_assignment/BAM/brain_lineage_rescued/brain_03.Neuron.bam
  01.cell_type_assignment/BAM/cell_type/brain_03.Astrocyte.bam
  01.cell_type_assignment/BAM/cell_type/brain_03.Exc_neuron.bam
  01.cell_type_assignment/BAM/cell_type/brain_03.Inh_neuron.bam
  01.cell_type_assignment/BAM/cell_type/brain_03.Microglia.bam
  01.cell_type_assignment/BAM/cell_type/brain_03.Oligodendrocyte.bam
  01.cell_type_assignment/BAM/cell_type/brain_03.OPC.bam
  01.cell_type_assignment/BED/brain_lineage_rescued/brain_03.Macroglia.CG.bed.gz
  01.cell_type_assignment/BED/brain_lineage_rescued/brain_03.Neuron.CG.bed.gz
  01.cell_type_assignment/BED/bulk/brain_03.Bulk.CG.bed.gz
  01.cell_type_assignment/BED/cell_type/brain_03.Astrocyte.CG.bed.gz
  01.cell_type_assignment/BED/cell_type/brain_03.Exc_neuron.CG.bed.gz
  01.cell_type_assignment/BED/cell_type/brain_03.Inh_neuron.CG.bed.gz
  01.cell_type_assignment/BED/cell_type/brain_03.Microglia.CG.bed.gz
  01.cell_type_assignment/BED/cell_type/brain_03.Oligodendrocyte.CG.bed.gz
  01.cell_type_assignment/BED/cell_type/brain_03.OPC.CG.bed.gz
  01.cell_type_assignment/cell_type_prob/brain_03.cell_type_prob.tsv
  01.cell_type_assignment/read_list/brain_lineage_rescued/brain_03.Macroglia.reads.txt
  01.cell_type_assignment/read_list/brain_lineage_rescued/brain_03.Neuron.reads.txt
  01.cell_type_assignment/read_list/cell_type/brain_03.Astrocyte.reads.txt
  01.cell_type_assignment/read_list/cell_type/brain_03.Exc_neuron.reads.txt
  01.cell_type_assignment/read_list/cell_type/brain_03.Inh_neuron.reads.txt
  01.cell_type_assignment/read_list/cell_type/brain_03.Microglia.reads.txt
  01.cell_type_assignment/read_list/cell_type/brain_03.Oligodendrocyte.reads.txt
  01.cell_type_assignment/read_list/cell_type/brain_03.OPC.reads.txt
  02.cell_type_quantification/brain_03.composition.tsv
  02.cell_type_quantification/brain_03.gw_methylation.tsv
  03.run_info/brain_03.pass_number_table.tsv
  03.run_info/OUTPUT_README.txt
  03.run_info/run_parameters.txt
