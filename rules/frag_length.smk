# ------------------------------ #
#  Beginning of Snakemake Rules  #
# ------------------------------ #

path_to_data = config['out_dir'] 

# TO FIX
# rule collect_hs_metrics:
#     input:
#         bam="path/to/HUCON_Capture_Pool_7_S236_L002_HPVUnmapped_Hg19_Aligned.sorted.bam",
#         ref="/cluster/tools/data/genomes/human/hg19/bwa/ucsc.hg19.fasta",
#         bait_intervals="/cluster/projects/scottgroup/people/jinfeng/data/bait_target_bed/HPV/Human4_STK11onHuman3.interval_list",
#         target_intervals="/cluster/projects/scottgroup/people/jinfeng/data/bait_target_bed/HPV/Human4_STK11onHuman3.interval_list"
#     output:
#         metrics_txt="path/to/HsMetrics/HUCON_Capture_Pool_7_S236_L002_hsMetrics.txt"
#     params:
#         picard_dir="path/to/picard_dir"
#     shell:
#         """
#         java -jar {params.picard_dir}/picard.jar CollectHsMetrics \
#         I={input.bam} O={output.metrics_txt} R={input.ref} \
#         BAIT_INTERVALS={input.bait_intervals} TARGET_INTERVALS={input.target_intervals}
#         """



rule collect_insert_size_metrics:
    input:
        # bam=lambda wildcards: get_sample_path("OPC_test", wildcards.sample)
        bam = path_to_data + "fragment_metrics/"+ "virusDuplex/fusion_aligned/" + "{sample}_fusion_aligned_sorted.bam"
    output:
        metrics_txt= path_to_data + "fragment_metrics/"+ "fragment_length/" + "{sample}_metrics.txt",
        histogram_pdf=path_to_data + "fragment_metrics/"+ "fragment_length/" + "{sample}_histogram_plot.pdf"
    params:
        picard_dir="path/to/picard_dir"
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        """
        java -jar /cluster/tools/software/picard/2.10.9/picard.jar CollectInsertSizeMetrics \
        I={input.bam} O={output.metrics_txt} H={output.histogram_pdf} M=0.5
        """


        
rule create_metrics_summary:
    input:
        metrics_files=
        lambda wildcards: expand(
        path_to_data + "fragment_metrics/"+ "fragment_length/" + "{sample}_metrics.txt",
        sample = get_all_samplesnames_list(None)
        )
    output:
        outfile=path_to_data + "fragment_metrics/"+ "fragment_length/" +"fragment_length_summary.csv"
    resources: cpus=1, mem_mb=15000, time_min=60
    shell:
        """
        # Create the header
        outfile={output.outfile}
        echo -n "filename,READ_PAIRS," > $outfile
        for i in $(seq 0 600); do
            echo -n "$i," >> $outfile
        done
        echo "" >> $outfile

        # Process each metrics file
        for file in {input.metrics_files}; do
            # Create an array with 601 zeroes
            arr=($(for i in $(seq 0 600); do echo -n "0 "; done))

            # Extract the filename
            filename=$(basename $file)

            # Extract the READ_PAIRS value
            read_pairs=$(awk '/READ_PAIRS/{{getline; print $7}}' $file)

            # Update the array with the actual data
            awk -v arr="${{arr[*]}}" -v filename="$filename" -v read_pairs="$read_pairs" \\
                'BEGIN{{split(arr, a, " ")}}
                $1 ~ /^[0-9]+$/ && $1 >= 0 && $1 <= 600 {{a[$1+1]=$2}}
                END{{printf("%s,%s,", filename, read_pairs); for(i=1; i<=601; i++) printf("%s,", a[i]); printf("\\n")}}' $file >> $outfile
        done
        """

