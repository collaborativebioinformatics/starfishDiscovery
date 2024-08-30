
# ROOM FOR IMPROVEMENT
# not flexible and requires specific extensions fasta and gff3
localrules: prep_path_files
rule prep_control_files:
    input:
        fna = config['fna_dir'],
        gff = config['gff_dir']
    output:
        fnas = 'results/controls/ome2assembly.txt',
        gffs = 'results/controls/ome2gff.txt'
    run:
        file_paths(input.fna, output.fnas, 'fasta')
        file_paths(input.gff, output.gffs, 'gff3')

rule build_blastdb:
    input:
        rules.prep_control_files.output.fnas,
    output:
        'results/blastdb/macpha6.assemblies.not',
        'results/blastdb/macpha6.assemblies.nto',
        'results/blastdb/macpha6.assemblies.nin',
        'results/blastdb/macpha6.assemblies.ndb',
        'results/blastdb/macpha6.assemblies.nog',
        'results/blastdb/macpha6.assemblies.nsq',
        'results/blastdb/macpha6.assemblies.nos',
        'results/blastdb/macpha6.assemblies.ntf',
        'results/blastdb/macpha6.assemblies.nhr',
        fna = 'results/blastdb/macpha6.assemblies.fna',
    params:
        out_dir = lambda wildcards, output: os.path.splitext(output.fna)[0]
    threads: 1
    resources:
        time = 10,
        mem_mb = 2000
    shell:
        '''
            cut -f2 {input} | xargs cat > {output.fna}

            makeblastdb \
                -in {output.fna} \
                -out {params.out_dir} \
                -parse_seqids \
                -dbtype nucl
        '''
# FIX ME
# seq-gc.sh and all other aux commands from starfish should be soruced
# from the container. temporarily using a local
rule calc_gc_content:
    input:
        rules.build_blastdb.output.fna
    output:
        'results/viz/macpha6.assemblies.gcContent_w1000.bed'
    params:
        aux_dir = config['star_aux']
    threads: 1
    resources:
        time = 10,
        mem_mb = 2000
    shell:
        '''
            {params.aux_dir}/seq-gc.sh \
                -Nbw 1000 \
                {input} \
                > {output}
        '''

# ROOM FOR IMPROVEMENT
# not flexible and requires specific extensions fasta and gff3
rule parse_eggnogs:
    input:
        config['ann_dir']
    output:
        'results/ann/macph6.gene2emap.txt'
    threads: 1
    resources:
        time = 10,
        mem_mb = 2000
    shell:
        '''
            cut -f1,12 {input}/*emapper.annotations \
                | grep -v  '#' \
                | grep -v -P '\t-' \
                | perl -pe 's/\t/\tEMAP\t/' \
                | grep -vP '\tNA' \
                > {output}
        '''

rule get_narrow_ortho:
    input:
        config['ann_dir']
    output:
        'results/ann/macph6.gene2og.txt'
    threads: 1
    resources:
        time = 10,
        mem_mb = 2000
    shell:
        '''
            cut -f1,10 {input}/*emapper.annotations \
                | grep -v '#' \
                | perl -pe 's/^([^\s]+?)\\t([^\|]+).+$/\\1\\t\\2/' \
                > {output}
        '''

rule convert_to_mcl:
    input:
        rules.get_narrow_ortho.output
    output:
        'results/ann/macph6.gene2og.mcl'
    params:
        aux_dir = config['star_aux'],
        out_dir = lambda wildcards, output: os.path.dirname(output[0])
    threads: 1
    resources:
        time = 10,
        mem_mb = 2000
    shell:
        '''
            {params.aux_dir}/geneOG2mclFormat.pl \
                -i {input} \
                -o {params.out_dir}
        '''

