# ROOM FOR IMPROVEMENT
# as i am unclear on what other elements one might want to annotated
# using this method, i've left it with a tyr specific annotation but
# included some structure to the yaml to allow for more flexibility
# FURTHER ROOM
# this is a simple version that expects the filt_intersect to exist but
# the docs make it clear that is not a guarantee. a checkpoint (instead
# of a rule) would be better here and will implement if time allows
rule sf_annotate_tyr:
    '''
    what do i do?
    '''
    input:
        fnas = rules.prep_control_files.output.fnas,
        gffs = rules.prep_control_files.output.gffs
    output:
        gff    = 'results/gene_finder/macpha6_tyr.filt_intersect.gff',
        sf_ids = 'results/gene_finder/macpha6_tyr.filt_intersect.ids'
    params:
        db_dir  = config['star_db'],
        prefix  = config['annotate']['tyr']['prefix'],
        out_dir = lambda wildcards, output: os.path.dirname(output.gff),
    threads: 2
    resources:
        time = 20,
        mem_mb = 2000
    shell:
        '''
			starfish annotate \
				-T {threads} \
				-x {params.prefix} \
				-a {input.fnas} \
				-g {input.gffs} \
				-p {params.db_dir}/YRsuperfams.p1-512.hmm \
				-P {params.db_dir}/YRsuperfamRefs.faa \
				-i tyr \
				-o {params.out_dir}
        '''

rule sad_old_gff:
    output:
        'results/macpha6.gff3'
    params:
        gffs = config['gff_dir']
    shell:
        '''
            cat {params.gffs}* > {output}
        '''

rule sf_consolidate:
    '''
    what about me?
    '''
    input:
        old_gff = rules.sad_old_gff.output,
        new_gff = rules.sf_annotate_tyr.output
    output:
        gff   = 'results/gene_finder/macpha6_tyr.filt_intersect.consolidated.gff',
        fpath = 'results/controls/ome2consolidatedGFF.txt',
    params:
        out_dir = lambda wildcards, output: os.path.dirname(output.gff)
    threads: 1
    resources:
        time = 10,
        mem_mb = 2000
    shell:
        '''
            starfish consolidate \
                -o {params.out_dir} \
                -g {input.old_gff} \
                -G {input.new_gff}

            realpath {output.gff} | perl -pe 's/^/macpha6\\t/' > {output.fpath}
        '''

# ROOM FOR IMPROVEMENT
# items in this rule (and others) like -m/--merge could be set in the config
# as opposed to in the rule
# FURTHER
# this is similar to annotate above where its specifically sketching tyr but
# could/should be more generic to annotate other elements in parallel
rule sf_sketch_tyr:
    '''
    same thing
    '''
    input:
        sf_ids = rules.sf_annotate_tyr.output.sf_ids,
        fpath  = rules.sf_consolidate.output.fpath
    output:
        bed     = 'results/gene_finder/macpha6.bed',
        mat     = 'results/gene_finder/macpha6.mat',
        tyr_bed = 'results/gene_finder/macpha6.tyr.bed',
    params:
        merge   = '10000',
        prefix  = 'macpha6',
        out_dir = lambda wildcards, output: os.path.dirname(output.bed),
    threads: 1
    resources:
        time = 10,
        mem_mb = 2000
    shell:
        '''
            starfish sketch \
                -m {params.merge} \
                -q {input.sf_ids} \
                -g {input.fpath} \
                -i s \
                -x {params.prefix} \
                -o {params.out_dir}

            grep -P '\\ttyr\\t' {output.bed} > {output.tyr_bed}
        '''

