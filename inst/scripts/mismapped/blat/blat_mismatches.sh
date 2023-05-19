
#BSUB -n 44
#BSUB -J align 
#BSUB -o log_%J.out
#BSUB -e log_%J.err
#BSUB -R "select[mem>50] rusage[mem=50] span[hosts=1]"
#BSUB -q rna

pblat \
    ../../dbases/GRCh38.primary_assembly.genome.fa \
    GRCm38.75.fasta \
    -threads=44 \
    -stepSize=5 \
    -repMatch=2253 \
    -minScore=0 \
    -minIdentity=0 \
    -out=pslx \
    -noHead \
    GRCm38.75.pslx

# code from gabay et al to parse blat output...
in_pslx=GRCm38.75.pslx
awk '$1>=68 && $2>0 && $5==0 {aln_blocks_sum=0; split($19,aln_blocks,","); for (i=1;i<=$18;i++) aln_blocks_sum+=aln_blocks[i]; if ($1/aln_blocks_sum >= 0.96) print $0}' $in_pslx \
    | awk 'BEGIN{FS=OFS="\t"; rev["A"]="T";rev["C"]="G";rev["G"]="C";rev["T"]="A"}{aln_length=$1; aln_strand=$9; n=split($10,query_name,"_"); query_chr=query_name[4]; query_strand=query_name[3]; if (query_strand=="-") query_strand_flag=1; else query_strand_flag=0; j=1; for (i=5; i<=n; i++) {split(query_name[i],query_intervals,"-"); query_intervals[1]+=1; 
        for (w=query_intervals[1+(1*query_strand_flag)];
        w!=query_intervals[2-(1*query_strand_flag)]+(1-(2*query_strand_flag));w=w+(1-(2*query_strand_flag))){query_coords[j]=w;j++}}
        n_aln=$18; split($19,aln_blocks,","); split($20,query_blocks,",");
        split($21,target_blocks,","); split($22,query_seqs,",");
        split($23,target_seqs,","); if (aln_strand=="-")
        aln_strand_flag=1; else aln_strand_flag=0; 
        for (cur_block=n_aln-(n_aln+1)*aln_strand_flag;cur_block!=0+(n_aln+1)*aln_strand_flag;
        cur_block=cur_block-1+2*aln_strand_flag){cur_aln_length=aln_blocks[cur_block];
        for (pos=1; pos<=cur_aln_length;pos++)
            {query_nuc=toupper(substr(query_seqs[cur_block],pos,1));
                target_nuc=toupper(substr(target_seqs[cur_block],pos,1));
                if (query_nuc!=target_nuc) {
                    if (aln_strand=="-"){query_nuc=rev[query_nuc];target_nuc=rev[target_nuc]} query_pos=pos+query_blocks[cur_block];
                    query_coord=query_coords[query_pos+(j-1-2*query_pos+1)*aln_strand_flag];
                    target_pos=pos+target_blocks[cur_block]-1; 
                    print query_chr,query_coord-1,query_coord, query_nuc target_nuc,"0", query_strand, $14 "_" target_pos "_" target_pos+1 "_" aln_strand}}}}' \
                      
