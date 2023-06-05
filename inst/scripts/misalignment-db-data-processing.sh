#BSUB -n 44
#BSUB -J align
#BSUB -o log_%J.out
#BSUB -e log_%J.err
#BSUB -R "select[mem>50] rusage[mem=50] span[hosts=1]"
#BSUB -q rna

# note requires bedtools, perl, and pblat in PATH
# uses ~44 CPU cores for processing

# Approach, awk code, and perl scripts are from https://github.com/a2iEditing/deNovo-Detect
# related to citation
# Gabay, O., Shoshan, Y., Kopel, E. et al. Landscape of adenosine-to-inosine RNA recoding across human tissues. Nat Commun 13, 1184 (2022). https://doi.org/10.1038/s41467-022-28841-4

# RefSeq_Curated_NM_CDS_010817.bed RefSeq_Curated_NM_Exons_010817.bed  RefSeq_Curated_NM_Gene_010817.bed
# were all downloaded from the table browser.

# a gene wide search is too expensive with this number of windows, not
# performing for now. 

# RefSeq_Curated_NM_Gene_010817.bed was modified to match formatting for CDS and Exons using:
# awk 'BEGIN {OFS=FS="\t"} strand=($6 == "+")?"f":"r" {print $1,$2,$3,$4"_gene_0_0_"$1"_"$2 + 1"_"strand, $5, $6}' > RefSeq_Curated_NM_Gene_010817_mod.bed

mkdir -p misalignment_db
cd misalignment_db

if [[ ! -f "hg38_all.fa" ]]; then
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.chroms.tar.gz
  tar -xzvf hg38.analysisSet.chroms.tar.gz
  cat ./hg38.analysisSet.chroms/*.fa > hg38_all.fa
  rm hg38.analysisSet.chroms.tar.gz
  rm -rf hg38.analysisSet.chroms/
fi

get_misalignment_db() {
    in_bed=$1
    fa=$2
    pslx=$3
    out_bed=$4
    threads=$5
    bedtools getfasta -bed $in_bed -fi $fa -s -bedOut 2> /dev/null \
      | awk 'BEGIN{FS=OFS="\t"}{split($4,a,"_"); print $1,$2,$3,a[1]"_"a[2],$5,$6,$7}' \
      | grep -i -v "NN" \
      | perl sequence_splicer_gene_coords.pl /dev/stdin tmp.fa


    perl sequence_chopper_with_coords.pl \
      tmp.fa \
      76 \
      19 \
      tmp_chopped.fa

    # setting minScore and minIdentity to 0 seems unnecessary if we will post filter on match length >= 68
    # changing to defaults of 30 and 90 respectfully
    pblat \
        $fa \
        tmp_chopped.fa \
        -threads=$threads \
        -stepSize=5 \
        -repMatch=2253 \
        -minScore=30 \
        -minIdentity=90 \
        -out=pslx \
        -noHead \
        $pslx


    awk '$1>=68 && $2>0 && $5==0 {aln_blocks_sum=0; split($19,aln_blocks,",");
        for (i=1;i<=$18;i++) aln_blocks_sum+=aln_blocks[i];
          if ($1/aln_blocks_sum >= 0.96) print $0}' $pslx \
        | awk 'BEGIN{FS=OFS="\t";
                     rev["A"]="T";
                     rev["C"]="G";
                     rev["G"]="C";
                     rev["T"]="A"}
               {aln_length=$1;
                aln_strand=$9;
                n=split($10,query_name,"_");
                query_chr=query_name[4];
                query_strand=query_name[3];

                if (query_strand=="-") query_strand_flag=1;
                else query_strand_flag=0;
                j=1;

                for (i=5; i<=n; i++) {split(query_name[i], query_intervals, "-");
                  query_intervals[1]+=1;
                  for (w=query_intervals[1+(1*query_strand_flag)];
                    w!=query_intervals[2-(1*query_strand_flag)]+(1-(2*query_strand_flag));
                    w=w+(1-(2*query_strand_flag))) {query_coords[j]=w;j++}
                }
                n_aln=$18;
                split($19,aln_blocks,",");
                split($20,query_blocks,",");
                split($21,target_blocks,",");
                split($22,query_seqs,",");
                split($23,target_seqs,",");
                if (aln_strand=="-")
                  aln_strand_flag=1;
                else
                  aln_strand_flag=0;

                for (cur_block=n_aln-(n_aln+1)*aln_strand_flag;cur_block!=0+(n_aln+1)*aln_strand_flag;
                     cur_block=cur_block-1+2*aln_strand_flag){cur_aln_length=aln_blocks[cur_block];
                for (pos=1; pos<=cur_aln_length;pos++)
                    {query_nuc=toupper(substr(query_seqs[cur_block],pos,1));
                    target_nuc=toupper(substr(target_seqs[cur_block],pos,1));
                    if (query_nuc!=target_nuc) {
                        if (aln_strand=="-"){query_nuc=rev[query_nuc];target_nuc=rev[target_nuc]} query_pos=pos+query_blocks[cur_block];
                        query_coord=query_coords[query_pos+(j-1-2*query_pos+1)*aln_strand_flag];
                        target_pos=pos+target_blocks[cur_block]-1;
                        print query_chr, query_coord-1, query_coord, query_nuc target_nuc,"0", query_strand, $14 "_" target_pos "_" target_pos+1 "_" aln_strand}}}}' \
        > $out_bed

        rm -f tmp.fa tmp_chopped.fa $pslx
}

get_misalignment_db \
  RefSeq_Curated_NM_CDS_010817.bed \
  hg38_all.fa \
  RefSeq_Curated_NM_CDS_010817_spliced_noNs_chopped_76_19.pslx \
  RefSeq_Curated_NM_CDS_010817_spliced_noNs_chopped_76_19_mms.bed \
  44

get_misalignment_db \
  RefSeq_Curated_NM_Exons_010817.bed \
  hg38_all.fa \
  RefSeq_Curated_NM_Exons_010817_spliced_noNs_chopped_76_19.pslx \
  RefSeq_Curated_NM_Exons_010817_spliced_noNs_chopped_76_19_mms.bed \
  44

