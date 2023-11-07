

macs2 callpeak -t CS1.cuttag.bam -c CSIgg.cuttag.bam --outdir out_CS1_vs_IgG_q0.01 --tempdir .  -g hs -f BAMPE -n CS1_vs_IgG --keep-dup all --qval 0.01 -B --SPMR > log.macs2.CS1_vs_IgG.q0.01 2>&1 &


##CEBPB cut&tag

for i in {2..3};do echo 'rep'$i;
   macs2 callpeak -t CS$i.cuttag.bam -c CSIgg.cuttag.bam --outdir out_CS${i}_vs_IgG_q0.01 --tempdir .  -g hs -f BAMPE -n CS${i}_vs_IgG --keep-dup all --qval 0.01 -B --SPMR > log.macs2.CS${i}_vs_IgG.q0.01 2>&1

done > log.run.CEBPB 2>&1 &

##FOSL2
for i in {1..3};do echo 'rep'$i;
   macs2 callpeak -t SF$i.cuttag.bam -c SFIgg.cuttag.bam --outdir out_SF${i}_vs_IgG_q0.01 --tempdir .  -g hs -f BAMPE -n SF${i}_vs_IgG --keep-dup all --qval 0.01 -B --SPMR > log.macs2.SF${i}_vs_IgG.q0.01 2>&1

done > log.run.FOSL2 2>&1 &












#############

macs2 callpeak  -t ../MS1.cuttag.bam  -c ../MSNC.cuttag.bam --outdir out_MS1_q0.05 --tempdir .  -g hs -f BAMPE -n MITF_STB1 --keep-dup all --qval 5e-2 -B --SPMR > log.macs2.MS1.q0.05 2>&1 &

#-f BAMPE --keep-dup all --nomodel --qval 5e-2 -B --SPMR --shift 100 --ext 200 


macs2 callpeak  -t ../MS2.cuttag.bam  -c ../MSNC.cuttag.bam --outdir out_MS2_q0.05 --tempdir .  -g hs -f BAMPE -n MITF_STB2 --keep-dup all --qval 5e-2 -B --SPMR > log.macs2.MS2.q0.05 2>&1 &



macs2 callpeak  -t ../MS3.cuttag.bam  -c ../MSNC.cuttag.bam --outdir out_MS3_q0.05 --tempdir .  -g hs -f BAMPE -n MITF_STB3 --keep-dup all --qval 5e-2 -B --SPMR > log.macs2.MS3.q0.05 2>&1 &















#######










#1 code in snapatac2-core/src/export.rs
"callpeak",
        "-f", "BED",
        "-t", bed_file.as_ref().to_str().unwrap(),
        "--keep-dup", "all",
        "--outdir", format!("{}", dir.path().display()).as_str(),
        "--qvalue", format!("{}", q_value).as_str(),
        "-g", format!("{}", (genome_size as f64 * 0.9).round()).as_str(),
        "--call-summits",
        "--nomodel", "--shift", "-100", "--extsize", "200",
        "--nolambda", #should we turn this off?
        "--tempdir", format!("{}", dir.path().display()).as_str(),
        
genome_size:  data.uns['reference_sequences']['reference_seq_length'].sum()


##2 snapatac call peak parameters
--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR
#  -B --bdg save extended fragment pileup
#  --SPMR SAVE signal per million reads for fragment pileup profiles


#######call macs manually with SnapATAC v2.2.0 (fragmens bed)/SnapATAC v1.0 macs parameter (Tn5 insertions bed)######

ls {1..9}_insertion.bed.gz|sort -h > bed_list.txt

##use SnapATAC v2.2.0 parameters with Tn5 insertions (--nolambda, --call-summits)
#cat bed_list.txt|while read a; do echo $a; prefix=`basename $a .bed.gz`; echo $prefix; mkdir out_${prefix};
#   macs2 callpeak -f BED -t $a --keep-dup all --outdir out_${prefix} --qvalue 0.05 -g 3088269832 --nomode --shift -100 --extsize 200 --nolambda -B --SPMR --call-summits > out_${prefix}.log 2>&1 & #call for insertion bed, in background
#done > log.run.macs2 2>&1

#use SnapATAC v1.0 parameters with Tn5 insertions  (0.05: 5e-2) 
cat bed_list.txt|while read a; do echo $a; prefix=`basename $a .bed.gz`; echo $prefix; mkdir out_${prefix};
   macs2 callpeak -f BED -t $a  --outdir out_${prefix} --tempdir . --qvalue 0.01 -g hs --nomode --shift -100 --extsize 200 -B --SPMR --call-summits > out_${prefix}.log 2>&1 & #call for insertion bed, in background
   #macs2 callpeak -f BED -t $a  --outdir out_${prefix}  --qvalue 5e-2 -g hs --nomode --shift 100 --extsize 200 -B --SPMR > out_${prefix}.log 2>&1 & #call for insertion bed, in background
done > log.run.macs2 2>&1

tail out_*_insertion.log |grep "Done"

##collect 
for i in {1..9}; do echo $i; ln -s out_${i}_insertion/NA_peaks.narrowPeak ./${i}_insertion_peaks.narrowPeak   ;   done
for i in {1..9}; do echo $i; ln -s out_${i}_insertion/NA_treat_pileup.bdg ./${i}_insertion_treat_pileup.bdg   ;   done

##merge to peak.combined.bed by r
##output peak logical table

diff -s <(awk '{print $1":"$2"-"$3}' peaks.combined.bed) <( awk '(NR >1){print $0}'  peaks.combined.logic.txt|cut -f 1 )
#identical


##merge within 100bp then extend to 200 if less than 200
bedtools merge -i summits.combined.bed -d 100 |awk -vOFS='\t' '{if($3-$2<200){ mid= int(0.5*($3+$2));  print $1,mid-100,mid+100}else{ print $0 }}' > summits.combined.merge100.extend200.bed #use this

bedtools merge -i summits.combined.bed -d 500 |awk -vOFS='\t' '{if($3-$2<200){ mid= int(0.5*($3+$2));  print $1,mid-100,mid+100}else{ print $0 }}' > summits.combined.merge500.extend200.bed



