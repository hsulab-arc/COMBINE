#!/bin/bash

# making tmp directory for outputs
mkdir -p tmp

# simplex basecall with dorado
dorado basecaller ~/dorado_models/dna_r10.4.1_e8.2_260bps_sup@v4.1.0 -r pod5s/ --emit-moves > tmp/unmapped_reads_with_moves.sam

# call duplex pairs from unmapped_reads_with_moves.sam
duplex_tools pair --output_dir tmp/pairs_from_bam tmp/unmapped_reads_with_moves.sam

# finding additional pairs from non-split pairs and splitting them
duplex_tools split_pairs tmp/unmapped_reads_with_moves.sam pod5s/ tmp/pod5s_splitduplex/

# adding all split pairs to file for duplex calling
cat tmp/pod5s_splitduplex/*_pair_ids.txt > tmp/split_duplex_pair_ids.txt

# adding all split pairs to file to filter simplex fastq
cat tmp/pod5s_splitduplex/*_split_duplex.txt > tmp/split_duplex.txt

# basecalling original duplex reads
dorado duplex ~/dorado_models/dna_r10.4.1_e8.2_260bps_sup@v4.1.0 -r pod5s/ --pairs tmp/pairs_from_bam/pair_ids_filtered.txt > tmp/duplex_orig.sam

# basecalling split duplex reads
dorado duplex ~/dorado_models/dna_r10.4.1_e8.2_260bps_sup@v4.1.0 -r tmp/pod5s_splitduplex/ --pairs tmp/split_duplex_pair_ids.txt > tmp/duplex_splitduplex.sam

echo "Converting sams to fastqs..."

# convert simplex sam to fastq
samtools fastq tmp/unmapped_reads_with_moves.sam > tmp/simplex.fastq

# convert original duplex to fastq
samtools fastq tmp/duplex_orig.sam > tmp/duplex_orig.fastq

# convert split duplex to fastq
samtools fastq tmp/duplex_splitduplex.sam > tmp/duplex_splitduplex.fastq

echo "Filtering fastqs..."

# saving all sequence IDs from duplex_orig.fastq to file for removal from simplex.fastq
awk 'NR%4==1{gsub("@",""); gsub(";","\n",$1); print $1}' tmp/duplex_orig.fastq > tmp/duplex_orig_IDs.txt

# saving all sequence IDs from duplex_splitduplex.fastq to file for lookup in split_duplex.txt
awk 'NR%4==1{gsub("@",""); gsub(";","\t",$1); print $1}' tmp/duplex_splitduplex.fastq > tmp/duplex_splitduplex_IDs.txt

# looking up the nonsplit read ID corresponding to each read in splitduplex_IDs
awk 'BEGIN{FS="\t"} NR==FNR{a[$1,$2]; next} ($2,$3) in a {print $1}' tmp/duplex_splitduplex_IDs.txt tmp/split_duplex.txt > tmp/duplex_nonsplit_IDs.txt

# concatenating the read IDs from removal
cat tmp/duplex_orig_IDs.txt tmp/duplex_nonsplit_IDs.txt > tmp/remove.txt

# saving all sequence IDs from simplex.fastq to file
awk 'NR%4==1{gsub("@",""); gsub(";","\n",$1); print $1}' tmp/simplex.fastq > tmp/simplex_IDs.txt

# getting complement of sequence IDs from simplex_IDs.txt
grep -v -f tmp/remove.txt tmp/simplex_IDs.txt > tmp/keep.txt

# filtering simplex.fastq to contain only reads that were not used in duplex calls
seqtk subseq tmp/simplex.fastq tmp/keep.txt > tmp/simplex_filt.fastq
mv -f tmp/simplex_filt.fastq tmp/simplex.fastq

# concatenating both duplex fastqs
cat tmp/duplex_orig.fastq tmp/duplex_splitduplex.fastq > tmp/duplex.fastq

# filtering for reads greater than given length
seqtk seq -L 500 tmp/simplex.fastq > tmp/simplex_filt.fastq
seqtk seq -L 500 tmp/duplex.fastq > tmp/duplex_filt.fastq
mv -f tmp/simplex_filt.fastq tmp/simplex.fastq
mv -f tmp/duplex_filt.fastq tmp/duplex.fastq

echo "Writing fastqs..."

# splitting filtered fastq into smaller files
mkdir -p fastqs && split -l 400000 tmp/simplex.fastq fastqs/fastq_simplex_ --numeric-suffixes=1 --suffix-length=3 && for f in fastqs/fastq_simplex_*; do mv "$f" "$f.fastq"; done
split -l 400000 tmp/duplex.fastq fastqs/fastq_duplex_ --numeric-suffixes=1 --suffix-length=3 && for f in fastqs/fastq_duplex_*; do mv "$f" "$f.fastq"; done

echo "Compressing fastqs..."

# compressing all fastqs
find ./fastqs -name '*.fastq' | parallel -j $(nproc) 'gzip {}'

# removing temporary directory
# rm -rf tmp

echo "Done."