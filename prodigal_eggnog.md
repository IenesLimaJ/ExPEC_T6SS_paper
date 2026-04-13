### Running Prodigal and EggNOG using the reference sequences as input 
Analysis performed on UGA-GACRC cluster 

1) Use the nucleotide sequence of each reference as input in Prodigal, to obtain the aminoacid sequence and GFF 

```
prodigal -i "$FASTA1" -a "$OUTDIR/$(basename "$FASTA1" .fasta)_proteins.faa" -d "$OUTDIR/$(basename "$FASTA1" .fasta)_genes.fna" -o "$OUTDIR/$(basename "$FASTA1" .fasta)_prodigal.gff" -p meta
```

NOTE: The GFF file labeled all regions as CDS and did not identied any gene. So, the aminoacid sequences was used as input for EggNog

2) EggNOG was applied in an attempt to identify the T6SS genes 

```
emapper.py --override -i ${INPUT1} --itype proteins -m diamond --data_dir ${DATABASE} -o t6ss_1 --output_dir ${OUTPUT1} --cpu 12
```

