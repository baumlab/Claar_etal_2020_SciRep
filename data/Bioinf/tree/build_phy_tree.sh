#working from KI_Platy
#Originally from Moorea_Sym pipeline

#data/Bioinf/clust/all_rep_set_rep_set.fasta = seqs
#data/tax_table.txt = ISDs
#column 9 is clade ID
#awk '{print $1}' data/tax_table.txt = denovo name
#awk '{print $8}' data/tax_table.txt = clade name

#subset ID lists by clade
awk '$9 ~ /^A/{ print $1; }' data/tax_table.txt > data/Bioinf/tree/A_tree_seq.ids
awk '$9 ~ /^B/{ print $1; }' data/tax_table.txt > data/Bioinf/tree/B_tree_seq.ids
awk '$9 ~ /^C/{ print $1; }' data/tax_table.txt > data/Bioinf/tree/C_tree_seq.ids
awk '$9 ~ /^D/{ print $1; }' data/tax_table.txt > data/Bioinf/tree/D_tree_seq.ids
awk '$9 ~ /^E/{ print $1; }' data/tax_table.txt > data/Bioinf/tree/E_tree_seq.ids
awk '$9 ~ /^F/{ print $1; }' data/tax_table.txt > data/Bioinf/tree/F_tree_seq.ids
awk '$9 ~ /^G/{ print $1; }' data/tax_table.txt > data/Bioinf/tree/G_tree_seq.ids
awk '$9 ~ /^H/{ print $1; }' data/tax_table.txt > data/Bioinf/tree/H_tree_seq.ids
awk '$9 ~ /^I/{ print $1; }' data/tax_table.txt > data/Bioinf/tree/I_tree_seq.ids

#filter seqs by id list
filter_fasta.py -f data/Bioinf/clust/all_rep_set_rep_set.fasta -s data/Bioinf/tree/A_tree_seq.ids -o data/Bioinf/tree/A_tree_seqs.fasta
filter_fasta.py -f data/Bioinf/clust/all_rep_set_rep_set.fasta -s data/Bioinf/tree/B_tree_seq.ids -o data/Bioinf/tree/B_tree_seqs.fasta
filter_fasta.py -f data/Bioinf/clust/all_rep_set_rep_set.fasta -s data/Bioinf/tree/C_tree_seq.ids -o data/Bioinf/tree/C_tree_seqs.fasta
filter_fasta.py -f data/Bioinf/clust/all_rep_set_rep_set.fasta -s data/Bioinf/tree/D_tree_seq.ids -o data/Bioinf/tree/D_tree_seqs.fasta
filter_fasta.py -f data/Bioinf/clust/all_rep_set_rep_set.fasta -s data/Bioinf/tree/E_tree_seq.ids -o data/Bioinf/tree/E_tree_seqs.fasta
filter_fasta.py -f data/Bioinf/clust/all_rep_set_rep_set.fasta -s data/Bioinf/tree/F_tree_seq.ids -o data/Bioinf/tree/F_tree_seqs.fasta
filter_fasta.py -f data/Bioinf/clust/all_rep_set_rep_set.fasta -s data/Bioinf/tree/G_tree_seq.ids -o data/Bioinf/tree/G_tree_seqs.fasta
filter_fasta.py -f data/Bioinf/clust/all_rep_set_rep_set.fasta -s data/Bioinf/tree/H_tree_seq.ids -o data/Bioinf/tree/H_tree_seqs.fasta
filter_fasta.py -f data/Bioinf/clust/all_rep_set_rep_set.fasta -s data/Bioinf/tree/I_tree_seq.ids -o data/Bioinf/tree/I_tree_seqs.fasta

#check seq number
grep -c ">" data/Bioinf/tree/A_tree_seqs.fasta
grep -c ">" data/Bioinf/tree/B_tree_seqs.fasta
grep -c ">" data/Bioinf/tree/C_tree_seqs.fasta
grep -c ">" data/Bioinf/tree/D_tree_seqs.fasta
grep -c ">" data/Bioinf/tree/E_tree_seqs.fasta
grep -c ">" data/Bioinf/tree/F_tree_seqs.fasta
grep -c ">" data/Bioinf/tree/G_tree_seqs.fasta
grep -c ">" data/Bioinf/tree/H_tree_seqs.fasta
grep -c ">" data/Bioinf/tree/I_tree_seqs.fasta

#align fasta files for each clade
align_seqs.py -i data/Bioinf/tree/A_tree_seqs.fasta -m muscle -o data/Bioinf/tree/
# align_seqs.py -i data/Bioinf/tree/B_tree_seqs.fasta -m muscle -o data/Bioinf/tree/
align_seqs.py -i data/Bioinf/tree/C_tree_seqs.fasta -m muscle -o data/Bioinf/tree/
align_seqs.py -i data/Bioinf/tree/D_tree_seqs.fasta -m muscle -o data/Bioinf/tree/
# align_seqs.py -i data/Bioinf/tree/E_tree_seqs.fasta -m muscle -o data/Bioinf/tree/
# align_seqs.py -i data/Bioinf/tree/F_tree_seqs.fasta -m muscle -o data/Bioinf/tree/
align_seqs.py -i data/Bioinf/tree/G_tree_seqs.fasta -m muscle -o data/Bioinf/tree/
# align_seqs.py -i data/Bioinf/tree/H_tree_seqs.fasta -m muscle -o data/Bioinf/tree/
# align_seqs.py -i data/Bioinf/tree/I_tree_seqs.fasta -m muscle -o data/Bioinf/tree/

#remove extra header info from alignments
sed 's/ .*//' data/Bioinf/tree/A_tree_seqs_aligned.fasta > data/Bioinf/tree/A_tree_seqs_aligned_clean.fasta
# sed 's/ .*//' data/Bioinf/tree/B_tree_seqs_aligned.fasta > data/Bioinf/tree/B_tree_seqs_aligned_clean.fasta
sed 's/ .*//' data/Bioinf/tree/C_tree_seqs_aligned.fasta > data/Bioinf/tree/C_tree_seqs_aligned_clean.fasta
sed 's/ .*//' data/Bioinf/tree/D_tree_seqs_aligned.fasta > data/Bioinf/tree/D_tree_seqs_aligned_clean.fasta
# sed 's/ .*//' data/Bioinf/tree/E_tree_seqs_aligned.fasta > data/Bioinf/tree/E_tree_seqs_aligned_clean.fasta
# sed 's/ .*//' data/Bioinf/tree/F_tree_seqs_aligned.fasta > data/Bioinf/tree/F_tree_seqs_aligned_clean.fasta
sed 's/ .*//' data/Bioinf/tree/G_tree_seqs_aligned.fasta > data/Bioinf/tree/G_tree_seqs_aligned_clean.fasta
# sed 's/ .*//' data/Bioinf/tree/H_tree_seqs_aligned.fasta > data/Bioinf/tree/H_tree_seqs_aligned_clean.fasta
# sed 's/ .*//' data/Bioinf/tree/I_tree_seqs_aligned.fasta > data/Bioinf/tree/I_tree_seqs_aligned_clean.fasta
