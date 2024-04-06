# Investigation-of-the-Regulatory-Element-of-MSP-genes-in-the-Caenorhabditis-family

Workflow:
1. retrieve upstream 100kb genPromoter.sh
   use the following line to delete sequences that are shorter than 5 bp:
    sed -E '/>/N;/\n[^>].{0,4}$/d' input_file > output_file
3. find motifs: runstreme.sh, optional: use ceqlogo.sh to gen logo
4. find similarity: runTomtom.sh
5. cluster motifs and gen heatmap based on Tomtom using heatmap.py
6. find GC-content, lengths, number of motifs using avg_org_motifs.py
7. df of GC-content vs. lengths
8. find pairwise divergence time using cophenetic.phylo() (and divide by 2)
9. use iTol and ape to create customized trees
10. df of similarity vs. divergence using div_time_similarity.py
11. find correlation of length vs GC using phy_regress.r
12. generate dataframe containing diff GC and diff length of all species pairs
13. find correlation of 3 relationshps using 
