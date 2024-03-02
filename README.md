# Investigation-of-the-Regulatory-Element-of-MSP-genes-in-the-Caenorhabditis-family

Workflow:
1. retrieve upstream 100kb genPromoter.sh
2. find motifs: runstreme.sh, optional: use ceqlogo.sh to gen logo
3. find similarity: runTomtom.sh
4. cluster motifs and gen heatmap based on Tomtom using heatmap.py
5. find GC-content, lengths, number of motifs using avg_org_motifs.py
6. df of GC-content vs. lengths
7. find pairwise divergence time using cophenetic.phylo() (and divide by 2)
8. use iTol and ape to create customized trees
9. df of similarity vs. divergence using div_time_similarity.py
10. find correlation of length vs GC using phy_regress.r
11. generate dataframe containing diff GC and diff length of all species pairs
12. find correlation of 3 relationshps using 
