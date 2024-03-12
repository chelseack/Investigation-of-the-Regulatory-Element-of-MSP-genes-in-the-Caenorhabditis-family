# make a df with cluster, GC of each cluster, length of each cluster, divergence and conservation
length_all=[]
GC_all=[]
clusters=['AC0001', 'AC0002', 'AC0008', 'AC0009', 'AC0013', 'AC0041', 'AC0053',
 'AC0062', 'AC0065', 'AC0074', 'AC0079', 'AC0083', 'AC0084', 'AC0086', 'AC0092',
 'AC0095', 'AC0104', 'AC0107', 'AC0108']
cons=[5/12,18/45,8/29,10/46,5/29,4/46,4/43,4/29,3/43,8/31,3/43,3/45,4/45,4/45,3/45,
      3/47,8/45,10/43,11/45]
divergence=[47.378,128.053,72.39,164.249,72.39,164.249,108.238,72.39,108.238,
            108.045,108.238,128.053,128.053,128.053,128.053,164.767,128.053,108.238,128.053]
#cluster = 'AC0095'
for cluster in clusters:
    df = motif_annot_df.groupby('cluster').get_group(cluster)
    
    pwm = process_cluster(df, tomtom)

    # Extracting values using regular expressions
    pattern = r'(\d+\.\d+)'
    data = {}
    for line in pwm.strip().split('\n'):
        parts = line.split(':')
        nucleotide = parts[0].strip()
        values = re.findall(pattern, parts[1])
        data[nucleotide] = [float(val) for val in values]
    
    df = pd.DataFrame(data)
    # Summing up values in the 'C' and 'G' columns
    sum_GC = df['C'].sum() + df['G'].sum()
    GC=sum_GC/len(df)
    print(GC)
    GC_all.append(GC)
    length_all.append(len(df))

df = pd.DataFrame({'Cluster':clusters, 'GC': GC_all, 'length': length_all, 
                   'Conservation':cons, 'Divergence':divergence})

# Plot all pairs of variables
sns.pairplot(df)
plt.show()
