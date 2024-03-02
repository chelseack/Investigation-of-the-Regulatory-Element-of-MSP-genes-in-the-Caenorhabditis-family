# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:33:05 2024

@author: Chelsea
"""
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
import matplotlib.pyplot as plt
import seaborn as sns
import os.path
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgridspec
import os
import glob
import logomaker as lm
#from genome_tools.plotting import sequence

tomtom = r"C:\Users\kenne\Downloads\EEB498\Tomtom\tomtom.tsv"
pfm = r"C:\Users\kenne\Downloads\EEB498\pfm"
EEB498 = r"C:\Users\kenne\Downloads\EEB498\cluster"

tomtom = pd.read_table(tomtom)

# Pivot table to create a similarity matrix
sim = tomtom.pivot_table(index='Query_ID', columns='Target_ID', values='E-value', fill_value=np.nan)
#print("3333", sim)
# Process the similarity matrix
x = sim.values
w = np.triu(x) + np.triu(x, 1).T
v = np.tril(x) + np.tril(x, -1).T
sim.iloc[:, :] = np.nanmin(np.dstack([w, v]), axis=2)
sim.fillna(100, inplace=True)
sim = -np.log10(sim)
sim[np.isinf(sim)] = 10
# #print("22222",sim)
# Cluster the square matrix
Z = linkage(sim, method='complete', metric='correlation')
cl = fcluster(Z, 0.7, criterion='distance')
o = dendrogram(Z, no_plot=True)['leaves']

#print(f'Number of motif clusters: {max(cl)}')

motif_annot_df = pd.DataFrame({'motif_id':sim.index, 'cluster':cl})
motif_annot_df['cluster'] = 'AC' + motif_annot_df['cluster'].astype(str).str.zfill(4)

pivoted_df = motif_annot_df.pivot(columns='cluster', values='motif_id')
pivoted_df = pivoted_df.reindex(sorted(pivoted_df.columns), axis=1)
pivoted_df = pivoted_df.rename_axis(index='motifs').reset_index()
#pivoted_df.to_csv("pivoted_df.tsv", index=False)
cluster_list = list(motif_annot_df['cluster'])
# print("11111",sim.iloc[0,0])
#sim.iloc[0,0].to_csv("clustermap.tsv")
data= r"C:\Users\kenne\Downloads\EEB498\matrix_heatmap.csv"
df = pd.read_csv(data, index_col=0)
labels_y = df.index
labels_x = df.columns
df_numeric = df.iloc[1:, 1:].astype(float)

# Heatmap plot 
sns.clustermap(df_numeric, row_cluster=False, col_cluster=False, cmap='Blues')
plt.xticks(ticks=range(len(labels_x)), labels=labels_x, rotation=45, ha='right')
plt.yticks(ticks=range(len(labels_y)), labels=labels_y)
plt.show()

# Heatmap plot 
#sns.clustermap(sim.iloc[o, o], vmin=0, vmax=5, row_cluster=False, col_cluster=False, cmap='Blues')

# Generate clusters and logo
cluster = pd.DataFrame()
def relative_info_content(pwm):
    p = pwm/np.sum(pwm, axis = 1)[:,np.newaxis]
    ic = 2+np.sum(p*np.nan_to_num(np.log2(p)), axis = 1)
    ric = p*ic[:,np.newaxis]
    return ric


def rev_compl(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))

import warnings
warnings.filterwarnings('ignore')

def abs_mean(x):
    return np.mean(np.abs(x))

def process_cluster(df, tomtom_df):  
    print("111111")
    print(df)
    motifs = df["motif_id"]
    motifs_index = motifs.index
    rows = (tomtom_df['Query_ID'].isin(motifs)) | (tomtom_df['Target_ID'].isin(motifs))
    all_pairwise_df = tomtom_df[rows]      

      # Seed motif has the best optimal offset in the group
    seed_motif = all_pairwise_df.groupby('Query_ID').agg({'Overlap': abs_mean}).sort_values('Overlap', ascending=False).index[0]
    seed_motif = all_pairwise_df.groupby('Query_ID').agg({'E-value': np.median}).sort_values('E-value', ascending=True).index[0]
    
    rows = (tomtom_df['Query_ID'] == seed_motif) & (tomtom_df['Target_ID'].isin(motifs))
    pairwise_df = tomtom_df[rows] 
    pivot_df = all_pairwise_df.pivot_table(index='Query_ID', columns='Target_ID', values='Optimal_offset')
    q = pivot_df.loc[seed_motif]
    qi = q[~np.isnan(q)].index
    
    query = []
    target = []
    offset = []
    orientation = []
    target_consensus = []
    query_consensus = []

    for m in motifs[~motifs.isin(pairwise_df['Target_ID'])]:
        t = pivot_df[m]
        ti = t[~np.isnan(t)].index
        try:
            common_motif = (qi&ti)[0]
            print("common motif: "+common_motif)
        except:
            print(f'ERROR: {m} no alignments available!')
            continue
  
        row_q = (all_pairwise_df['Query_ID'] == seed_motif) & (all_pairwise_df['Target_ID'] == common_motif)
        row_t = (all_pairwise_df['Query_ID'] == common_motif) & (all_pairwise_df['Target_ID'] == m)

        
        offset_q = all_pairwise_df[row_q]['Optimal_offset'].iloc[0]
        offset_t = all_pairwise_df[row_t]['Optimal_offset'].iloc[0]
        
        orientation_q = all_pairwise_df[row_q]['Orientation'].iloc[0]
        orientation_t = all_pairwise_df[row_t]['Orientation'].iloc[0]
        
        consensus_q = all_pairwise_df[row_q]['Query_consensus'].iloc[0]
        
        consensus_tq = all_pairwise_df[row_t]['Query_consensus'].iloc[0]
        consensus_tt = all_pairwise_df[row_t]['Target_consensus'].iloc[0]
        
        offset_p = len(consensus_tt) - offset_t - len(consensus_tq)
        
        target.append(m)
        query_consensus.append(consensus_q)
        print(query_consensus)
        print(target_consensus)

        if orientation_t == orientation_q:
            orientation.append('+')
            target_consensus.append(consensus_tt)
            
            if orientation_t == '+':
                offset.append(offset_q+offset_t)
            else:
                offset.append(offset_p+offset_q)
        else:            
            orientation.append('-')
            target_consensus.append(rev_compl(consensus_tt))
            
            if orientation_q == '-':        
                offset.append(offset_p+offset_q)
            else:
                offset.append(offset_q+offset_t)
        
    z = pd.DataFrame({
        'Query_ID': seed_motif,
        'Target_ID': target,
        'Optimal_offset': offset,
        'p-value': 0, 
        'E-value': 0,
        'q-value': 1,
        'Overlap': 0,
        'Query_consensus': query_consensus,
        'Target_consensus': target_consensus,
        'Orientation': orientation,
    })
    
    
    if len(z) > 0:
        pairwise_df = pd.concat([pairwise_df, z])
    
    w = pairwise_df['Target_consensus'].str.len()
    left = min(-pairwise_df['Optimal_offset'])
    l_offset = -left - pairwise_df['Optimal_offset']
    right = max(l_offset + w)
    r_offset = right - w - l_offset
    alignment_df = pairwise_df.drop(['Query_ID', 'Optimal_offset', 'p-value', 'E-value', 'q-value', 'Overlap', 'Query_consensus'], axis=1)
    alignment_df.loc[:,'w'] = w
    alignment_df.loc[:,'l_offset'] = l_offset
    alignment_df.loc[:,'r_offset'] = r_offset
    alignment_df.columns = ['motif', 'consensus', 'strand', 'w', 'l_offset', 'r_offset']
    
    alignment_df.reset_index(drop=True, inplace=True)

    alignment_df = alignment_df.merge(df.reset_index(), left_on='motif', right_on='motif_id')
    
    alignment_df.sort_values(by=['w'], inplace=True)
    alignment_df.reset_index(inplace=True)
    # print("alignment df")
    # print(alignment_df)
    n = len(alignment_df)
    l = min(alignment_df['l_offset'])
    r = max(alignment_df['r_offset'] + alignment_df['w'])
    w = r - l
    summed_pwm = np.zeros((4, int(w), n))
    print("alignment df ")
    print(alignment_df)
    for i, row in alignment_df.iterrows():
        #print(row['motif'])
        motif_id = row['motif']
        rc = row['strand'] == '-'
        left = int(row['l_offset'])
        width = int(row['w'])
        file_names = os.listdir(pfm)
        motif_pfm1 = [os.path.join(pfm, file_name) for file_name in file_names if motif_id in file_name]
        if motif_pfm1:
            motif_pfm = motif_pfm1[0]
            pwm = np.loadtxt(motif_pfm, skiprows=1, dtype = str)
            #pwm = pwm.T
        else:
            print("due to pwm found before asignment",cluster, motif_id)
        if rc:
            pwm = pwm[::-1,::-1]
        else:
            #print("due to pwm found before asignment",cluster, motif_id)
            print("1324")
        
        extended_pwm = np.ones((4, int(w))) * 0.25
        transposed_pwm = list(zip(*pwm))
        pwm1 = [list(sublist) for sublist in transposed_pwm]
        try:
            extended_pwm[:,left:left+width] = pwm1
            # print(motif_id)
            # print(motif_pfm1)
            # print(pwm)
        except:
            print(motif_id)
            # print(motif_pfm1)
            # print(pwm)
            # print(transposed_pwm)
            # print(pwm1)
            # print(extended_pwm[:,left:left+width])
            
        summed_pwm[:,:,i] += extended_pwm
        # print("summed pwm ")
        # print(summed_pwm)
    avg_pwm = np.nanmean(summed_pwm, axis=2).T
    
    ic = relative_info_content(avg_pwm)
    total_ic = ic.sum(axis=1)

    cdf = np.cumsum(total_ic)/np.sum(total_ic)
    s = np.where(cdf > 0.05)[0][0]
    e = np.where(cdf > 0.95)[0][0] + 1    

    avg_pwm = avg_pwm[s:e,:]
    #print("avg pwm ")
    #print(avg_pwm)
    ## plot
    
    fig = plt.figure()
    fig.set_size_inches((w+2)*.125+2, (n+1)*0.5+1)
    
    gs = mgridspec.GridSpec(n+1, 1)
    logos = []
    id_motif = []
    cluster_species={}
    new_folder_path = os.path.join(EEB498, cluster)
    os.makedirs(new_folder_path, exist_ok=True)

    for i, row in alignment_df.iterrows():
        ax = fig.add_subplot(gs[i+1, :])
        
        motif_id = row['motif']
        rc = row['strand'] == '-'
        left = row['l_offset']
        width = row['w']

        file_names = os.listdir(pfm)
        motif_pfm1 = [os.path.join(pfm, file_name) for file_name in file_names if motif_id in file_name]
        try:
            file_name = os.path.splitext(os.path.basename(motif_pfm1[0]))[0]    
        except IndexError():
            print(motif_pfm1[0])
        if motif_pfm1:
            id_motif.append(file_name)
            motif_pfm = motif_pfm1[0]
            pwm = np.loadtxt(motif_pfm, skiprows=1, dtype = str)
    
        if rc:
            pwm = pwm[::-1,::-1]
        
        array_data = np.array(pwm).astype(float)
        # print(array_data)
        # print(pwm)
        # print(pwm.T)
        df = pd.DataFrame(array_data, columns=['A', 'C', 'G', 'T'])
        df['pos'] = range(len(df))
        df.set_index('pos', inplace=True)
        logo = lm.Logo(df)
        # logo.ax.set_xlabel(file_name, fontsize = 15)

        ax.axvspan(l-1, s, fc='lightgrey', alpha=0.5)
        ax.axvspan(e, r+1, fc='lightgrey', alpha=0.5)
        
        ax.set_xlim(left=l-1, right=r+1)
        ax.set_ylim(bottom=0, top=2.1)

        ax.xaxis.set_visible(False)
        ax.set_yticks([])
        
        #source_id = str(row['source_id'])
        #source_id = source_id[:10] + '...' if len(source_id) > 10 else source_id
        #tf_name = str(row['tf_name'])  + '(' + str(row['motif_type']) + ')'

        #ax.set_ylabel(tf_name + '\n (' + source_id + ')', rotation=0, ha='right', va='center', fontname="IBM Plex Mono", fontsize='medium')
    arche_logo = []
    # Archetype motif
    ax = fig.add_subplot(gs[0,:])

    avg_pwm_nonnp = np.array(avg_pwm)
    df_avg_pwm = pd.DataFrame(avg_pwm_nonnp, columns=['A', 'C', 'G', 'T'])
    df_avg_pwm.index.name = 'pos'
    out_df_avg_pwm = pd.DataFrame(df_avg_pwm)
    txt=cluster+".txt"
    out_df_avg_pwm.to_csv(os.path.join(new_folder_path, txt), index=False)
    logo = lm.Logo(out_df_avg_pwm, show_spines=None)
    #logo.ax.set_title(id_motif, fontsize = 15)
    #plt.savefig(os.path.join(new_folder_path, str(cluster)))

    ax.set_xlim(left=l-1, right=r+1)
    ax.set_ylim(bottom=0, top=2.1)
    ax.xaxis.set_visible(False)
    ax.set_yticks([])

    ax.axvspan(s, e, fc='none', ec='r', lw=2, clip_on=False)
    [ax.spines[loc].set_visible(False) for loc in ['top', 'bottom', 'left', 'right']]

    ax.set_ylabel('Archetype\nconsensus', rotation=0, ha='right', va='center', fontsize='large', fontweight='bold', color='r')
    
    cluster_id = str(alignment_df['cluster'][0])
    #gene_family = alignment_df['tf_name'].str.replace('[\-0-9]+$', '').value_counts().index[:2].str.cat(sep='/')
    #dbd = alignment_df['family_name'].astype(str).value_counts().index[0].replace(' ', '_')
    #cluster_name = cluster_id + ':' + gene_family + ':' + dbd

    figw, figh = fig.get_size_inches()
    height_frac = (figh-0.75)/figh
    
   
    
    gs.update(left=1-((figw-1.75)/figw), right=(figw-0.25)/figw, top=(figh-0.75)/figh, bottom=1-((figh-0.25)/figh))
    
    fig.suptitle(cluster_id.upper(), fontweight='bold', fontsize='large', y=1-(.5/figh), va='center')
    #plt.savefig(f'results/clusters/{cluster_id}.pdf')
    #plt.savefig(f'results/clusters/{cluster_id}.png')
    
    #
    w = avg_pwm.shape[0]
    
    fig = plt.figure()
    #fig.set_size_inches(w*0.125+0.5, 0.75)
    fig.set_size_inches(w*0.125, 0.5)
    
    figw, figh = fig.get_size_inches() 
    
    gs = mgridspec.GridSpec(1, 1)
    #gs.update(left=1-((figw-0.25)/figw), right=(figw-0.25)/figw, top=1-(0.25/figh), bottom=0)
    gs.update(left=0, right=1, top=1, bottom=0)

    ax = fig.add_subplot(gs[:,:])
    
    logo = lm.Logo(out_df_avg_pwm)
    
    ax.set_xlim(left=0, right=w)
    ax.set_ylim(bottom=0, top=2.1)
    
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    
    [ax.spines[loc].set_visible(False) for loc in ['top', 'bottom', 'left', 'right']]
    
    #plt.savefig(f'results/clusters/logos/{cluster_id}.pdf')
    #plt.savefig(f'results/clusters/logos/{cluster_id}.png')
    
    #header_line =  cluster_name + '\n'
    mat = pd.DataFrame(avg_pwm.T, index=['A:', 'C:', 'G:', 'T:']).to_string(header=False)
    return mat

cluster = 'AC0095'

df = motif_annot_df.groupby('cluster').get_group(cluster)

pwm = process_cluster(df, tomtom)
