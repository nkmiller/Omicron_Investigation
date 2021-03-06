#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import csv
import networkx as nx

"""
Networking Data
"""
# set path to CRM Matrices folder
CSV_path = './CRM_Matrices/'

# Helpers from https://stackoverflow.com/questions/14247586/how-to-select-rows-with-one-or-more-nulls-from-a-pandas-dataframe-without-listin 
def row_nan_sums(df):
    sums = []
    for row in df.values:
        sum = 0
        for el in row:
            if el > 0.05: 
                sum+=1
        sums.append(sum)
    return sums
def query_k_plus_sums(df, k, include_list):
    sums = row_nan_sums(df)
    indices = []
    i = 0
    for sum in sums:
        current_resi = df.index[i]
        if (sum >= k) or (current_resi in include_list):
            indices.append(True)
        else:
            indices.append(False)
        i += 1
    return indices


SIN_Matrices = ['6xdg_REGN10933_SINX_mat.csv','6xdg_REGN10987_SINX_mat.csv','7b3o_SINX_mat.csv',
                '7byr_SINX_mat.csv','7bz5_SINX_mat.csv','7c01_SINX_mat.csv','7c8v_SINX_mat.csv',
                '7c8w_SINX_mat.csv','7cah_SINX_mat.csv','7cdi_SINX_mat.csv','7cdj_SINX_mat.csv',
                '7ch4_SINX_mat.csv','7ch5_SINX_mat.csv','7chh_SINX_mat.csv','7jmo_SINX_mat.csv',
                '7jmw_SINX_mat.csv','7jva_SINX_mat.csv','7jvb_SINX_mat.csv','7bwj_SINX_mat.csv',
                '7jx3_S2H14_SINX_mat.csv','7jx3_S304_SINX_mat.csv','7jx3_S309_SINX_mat.csv','7k8s_SINX_mat.csv',
                '7k8u_SINX_mat.csv','7k8v_SINX_mat.csv','7k8w_SINX_mat.csv','7k8z_SINX_mat.csv',
                '7k9z_F298_SINX_mat.csv','7k45_SINX_mat.csv','7k90_SINX_mat.csv','6xc7_SINX_mat.csv',
                '6zer_SINX_mat.csv','7k43_SINX_mat.csv','7ld1_SINX_mat.csv','7kmi_SINX_mat.csv','7kmh_SINX_mat.csv',
                '7kzb_SINX_mat.csv','7dpm_SINX_mat.csv','7bel_SINX_mat.csv','7kmg_SINX_mat.csv', 
                '7L7E_AZD1061_SINX_mat.csv','7L7E_AZD8895_SINX_mat.csv','7RAL_S2X259_SINX_mat.csv',
                '7mzk_SINX_mat.csv','7mkm_SINX_mat.csv','7ezv_SINX_mat.csv','7lm8_SINX_mat.csv',
                '7det_SINX_mat.csv','7ey5_SINX_mat.csv'
                ]

matrix_to_mAb = {'6xdg_REGN10933_SINX_mat.csv':'REGN10933','6xdg_REGN10987_SINX_mat.csv':'REGN10987',
                 '7b3o_SINX_mat.csv':'COR-101','7c8v_SINX_mat.csv':'SR4','7bwj_SINX_mat.csv':'BRII-196',
                '7byr_SINX_mat.csv':'BD23','7bz5_SINX_mat.csv':'B38','7c01_SINX_mat.csv':'LY-CoV016',
                '7c8w_SINX_mat.csv':'MR17','7cah_SINX_mat.csv':'H014','7cdi_SINX_mat.csv':'P2C-1F11',
                '7cdj_SINX_mat.csv':'P2C-1A3','7jmo_SINX_mat.csv':'COVA2-04','7k8s_SINX_mat.csv':'C002',
                '7ch4_SINX_mat.csv':'BD-604','7ch5_SINX_mat.csv':'BD-629','7chh_SINX_mat.csv':'BD-368-2',
                '7jmw_SINX_mat.csv':'COVA1-16','7jva_SINX_mat.csv':'S2A4','7jvb_SINX_mat.csv':'nb20',
                '7jx3_S2H14_SINX_mat.csv':'S2H14','7jx3_S304_SINX_mat.csv':'S304','7jx3_S309_SINX_mat.csv':'S309',
                '7k8u_SINX_mat.csv':'C104','7k8v_SINX_mat.csv':'C110','7k8w_SINX_mat.csv':'C119',
                '7k9z_F298_SINX_mat.csv':'298','7k45_SINX_mat.csv':'S2E12','7k90_SINX_mat.csv':'C144',
                '6xc7_SINX_mat.csv':'CR3022','7kmg_SINX_mat.csv':'LY-CoV555','7k8z_SINX_mat.csv':'C135',
                '6zer_SINX_mat.csv':'EY6A','7k43_SINX_mat.csv':'S2M11','7RAL_S2X259_SINX_mat.csv':'S2X259',
                '7ld1_SINX_mat.csv':'DH1047','7kmi_SINX_mat.csv':'LY-CoV481','7kmh_SINX_mat.csv':'LY-CoV488',
                '7dpm_SINX_mat.csv':'MW06','7bel_SINX_mat.csv':'COVOX-45','7kzb_SINX_mat.csv':'CR3014-C8',
                '7L7E_AZD1061_SINX_mat.csv':'AZD1061','7L7E_AZD8895_SINX_mat.csv':'AZD8895',
                 '7mzk_SINX_mat.csv':'PDI 96','7mkm_SINX_mat.csv':'SARS2-38','7ezv_SINX_mat.csv':'BD-812',
                 '7lm8_SINX_mat.csv':'CV38-142','7det_SINX_mat.csv':'PR961','7ey5_SINX_mat.csv':'BD-821'}

df_paratope_networking = pd.DataFrame() 
for csv in SIN_Matrices:
    temp_df = pd.read_csv(CSV_path + csv) # convert csv to df
    store_df_paratope = {} # library to store network data
    temp_df = temp_df.set_index('Observations') # set index
    temp_df.columns = temp_df.index # convert AA indices to col indices
    for row in temp_df.index: #iterate through the rows
        resnum = row[3:-2] # numerical res only
        resnum = ''.join(filter(str.isdigit, row))
        # check if RBD residue
        if (331 <= int(resnum) <= 531):
            paratope_sum = 0 # track connection to paratope
            epitope_sum = 0  # track connection within epitope 
            for col in temp_df.columns:
                if (col[-1] != row[-1]): # if interaction is between RBD and a different chain = mAbs
                    paratope_sum += temp_df[row][col]
                else:  #otherwise, if interaction within the RBD
                    epitope_sum += temp_df[row][col]
            if paratope_sum > 0: # else if paratope+epitope and we have at least one paratope connection
                store_df_paratope[str(resnum)] = paratope_sum
            else: #if no connection to paratope, store zero
                store_df_paratope[str(resnum)] = 0
    df_paratope_networking[matrix_to_mAb[csv]] = pd.DataFrame.from_dict(store_df_paratope,orient='index')[0]
df_paratope_networking = df_paratope_networking.fillna(0)
df_paratope_networking_full_data = df_paratope_networking.fillna(0) 

glycine_pos = [339,381,404,413,416,431,446,447,476,482,485,496,502,504] # for WT RBD epitopes
df_paratope_networking_normalized =(df_paratope_networking-df_paratope_networking.min())/(df_paratope_networking.max()-df_paratope_networking.min())
df_paratope_networking_indices = df_paratope_networking_normalized[query_k_plus_sums(df_paratope_networking_normalized,1, [])]
indirect_paratope_indices = df_paratope_networking_indices.index

#repeat previous for indirect networking
df_paratope_networking_linkage = pd.DataFrame() 
df_paratope_networking_linkage_count = pd.DataFrame() 

df_Omicron_S309_networking = pd.DataFrame(np.zeros((len(df_paratope_networking_normalized.index), len(df_paratope_networking_normalized.index)))) # track linkage networking from Omicron sites to all direct networked sites on an Ab paratope
df_Omicron_S309_networking.index = df_paratope_networking_normalized.index
df_Omicron_S309_networking.columns = df_paratope_networking_normalized.index
df_Omicron_ADG2_networking = df_Omicron_S309_networking.copy() 
df_Omicron_AZD1061_networking = df_Omicron_S309_networking.copy()
df_Omicron_AZD8895_networking = df_Omicron_S309_networking.copy()

for csv in SIN_Matrices:
    temp_df = pd.read_csv(CSV_path + csv) # convert csv to df
    store_df_paratope_linkage = {}
    store_df_paratope_linkage_count = {}
    temp_df = temp_df.set_index('Observations') # set index
    temp_df.columns = temp_df.index # convert AA indices to col indices
    for row in temp_df.index: #iterate through the rows
        resnum = row[3:-2] # numerical res only
        resnum = ''.join(filter(str.isdigit, row)) 
        # check if RBD
        if (331 <= int(resnum) <= 531):
            paratope_sum = 0 # track connection to paratope
            epitope_sum = 0  # track connection to epitope 
            epitope_sum_count = 0
            for col in temp_df.columns:
                resnum_interact = col[3:-2] # numerical res only
                if (col[-1] != row[-1]): # if interaction is between RBD and a different chain = mAbs
                    paratope_sum += temp_df[row][col]
                else:  #otherwise, if interaction within the RBD
                     if resnum_interact in indirect_paratope_indices[df_paratope_networking_indices[matrix_to_mAb[csv]] != 0]:# if residue is a paratope networked residue
                        epitope_sum += temp_df[row][col]
                        if temp_df[row][col] > 0: # if there's networking between sites
                            epitope_sum_count += 1 #count as an indirect interaction paratope
                        # if a residue is networked to the given paratope, we get the within RBD networking to that site
                        if 'S309' in csv:
                            df_Omicron_S309_networking.loc[row[3:-2],col[3:-2]] =  temp_df[row][col]
                        elif 'AZD1061' in csv:
                            df_Omicron_AZD1061_networking.loc[row[3:-2],col[3:-2]] =  temp_df[row][col]
                        elif 'AZD8895' in csv:
                            df_Omicron_AZD8895_networking.loc[row[3:-2],col[3:-2]] =  temp_df[row][col]
                     ADG2_critical_network = ['402','403','404','405','406','453','495','496','497','501','502','503','504','505','506'] #The 1st order network from the four ADG2 critical sites
                     if (resnum_interact in ADG2_critical_network) and ('S309' in csv): #use S309 PDB for ADG-2 within RBD networking as S309 binds away from ADG2 epitope as inidcated by Figure S7 Rappazzo et al
                          df_Omicron_ADG2_networking.loc[row[3:-2],col[3:-2]] =  temp_df[row][col]
            if (epitope_sum > 0) or (int(resnum) in glycine_pos):
                store_df_paratope_linkage[str(resnum)] = epitope_sum
                store_df_paratope_linkage_count[str(resnum)] = epitope_sum_count
            else: #if no connection to paratope, store zero
                store_df_paratope_linkage[str(resnum)] = 0 # store if we have a noznero interaction with epitope residue
                store_df_paratope_linkage_count[str(resnum)] = 0
    df_paratope_networking_linkage[matrix_to_mAb[csv]] = pd.DataFrame.from_dict(store_df_paratope_linkage,orient='index')[0]
    df_paratope_networking_linkage_count[matrix_to_mAb[csv]] = pd.DataFrame.from_dict(store_df_paratope_linkage_count,orient='index')[0]

df_paratope_networking_linkage_full_data = df_paratope_networking_linkage.fillna(0) 
df_paratope_networking_linkage_count_full_data = df_paratope_networking_linkage_count.fillna(0)

#glycine networking, within 5A of G ac
adjacent_resi = {}
adjacent_resi[339] = [337,338,340,341,342,343]
adjacent_resi[381] = [380, 382, 430]
adjacent_resi[404] = [402,403,405,406,407]
adjacent_resi[413] = [411,412,414,424,427]
adjacent_resi[416] = [409,415,416,417,418,419,420]
adjacent_resi[431] = [380,429,430,432,513,514,515]
adjacent_resi[446] = [498,444,445,447]
adjacent_resi[447] = [444,445,446,448,449,497,498]
adjacent_resi[476] = [474,475,477,478,487]
adjacent_resi[482] = [472,480,481,483]
adjacent_resi[485] = [483,484,486,488]
adjacent_resi[496] = [448,449,494,495,497,498]
adjacent_resi[502] = [500,501,503,504,505,506]
adjacent_resi[504] = [502,503,505,506]

def glycine_network(linkage_data):
    temp_linkage = linkage_data.copy().fillna(0)
    linkage_data2 = linkage_data.copy().fillna(0)
    for i in linkage_data.index:
        if int(i) in glycine_pos:
            adjacent_list = adjacent_resi[int(i)]
            temp_score = temp_linkage.loc[i] # raw value
            for gly in adjacent_list:
                temp_score += temp_linkage.loc[str(gly)] / (len(adjacent_list)+1)
            linkage_data2.loc[str(i)] = temp_score 
    return linkage_data2

# glycine correction for indirect networking
df_paratope_networking_linkage_full_data = glycine_network(df_paratope_networking_linkage_full_data)
df_paratope_networking_linkage_count_full_data = glycine_network(df_paratope_networking_linkage_count_full_data)
df_Omicron_S309_networking2 = glycine_network(df_Omicron_S309_networking)
df_Omicron_ADG2_networking2 = glycine_network(df_Omicron_ADG2_networking)
df_Omicron_AZD1061_networking2 = glycine_network(df_Omicron_AZD1061_networking)
df_Omicron_AZD8895_networking2 = glycine_network(df_Omicron_AZD8895_networking)

"""
Full heat map
"""
# mAb Heatmap for Direct AND Indirect
df_paratope_networking_full_data_temp = df_paratope_networking_full_data + df_paratope_networking_linkage_full_data
df_paratope_networking_full_data_normalized = (df_paratope_networking_full_data_temp-df_paratope_networking_full_data_temp.min())/(df_paratope_networking_full_data_temp.max()-df_paratope_networking_full_data_temp.min())
df_paratope_networking_normalized_enriched = df_paratope_networking_full_data_normalized[query_k_plus_sums(df_paratope_networking_full_data_normalized,3, [])]
sns.set(font_scale=0.8)
g = sns.clustermap(df_paratope_networking_normalized_enriched.fillna(0), dendrogram_ratio=0.2, metric='canberra', figsize=(12, 18), row_cluster=True, cmap="rocket_r",yticklabels=True)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 5)
g.ax_heatmap.invert_yaxis()
g.savefig('CRM Fig1.png', dpi=300)

"""
Variants Networks
"""
def graph_network(df_networking, muts, ADG2):
    
    G = nx.Graph()
    used_node = []
    # index (rows) are sites networking to sites directly networked to paratope
    for resi in df_networking.index: 
        # columns are sites directly networked to the paratope
        for resi2 in df_networking.columns: 
            try:
                networking = 0.5*df_networking.loc[resi,resi2] + 0.5*df_networking.loc[resi2,resi]
            except:
                networking = 0
            if networking:
                G.add_weighted_edges_from([(resi, resi2, (networking) )]) 
                if resi not in used_node:
                    used_node.append(resi)
                if resi2 not in used_node:
                    used_node.append(resi2)
                        
    node_order = used_node
    # size represents direct networking
    node_antigenicity = []
    for resi in node_order:
        node_antigenicity.append(1)
        
        
    node_colors = []
    df_networking_x = df_networking.loc[~(df_networking==0).all(axis=1)].copy()
    df_networking_x = df_networking_x.loc[:, (df_networking_x != 0).any(axis=0)].copy()
    if ADG2:
        ADG2_critical = ['405','502','504','505']
        for node in node_order:
            if (int(node) in muts) and (node in ADG2_critical):
                node_colors.append('lightsalmon')
            elif int(node) in muts:
                node_colors.append('plum')
            elif node in ADG2_critical: # no public structure, only use 4 critical knockout sites for direct network
                node_colors.append('palegreen')#palegreen
            else:
                node_colors.append('lightblue')#lightblue
    else:
        for node in node_order:
            if int(node) in muts:
                node_colors.append('plum')
            elif node in df_networking_x.columns:
                node_colors.append('palegreen')
            else:
                node_colors.append('lightblue')#lightblue
        
    #plotting 
    start_seed = 0
    posx = nx.spring_layout(G, iterations = 1000, weight='weight',seed=start_seed,k=10/np.sqrt(G.order()))
    fig_network = plt.figure(figsize = (10,10))
    NX = nx.draw(G, with_labels=True, nodelist = node_order, width=0.5, node_size= [i * 1000 for i in node_antigenicity], font_weight = 'bold', font_family = 'arial', font_size = 12, node_color=node_colors, edgecolors='black', pos= posx) 
    return fig_network

# plot for therapeutic mAbs
Omicron_muts = [339,371,373,375,417,440,446,477,478,484,493,496,498,501,505]
fig_network = graph_network(df_Omicron_S309_networking2, Omicron_muts, False)
fig_network.savefig('CRM Fig S1.png', dpi=300)
fig_network = graph_network(df_Omicron_ADG2_networking2, Omicron_muts, True)
fig_network.savefig('CRM Fig S2.png', dpi=300)
fig_network = graph_network(df_Omicron_AZD8895_networking2, Omicron_muts, False)
fig_network.savefig('CRM Fig S3.png', dpi=300)
fig_network = graph_network(df_Omicron_AZD1061_networking2, Omicron_muts, False)
fig_network.savefig('CRM Fig S4.png', dpi=300)


"""
OMICRON anitgenic residues only
"""
df_i = df_paratope_networking_linkage_full_data  # indirect 
df_d = df_paratope_networking_full_data  # direct

temp_beta_i = df_i.loc['417'] + df_i.loc['484'] + df_i.loc['501']
temp_beta_d = df_d.loc['417'] + df_d.loc['484'] + df_d.loc['501']
temp_beta_i.name = 'Beta'#'417+484+501'
temp_beta_d.name = 'Beta'#'417+484+501'
df_d_VOC = df_d.append(temp_beta_d)
df_i_VOC = df_i.append(temp_beta_i)

temp_Omi_d = df_d.loc['371'] + df_d.loc['375'] + df_d.loc['373'] + df_d.loc['339'] + df_d.loc['417'] + df_d.loc['484'] + df_d.loc['501'] + df_d.loc['440'] + df_d.loc['446'] + df_d.loc['493'] + df_d.loc['496'] + df_d.loc['498'] + df_d.loc['505'] + df_d.loc['477'] + df_d.loc['478'] 
temp_Omi_d.name = 'Omicron' 
df_d_VOC = df_d_VOC.append(temp_Omi_d)
temp_Omi_i = df_i.loc['371'] + df_i.loc['375'] + df_i.loc['373'] + df_i.loc['339'] + df_i.loc['417'] + df_i.loc['484'] + df_i.loc['501'] + df_i.loc['440'] + df_i.loc['446'] + df_i.loc['493'] + df_i.loc['496'] + df_i.loc['498'] + df_i.loc['505'] + df_i.loc['477'] + df_i.loc['478'] 
temp_Omi_i.name = 'Omicron' 
df_i_VOC = df_i_VOC.append(temp_Omi_i)

temp_Delta_d = df_d.loc['452'] + df_d.loc['478'] 
temp_Delta_i = df_i.loc['452'] + df_i.loc['478'] 
temp_Delta_d.name = 'Delta'
temp_Delta_i.name = 'Delta'
df_d_VOC = df_d_VOC.append(temp_Delta_d)
df_i_VOC = df_i_VOC.append(temp_Delta_i)

temp_PMS20_d = df_d.loc['445']+ df_d.loc['455'] + df_d.loc['475'] + df_d.loc['417']+ df_d.loc['440'] + df_d.loc['484'] + df_d.loc['501']  + df_d.loc['346']
temp_PMS20_i = df_i.loc['445']+ df_i.loc['455'] + df_i.loc['475'] + df_i.loc['417']+ df_i.loc['440'] + df_i.loc['484'] + df_i.loc['501']  + df_i.loc['346']
temp_PMS20_d.name = 'PMS20'
temp_PMS20_i.name = 'PMS20'
df_d_VOC = df_d_VOC.append(temp_PMS20_d)
df_i_VOC = df_i_VOC.append(temp_PMS20_i)

"""
Variants and clinical mAbs only
"""
variants_antigenic_combo = ['Omicron','Beta','Delta','PMS20']
clinical_mAbs = ['REGN10933','REGN10987','LY-CoV016','LY-CoV555','S309','AZD1061','AZD8895']
df_d_VOC_mAb = df_d_VOC[clinical_mAbs]
df_i_VOC_mAb = df_i_VOC[clinical_mAbs]

# get enriched set of rows / residues of interest
df_d_VOC_mAb_enrich = query_k_plus_sums(df_d_VOC_mAb,40, variants_antigenic_combo)
df_i_VOC_mAb_enrich = query_k_plus_sums(df_i_VOC_mAb,40, variants_antigenic_combo)
df_d_VOC_mAb2 = df_d_VOC_mAb[df_d_VOC_mAb_enrich]
df_i_VOC_mAb2 = df_i_VOC_mAb[df_i_VOC_mAb_enrich]
#remove zero columns
df_d_VOC_mAb2 = df_d_VOC_mAb2.loc[:, (df_d_VOC_mAb2 != 0).any(axis=0)].fillna(0)
df_i_VOC_mAb2 = df_i_VOC_mAb2.loc[:, (df_i_VOC_mAb2 != 0).any(axis=0)].fillna(0)
# total networking
df_VOC_mAb2 = df_d_VOC_mAb2 + df_i_VOC_mAb2
# normalize
df_d_VOC_mAb2 = (df_d_VOC_mAb2-np.min(df_d_VOC_mAb2.values))/(np.max(df_d_VOC_mAb2.values)-np.min(df_d_VOC_mAb2.values))
df_i_VOC_mAb2 = (df_i_VOC_mAb2-np.min(df_i_VOC_mAb2.values))/(np.max(df_i_VOC_mAb2.values)-np.min(df_i_VOC_mAb2.values))
df_VOC_mAb2 = (df_VOC_mAb2-np.min(df_VOC_mAb2.values))/(np.max(df_VOC_mAb2.values)-np.min(df_VOC_mAb2.values))


# Plot
g = sns.clustermap(df_d_VOC_mAb2, metric='canberra', figsize=(3, 5), row_cluster=True, cmap="rocket_r",yticklabels=True)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 10)
g.savefig('CRM Fig2A.png', dpi=300)
g = sns.clustermap(df_i_VOC_mAb2,  metric='canberra', figsize=(3, 5), row_cluster=True, cmap="rocket_r",yticklabels=True)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 10)
g.savefig('CRM Fig2B.png', dpi=300)
g = sns.clustermap(df_VOC_mAb2,  metric='canberra', figsize=(3, 5), row_cluster=True, cmap="rocket_r",yticklabels=True)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 10)
g.savefig('CRM Fig2C.png', dpi=300)
