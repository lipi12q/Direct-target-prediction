import pandas as pd
import numpy as np
import  os
import re
from operator import *


# compound induced gene expression data
cp_info_combine_1 = pd.read_csv("D:/Py/Cmap/Cmap/cmapPy/CMAP_target/cp_info_combine_1.csv") # This data were collated previously


#genetic gene expression profile data
oe_gene_data = pd.read_csv("D:/Py/Cmap/oe_gene_data.csv")
unique_oe_gene = list(set(oe_gene_data.iloc[:,0]))#4040

xpr_gene_data = pd.read_csv("D:/Py/Cmap/xpr_gene_data.csv")
unique_xpr_gene = list(set(xpr_gene_data.iloc[:,0]))#5157

sh_cgs_gene_data = pd.read_csv("D:/Py/Cmap/sh_cgs_gene_data.csv")
sh_cgs_gene_data = sh_cgs_gene_data.iloc[:4345,:]# after 4345 were removed
unique_sh_cgs_gene = list(set(sh_cgs_gene_data.iloc[:,0]))#4345

#oe数据中有许多基因为假基因或者多个基因名以-连接
target_profile_from_info = []
for gene in unique_gene:
    pro_from = []
    if gene in list(cp_info_combine_1['target']):
        pro_from.append(1)
    else:
        pro_from.append(0)
    if gene in unique_sh_cgs_gene:
        pro_from.append(1)
    else:
        pro_from.append(0)
    if gene in unique_xpr_gene:
        pro_from.append(1)
    else:
        pro_from.append(0)
    if gene in unique_oe_gene:
        pro_from.append(1)
    else:
        pro_from.append(0)
    target_profile_from_info.append(pro_from)
target_profile_from_info = pd.DataFrame(target_profile_from_info)
target_profile_from_info.columns = ['cp','sh','xpr','oe']
target_profile_from_info.index = unique_gene
target_profile_from_info.to_csv("D:/Py/Cmap/Cmap/cmapPy/CMAP_target/target_profile_from_info.csv")



# using R package HGNChelper translates gene symbols into standard symbols，in the file of target_set_name_transinfo.csv
target_set_name1 = pd.read_csv("D:/Py/Cmap/Cmap/cmapPy/CMAP_target/target_set_name_transinfo.csv")
target_set_name1.index = list(target_set_name1['target_set_primary_name'])

#map gene symbols
cp_info_combine_1 = pd.read_csv("D:/Py/Cmap/Cmap/cmapPy/CMAP_target/cp_info_combine_1.csv")

x = target_set_name1.loc[list(cp_info_combine_1['target']),:]
x1 = eq(list(cp_info_combine_1['target']),list(x['target_set_primary_name']))

cp_info_combine_1['target_new'] = list(x['Suggested.Symbol'])
cp_info_combine_1.to_csv("D:/Py/Cmap/Cmap/cmapPy/CMAP_target/cp_info_combine_1.csv")
cp_info_combine_2 = cp_info_combine_1.dropna(subset='target_new')#16802 -- 16801


x = target_set_name1.loc[list(oe_gene_data.iloc[:,0]),:]#4040
print(eq(list(oe_gene_data.iloc[:,0]),list(x['target_set_primary_name'])))
oe_gene_data.iloc[:,0] = list(x['Suggested.Symbol'])
x = x.dropna(subset='Suggested.Symbol')#3590
oe_gene_data = oe_gene_data.loc[[True if i in list(x['Suggested.Symbol']) else False for i in list(oe_gene_data.iloc[:,0])],:]

x = target_set_name1.loc[list(xpr_gene_data.iloc[:,0]),:]#5157
print(eq(list(xpr_gene_data.iloc[:,0]),list(x['target_set_primary_name'])))
xpr_gene_data.iloc[:,0] = list(x['Suggested.Symbol'])#
x = x.dropna(subset='Suggested.Symbol')#5140
xpr_gene_data = xpr_gene_data.loc[[True if i in list(x['Suggested.Symbol']) else False for i in list(xpr_gene_data.iloc[:,0])],:]

x = target_set_name1.loc[list(sh_cgs_gene_data.iloc[:,0]),:]#4345
print(eq(list(sh_cgs_gene_data.iloc[:,0]),list(x['target_set_primary_name'])))
sh_cgs_gene_data.iloc[:,0] = list(x['Suggested.Symbol'])#
x = x.dropna(subset='Suggested.Symbol')#4326
sh_cgs_gene_data = sh_cgs_gene_data.loc[[True if i in list(x['Suggested.Symbol']) else False for i in list(sh_cgs_gene_data.iloc[:,0])],:]

unique_oe_gene = list(set(oe_gene_data.iloc[:,0]))#4040 3571
unique_xpr_gene = list(set(xpr_gene_data.iloc[:,0]))#5157 5139
unique_sh_cgs_gene = list(set(sh_cgs_gene_data.iloc[:,0]))#4345 4324

unique_gene = list(set(cp_info_combine_1['target_new']))#1461
unique_gene.extend(unique_sh_cgs_gene)#5785
unique_gene.extend(unique_xpr_gene)#10924
unique_gene.extend(unique_oe_gene)#14495
unique_gene = list(set(unique_gene))#8234



#oe gene expression values are negatively handled
oe_gene_data1 = -oe_gene_data.iloc[:,1:]
oe_gene_data1.index = oe_gene_data.iloc[:,0]

#consistent columns names
col_name = cp_gene_data_1.columns
xpr_gene_data = xpr_gene_data[col_name]
oe_gene_data1 = oe_gene_data1[col_name[1:]]
sh_cgs_gene_data = sh_cgs_gene_data[col_name]


target_profile_cmap = pd.DataFrame()
for gene in unique_gene:
    gene_profile = pd.DataFrame()
    gene_index =[]
    if gene in list(cp_info_combine_2['target_new']):
        cp_dt_info = cp_info_combine_2[cp_info_combine_2['target_new'] == gene]
        cp_dt_active = list(set(cp_dt_info.loc[cp_dt_info['pert_ult'] == 1,'cmap_name']))
        cpt_active = cp_gene_data_1.loc[[True if g in cp_dt_active else False for g in list(cp_gene_data_1.iloc[:, 0])],:]
        cpt_active_index = cp_gene_data_1.iloc[[True if g in cp_dt_active else False for g in list(cp_gene_data_1.iloc[:, 0])],0]
        cpt_active = cpt_active.iloc[:,1:]
        cpt_active1 = - cpt_active
        cp_dt_inhibit = list(set(cp_dt_info.loc[cp_dt_info['pert_ult'] == 0,'cmap_name']))
        cpt_inhibit = cp_gene_data_1.loc[[True if g in cp_dt_inhibit else False for g in list(cp_gene_data_1.iloc[:, 0])],:]
        cpt_inhibit_index = cp_gene_data_1.iloc[[True if g in cp_dt_inhibit else False for g in list(cp_gene_data_1.iloc[:, 0])], 0]
        cpt_inhibit1 = cpt_inhibit.iloc[:,1:]
        cpt = pd.concat([cpt_active1,cpt_inhibit1],axis=0)
        gene_profile = pd.concat([gene_profile,cpt],axis=0)
        gene_index = list(cpt_active_index) + list(cpt_inhibit_index)
    if gene in unique_sh_cgs_gene:
        sht = sh_cgs_gene_data.loc[[True if g == gene else False for g in list(sh_cgs_gene_data.iloc[:, 0])], :]
        sht = sht.iloc[:, 1:]
        sht.index = ['sht']*sht.shape[0]
        gene_profile = pd.concat([gene_profile,sht],axis=0)
        gene_index.append('sht')
    if gene in unique_xpr_gene:
        xprt = xpr_gene_data.loc[[True if g == gene else False for g in list(xpr_gene_data.iloc[:, 0])], :]
        xprt = xprt.iloc[:, 1:]
        xprt.index = ['xprt']*xprt.shape[0]
        gene_profile = pd.concat([gene_profile,xprt],axis=0)
        gene_index.append('xprt')
    if gene in unique_oe_gene:
        oet = oe_gene_data1.loc[[True if g == gene else False for g in list(oe_gene_data1.index)], :]
        oet.index = ['oet']*oet.shape[0]
        gene_profile = pd.concat([gene_profile,oet],axis=0)
        gene_index.append('oet')

    #gene_profile.index = gene_index
    target_corr = np.corrcoef(gene_profile)
    if gene_profile.shape[0] == 1:
        target_profile1 = gene_profile.T
    else:
        for i in range(target_corr.shape[0]):
            for j in range(target_corr.shape[1]):
                if target_corr[i][j] < 0:
                    target_corr[i][j] = 0.01
                elif i == j:
                    target_corr[i][j] = 0
                else:
                    pass

        corr_sum = np.sum(target_corr, axis=0)
        weights_sum = np.sum(corr_sum)
        weights = [value / weights_sum for value in corr_sum]
        weights_mat = gene_profile.T * weights
        target_profile1 = np.sum(weights_mat, axis=1)
    target_profile_cmap = pd.concat([target_profile_cmap, target_profile1], axis=1)

target_profile_cmap.columns = unique_gene


cp_geneinfo_beta = pd.read_csv("D:/Py/Cmap/geneinfo_beta.txt",sep="\t")
cp_geneinfo_beta = cp_geneinfo_beta[cp_geneinfo_beta['feature_space'] != 'inferred']
cp_geneinfo_beta['gene_id'] = [str(i) for i in list(cp_geneinfo_beta['gene_id'])]

target_profile_cmap1 = target_profile_cmap.loc[list(cp_geneinfo_beta['gene_id']),:]
target_profile_cmap1.index = cp_geneinfo_beta['gene_symbol']
target_profile_cmap1 = target_profile_cmap1.iloc[:,1:]
len(set(cp_geneinfo_beta['gene_symbol']))

target_profile_cmap1.to_csv("D:/Py/Cmap/Cmap/cmapPy/CMAP_target/target_profile_cmap.csv")#cmap gene expression profiles for targets
target_profile_cmap1 = pd.read_csv("D:/Py/Cmap/Cmap/cmapPy/CMAP_target/target_profile_cmap.csv",index_col=0)



#target = target_profile_cmap1.columns[1]
target_up_symbol = []
target_dn_symbol = []
for target in target_profile_cmap1.keys():
    t1 = target_profile_cmap1[target]
    t1 = t1.sort_values()
    t_head = list(t1.head(350).index)#increasing order
    target_dn_symbol.append(t_head)
    t_tail = list(t1.tail(350).index)
    target_up_symbol.append(t_tail)

target_up_symbol = pd.DataFrame(target_up_symbol)
target_dn_symbol = pd.DataFrame(target_dn_symbol)
target_dn_symbol.index = target_profile_cmap1.keys()
target_up_symbol.index = target_profile_cmap1.keys()
target_dn_symbol.to_csv("D:/Py/Cmap/Cmap/cmapPy/CMAP_target/target_dn_symbol.csv")
target_up_symbol.to_csv("D:/Py/Cmap/Cmap/cmapPy/CMAP_target/target_up_symbol.csv")




