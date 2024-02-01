# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:41:31 2024

@author: jmg21
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import Bio as Bio

##DEFINING AND VIEWING TOTAL DESEQ2 OUTPUTS
#Define the three normalised data sets from DESeq2

vsd = pd.read_csv('vsd.csv', index_col = 0)      
rld = pd.read_csv('rld.csv', index_col = 0)
ntd = pd.read_csv('ntd.csv', index_col = 0)


#Normalising the vsd data set further, to remove outliers (taking row means and sd)
gene_mean = vsd.mean(axis = 1)
gene_sd = vsd.std(axis = 1)

#Ranking the gene means (from the vsd data) and indexing them, in ascending order (of expression level)
gene_rank = gene_mean.rank()                
sort_index = np.argsort(gene_mean) 

#Plotting the variance in gene expression, in ascending order.
plt.plot(gene_rank, gene_sd,'.', alpha = 0.08)
plt.xlabel('Ranked gene number')
plt.ylabel('Mean??')

#Providing a rolling average
sorted_sd = gene_sd[sort_index]
rolling_sd = sorted_sd.rolling(1000, center = True).median()
plt.plot(np.arange(len(gene_rank)), rolling_sd, 'r-')

plt.show()

#Plot the three different data types
def plotmeansd(data, title):   
 
    gene_mean = data.mean(axis = 1)      # row means and sd
    gene_sd = data.std(axis = 1) 
    gene_rank = gene_mean.rank()        # gene rank (by mean counts)
    sort_index = np.argsort(gene_mean)          # index for sorting by rank
    
    plt.plot(gene_rank, gene_sd,'.', alpha = 0.1)
    sorted_sd = gene_sd[sort_index]
    rolling_sd = sorted_sd.rolling(1000, center = True).median()
    plt.plot(np.arange(len(gene_rank)), rolling_sd, 'r-')
    
    plt.xlabel('gene rank(mean)')
    plt.ylabel('standard deviation')
    plt.ylim(0,6)
    plt.title(title)
    
    plt.show()


plotmeansd(ntd, 'normalised transform')
plotmeansd(rld, 'regularised logarithm')
plotmeansd(vsd, 'variance stabilising transformation')



##PERFORMING PCA 
#Transpose the vsd data set (putting samples as rows, genes as columns). 
from Bio.Cluster import pca

tvsd = np.transpose(vsd)                 
ctvsd = tvsd - tvsd.mean(axis=0)    #centre the data by subtracting mean expression data 

#puts the transformed vsd into lower dimensional space (does the PCA here)    
columnmean, coordinates, components, eigenvalues = pca(ctvsd)

#Give the ratio (%) of the total variance which is spread across each PC
sume = sum(eigenvalues)
ratio = eigenvalues / sume

#Give the % of total variance accounted for by each PC
cumratio = np.cumsum(ratio)

#Plot the variance against cumulative variance
t = range(1,28)

fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('PC')
ax1.set_ylabel('% Variance', color=color)
ax1.plot(t, cumratio, color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xticks(t[:5])

ax2 = ax1.twinx()  
color = 'tab:blue'
ax2.set_ylabel('Contribution to variance', color=color)
ax2.plot(t, ratio, color=color)
ax2.tick_params(axis='y', labelcolor=color)
#fig.tight_layout()

plt.show()

#The point: PC1 has the widest spread along it (variance, in red). 
#Roughly 50% of the total variation is accounted for by the first 4 PCs

#Read the design data (up until now, have just been working with raw read counts)
design = pd.read_csv('design.csv', index_col = 0)
groups = design['celltype'] + '-' + design['condition']

#Using Seaborn to visualise the data. This is just vsd data labelled using design.csv
import seaborn as sns

sns.scatterplot(x = coordinates[:,0], y = coordinates[:,1], s = 100, 
                hue = groups, palette = sns.color_palette('Paired')[0:10])

plt.legend(bbox_to_anchor = (1,1))
plt.show()




##CLUSTERING. We are no longer using the PCA transformed data
#Use k-means clustering to group the vsd data, labelled using design.csv
from Bio.Cluster import kcluster
clustergraph = clusterid, error, nfound = kcluster(tvsd, dist = 'e', nclusters = 10)

sns.scatterplot(x = coordinates[:,0], y = coordinates[:,1], s = 100, 
                hue = clusterid, palette = sns.color_palette()[0:max(clusterid)+1])

plt.legend(bbox_to_anchor = (1,1))
plt.show()

#A for loop to change the number of clusters
cluster_range = range(2,12)

for n_clusters in cluster_range:
    clusterid, error, nfound = kcluster(tvsd, dist = 'e', nclusters = n_clusters)

    sns.scatterplot(x = coordinates[:,0], y = coordinates[:,1], s = 100, 
                hue = clusterid, palette = sns.color_palette()[0:max(clusterid)+1])
    
    plt.legend(bbox_to_anchor = (1,1))
    plt.show()

#A for loop to produce an Elbow plot
cluster_numbers = []
errors = []

for n_clusters in cluster_range:
    clusterid, error, _ = kcluster(tvsd, nclusters=n_clusters, method='a', dist='e')

    cluster_numbers.append(n_clusters)
    errors.append(error)

sns.lineplot(x=cluster_numbers, y=errors, marker='o')

plt.xlabel('Number of Clusters')
plt.ylabel('Clustering Error')
plt.title('Error vs. Number of Clusters')

plt.show()




##HIERARCHICAL CLUSTERING
import scipy.cluster.hierarchy as sch
Z = sch.linkage(tvsd, method = 'complete')

#Using sch, produce a dendrogram which clusters the labelled vsd data
plt.rcParams['figure.figsize'] = [10,10]
labels = groups.values.tolist()
dend = sch.dendrogram(Z, orientation = 'left', color_threshold = 110, labels = labels)

plt.show()




##DE ANALYSIS
#Now looking at the differences in subsets, eg lung disease vs control.

#Access DE files produced by DEseq2
DElung = pd.read_csv('DElung.csv', index_col = 0) 
DEepi = pd.read_csv('DEepi.csv', index_col = 0) 
DEneu = pd.read_csv('DEneu.csv', index_col = 0) 

#Retrieve row and column labels
DEepi.columns
DEepi.index
DEneu.columns
DEneu.index
DElung.columns
DElung.index

#make a list of genes that have a fold change greater than 2, and significance of 95%
degenesepi = (abs(DEepi['log2FoldChange']) > 2) & (DEepi['padj'] < 0.05)
degenesneu = (abs(DEneu['log2FoldChange']) > 2) & (DEneu['padj'] < 0.05)
degeneslung = (abs(DElung['log2FoldChange']) > 2) & (DElung['padj'] < 0.05)

#use this list to produce a volcano plot (showing DE genes, either up or down regulated. )
plt.rcParams['figure.figsize'] = [6,8]            
plt.scatter(DEepi['log2FoldChange'], -np.log10(DEepi['padj']), 
            s = 3, c = degenesepi, cmap = 'bwr')

plt.xlabel('log2 fold change')
plt.ylabel('-log10 (p-value)')

plt.show()

#use this list to produce a volcano plot (showing DE genes, either up or down regulated. )
plt.rcParams['figure.figsize'] = [6,8]            
plt.scatter(DEneu['log2FoldChange'], -np.log10(DEepi['padj']), 
            s = 3, c = degenesneu, cmap = 'bwr')

plt.xlabel('log2 fold change')
plt.ylabel('-log10 (p-value)')

plt.show()

#use this list to produce a volcano plot (showing DE genes, either up or down regulated. )
plt.rcParams['figure.figsize'] = [6,8]            
plt.scatter(DElung['log2FoldChange'], -np.log10(DEepi['padj']), 
            s = 3, c = degeneslung, cmap = 'bwr')

plt.xlabel('log2 fold change')
plt.ylabel('-log10 (p-value)')

plt.show()

#a function to take one DE dataset, identify DE genes (based on custom fold change and padj), 
#and export up and down regulated genes to text files
def analyse_de(data, downfilename, upfilename, lfc_threshold=2, padj_threshold=0.05):
    
    # select up-regulated genes (use .index to extract the row (gene) names)
    up = (data['log2FoldChange'] > lfc_threshold) & (data['padj'] < padj_threshold)
    up_genes = data.index[up]
    
    down = (data['log2FoldChange'] <- lfc_threshold) & (data['padj'] < padj_threshold)
    down_genes = data.index[down]
    
    print('number of down regulated genes:' , len(down_genes))
    print('number of up regulated genes:' , len(up_genes))
    
    
    with open(upfilename, 'w') as file:
        for gene_name in up_genes:
            file.write(f'{gene_name}\n')
            
    with open(downfilename, 'w') as file:
        for gene_name in down_genes:
            file.write(f'{gene_name}\n')
            
    # Volcano plot
    plt.scatter(data['log2FoldChange'], -np.log10(data['padj']), c='gray', label='Non-DE Genes', alpha=0.5)
    plt.scatter(data.loc[up_genes, 'log2FoldChange'], -np.log10(data.loc[up_genes, 'padj']), c='red', label='Up-regulated Genes')
    plt.scatter(data.loc[down_genes, 'log2FoldChange'], -np.log10(data.loc[down_genes, 'padj']), c='blue', label='Down-regulated Genes')
    plt.xlabel('Log Fold Change')
    plt.ylabel('-log10(P-Value)')
    plt.title('Total Lung')
    plt.axhline(-np.log10(padj_threshold), color='gray', linestyle='--', label='P-Value Threshold')
    plt.axvline(lfc_threshold, color='gray', linestyle='--', label='LogFC Threshold')
    plt.legend()
    
    plt.show()
    
    return len(up_genes), len(down_genes)
    
upfilename = 'DEEPIUP.txt'
downfilename = 'DEEPIDOWN.txt'
data = DEepi
custom_lfc = 2
custom_padj = 0.05

analyse_de(data, downfilename, upfilename, lfc_threshold = custom_lfc, padj_threshold = custom_padj)

upfilename = 'DENEUUP.txt'
downfilename = 'DENEUDOWN.txt'
data = DEneu
custom_lfc = 2
custom_padj = 0.05

analyse_de(data, downfilename, upfilename, lfc_threshold = custom_lfc, padj_threshold = custom_padj)

upfilename = 'DELUNGUP.txt'
downfilename = 'DELUNGDOWN.txt'
data = DElung
custom_lfc = 2
custom_padj = 0.05

analyse_de(data, downfilename, upfilename, lfc_threshold = custom_lfc, padj_threshold = custom_padj)

#this function produces two lists per dataset, up and down regulated genes. These lists can be compared for similarity
#Jaccard similarity test
def jaccard_similarity(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    if union == 0:
        return 0.0  # Avoid division by zero
    
    return intersection / union

#comparing DE genes between epithelia and neutrophils
upfilename1 = 'DENEUUP.txt'
downfilename1 = 'DENEUDOWN.txt'
data1 = DElung

upfilename2 = 'DEEPIUP.txt'
downfilename2 = 'DEEPIDOWN.txt'
data2 = DEneu

# Analyze the first dataset
up_count1, down_count1 = analyse_de(data1, downfilename1, upfilename1, lfc_threshold=custom_lfc, padj_threshold=custom_padj)
up_genes1 = set(pd.read_csv(upfilename1, header=None)[0])
down_genes1 = set(pd.read_csv(downfilename1, header=None)[0])

# Analyze the second dataset
up_count2, down_count2 = analyse_de(data2, downfilename2, upfilename2, lfc_threshold=custom_lfc, padj_threshold=custom_padj)
up_genes2 = set(pd.read_csv(upfilename2, header=None)[0])
down_genes2 = set(pd.read_csv(downfilename2, header=None)[0])

# Compare similarity using Jaccard similarity
up_similarity = jaccard_similarity(up_genes1, up_genes2)
down_similarity = jaccard_similarity(down_genes1, down_genes2)

print(f'Jaccard Similarity for Up-regulated Genes: {up_similarity:.2%}')
print(f'Jaccard Similarity for Down-regulated Genes: {down_similarity:.2%}')




##FINAL STEP: save a reference list of all genes used for enrichment analysis

filename = 'Myreference.txt'

with open(filename, 'w') as file:
    for gene_name in vsd.index:
        file.write(f'{gene_name}\n')

      
