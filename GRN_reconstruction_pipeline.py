#This script loads a single cell RNA/nucleus-seq CellxGene matrix and creates cell type-specifc GRNs as eventual output. In between these steps, we normalize the data, filters cells/genes, run DecoupleR to infer Transcription Factor (TF) activity and define cutoffs for TF/gene activity & expression in a cell type-specifc way. 

# Upload relevant packages 
import scanpy as sc
import pandas as pd
import numpy as np
import decoupler as dc
import os

# Load data.
adata = sc.read("")

# Make the index unique. 
adata.var_names_make_unique()
adata.obs_names_make_unique()
print("adata loaded\n", adata)

        
# Make 'obs' attribute categorical. This can include any obsevrational categroy like cell type and age of sample.
adata.obs['ClusterName'] = adata.obs['ClusterName'].astype('category') 


# Quality Control
        
# Filter cells with high mitochondrial gene expression, then remove mitochondrial, ribosomal and hemoglobin genes.

# Mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('mt-') 
# Ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith('Rp')
# Hemoglobin genes.
adata.var['hb'] = adata.var_names.str.contains(("^Hb.*-"))
        
# Calculate qc metrics for mt, ribo and hb genes. 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)
        
# Filter for percent mito
adata = adata[adata.obs['pct_counts_mt'] < 20, :]
print("adata post pct mito removal\n", adata)

# Remove Malat1 lncRNA, mitochondrial genes, ribosomal genes and blood genes. 
Malat1 = adata.var_names.str.startswith('Malat1')
mito_genes = adata.var_names.str.startswith('mt-')
hb_genes = adata.var_names.str.contains(("^Hb.*-"))
ribo_genes = adata.var_names.str.startswith('Rp')
        
        
remove = np.add(mito_genes, Malat1)
remove = np.add(remove, hb_genes)
remove = np.add(remove, ribo_genes)
keep = np.invert(remove)
        
adata = adata[:,keep]
print("adata post removal of mito, ribo and blood genes\n", adata)

# At this early juncture, we can introduce a 'for loop' for individual analysis across age categories or bypass the loop to analyze cells/cell-types across all time periods together. This script will proceed through the whole brain analysis. In order to perform an age delineated analysis, the following can be adopted, for example: 

# Define the categories you want to loop over (assuming they are in the 'Age' column)
#categories = adata.obs['Age'].unique()

#print("final categories to loop over\n", categories)

# Loop over each category and perform the analysis. 
#for category in categories:
    #print(f"Processing category: {category}")
    
    # Subset the data to include only cells belonging to this category
    #adata_cat = adata[adata.obs['Age'] == category, :]
#... This can then proceed by incorperating the entirety of the remaining script into the loop. The only difference being that instances of 'adata' should be changed to 'adata_cat', i.e. the current subset. 

# Most processed datasets will still have cell types labeled undefined or loq quality, these should be removed.
adata = adata[adata.obs['ClusterName'] != 'Undefined', :] # Substitute Undefined
print("adata after filtering undefined\n", adata)


# Remove cell types with less than 5 cells
adata = adata[adata.obs.groupby('ClusterName').ClusterName.transform(len) >= 5, :]
print("adata after removing cell-types < 5 cells\n", adata)        
print(adata)
        
# Remove cells with less than 200 genes, along with genes detected in less than 5 cells. 
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=5)
print("adata after filtering for min cells/genes\n, adata")


# Normalize the data
sc.pp.normalize_total(adata)

# Log-transform the data
sc.pp.log1p(adata)

# (Optional) Store the cell types and their corresponding number of cells
celltype_counts = adata.obs['ClusterName'].value_counts()
df = pd.DataFrame(celltype_counts)
# Sort the dataframe in descending order by the count column
df = df.sort_values(by=['ClusterName'], ascending=False)
# Save the dataframe to a CSV file in the specified directory
os.chdir( "")
df.to_csv("Celltype_counts.csv")
       
# Second Major step: Leverage the DecoupleR Multivariate Linear Model (MLM) against the adata matrix to infer TF activities.

# Aqcuire source-target DoRothEA network, lvls 1-3 of confidence (ABC). 
net = dc.get_dorothea(organism='Mouse', levels=['A','B','C']) # Specify mouse/human

# Inititlize MLM where the covariates are defined as the 'source'-'target' entries in the DoRothEA network. By default DecoupleR will only include TFs with > 5 targets for analysis. 
dc.run_mlm(mat=adata, net=net, use_raw=False, source='source', target='target', weight='weight', verbose=True)

#Extract activities from DecoupleR output, this output is an adata object. .var contains the TFs, and .obs contains cells, .X atttribute contains the TF activities. 
acts = dc.get_acts(adata, obsm_key='mlm_estimate')
        
# The MLM can infer active TFs with no read in the matrix, these will be filtered. 
# Define all TFs identified by the MLM
sources_mlm = adata.obsm['mlm_estimate'].columns.astype(str).values.tolist()
# Define set of all trascripts in the original data
Allgenes = adata.var_names.astype(str)
# Keep TFs with a read in the original matrix
source_in_raw = set(Allgenes).intersection(set(sources_mlm))
source_in_raw = list(source_in_raw)
# Subset `acts matrix` to include only genes in `source_in_raw`, i.e TFs with a read in original adata matrix.
acts = acts[:, source_in_raw]

# Filter cells with no TF expression and TFs expressed in less than 20 cells.
sc.pp.filter_genes(acts, min_cells=20)
sc.pp.filter_cells(acts, min_genes=1)
print("acts after filtering"\n, acts)

# (Optional) save to inspect
#os.chdir("")
#acts.write(".h5ad")


# Summarise activities based on mean expression per celltype
mean_acts = dc.summarize_acts(acts, groupby='ClusterName', min_std=0)
print(mean_acts)

# Delete clutter
del acts

# Calculate the 25th percentile of TF activity based on non-zero values for each celltype
def percentile_nonzero(x):
    x_nonzero = x[x > 0]
    if len(x_nonzero) > 0:
        return x_nonzero.quantile(q=0.25)
    else:
        return 0
           
percentiles = mean_acts.apply(percentile_nonzero, axis=1)
print(percentiles)

# Loop over each row (celltype) and set any expression values below the 25th percentile to 0
threshold = []
for celltype in mean_acts.index:
    threshold = percentiles[celltype]
    mean_acts.loc[celltype] = mean_acts.loc[celltype].apply(lambda x: x if (x >= threshold and x >= 0) else 0)

# Remove any columns with a sum of 0 (These TFs are filtered)
mean_acts = mean_acts.loc[:, (mean_acts != 0).sum() > 0]

# We have a new list of Regulons which passed activity QC
RegsinRaw = mean_acts.columns.tolist()

# Slight tangent: Now we want filter the DecoupleR network table based on targets/TFs which are present in adata and passed network QC. We also remove negative (inhibitory interactions) from the DoRothEA network. Skip forward for TF expression QC. 
# First we filter the target column.
# We remove targets in the DecoupleR network not present in our data.
index_list_1 = []
for target, entry in enumerate(net.target, start=0):
    if entry in list(adata.var.index.values):
        index_list_1.append(target)
        
net_fix_target = net.iloc[index_list_1]

# Now we filter the source column, using target corrected table 
index_list_2 = []
        
for source, entry in enumerate(net_fix_target.source, start=0):
    if entry in list(RegsinRaw):
        index_list_2.append(source)
        
net_fixed = net_fix_target.iloc[index_list_2]
#Now both source and targets entries exist in our adata matrix.

# Remove rows where weight is negative, we only want positive interactions.
net_fixed = net_fixed[net_fixed['weight'] >= 0]
# Optional save the dataframe to the specified directory

        
# Continuing TF activity QC
# Create boolean True/False of each TF activity (if it passed percentile cutoff) for each obs category.   
Raw_acts_regs = mean_acts > 0

# Create dictionary for each celltype with true/flase value for each TF based on above activity cutoff (25% celltype-specifc).Format: Key=celltype, entry='TF:True/False'  
Activities_boolean_dict = {}
        
for cluster, entry in zip(list(Raw_acts_regs.index),Raw_acts_regs.iloc):
    temp_dict2 = dict(zip(entry.index,list(entry)))
    Activities_boolean_dict[cluster]=temp_dict2

#To perform cutoff of TF expression in raw adata, subset these TFs. 
subset_matrix = adata[:,RegsinRaw]
clusters = subset_matrix.obs['ClusterName'].cat.categories
exprarray = subset_matrix[:,RegsinRaw].X.toarray()
Raw_expr_regs = pd.DataFrame(exprarray,columns=RegsinRaw,index=adata.obs['ClusterName'])
            
            
# Group the data by obs attribute and calculate the mean for each TF in each cell type
Raw_expr_regs = Raw_expr_regs.groupby('ClusterName').mean()
            
# Calculate the 25th percentile of TF expression based on non-zero values
def percentile_nonzero(x):
    x_nonzero = x[x > 0]
    if len(x_nonzero) > 0:
        return x_nonzero.quantile(q=0.25)
    else:
        return 0
        
percentiles = Raw_expr_regs.apply(percentile_nonzero, axis=1)
    
# Remove any TF below the specified expression threshold for each celltype.   
threshold = []
for celltype in Raw_expr_regs.index:
    threshold = percentiles[celltype]
    Raw_expr_regs.loc[celltype] = Raw_expr_regs.loc[celltype].apply(lambda x: x if x >= threshold else 0)
            
        
#create boolean true/false for presence of TF in each celtype
Reg_expr_fraction_50 = Raw_expr_regs > 0
    
# Create dictionary for each celltype with true/flase value for each TF based on above expression cutoff. Format: Key=celltype, entry='TF:True/False'  
Regsinraw_expr_dict = {}
            
for cluster, entry in zip(list(Reg_expr_fraction_50.index),Reg_expr_fraction_50.iloc):
    temp_dict = dict(zip(entry.index,list(entry)))
    Regsinraw_expr_dict[cluster]=temp_dict
    
# Delete clutter    
del Reg_expr_fraction_50
del exprarray
del subset_matrix
del Raw_expr_regs
del clusters
del Raw_acts_regs

# Convert both dictionaries to dataframe to perform AND operation for source activity + expression.
Regsinraw_expr_df = pd.DataFrame.from_dict(Regsinraw_expr_dict)
Activities_boolean_df = pd.DataFrame.from_dict(Activities_boolean_dict)
    
# Make sure indices align before performing AND operation
align_columns = Regsinraw_expr_df.columns.intersection(Activities_boolean_df.columns)
Regsinraw_expr_df = Regsinraw_expr_df.reindex(columns=align_columns)
Activities_boolean_df = Activities_boolean_df.reindex(columns=align_columns)
    
# Perform AND operation.
Filtered_Regs = np.logical_and(Regsinraw_expr_df, Activities_boolean_df)
    
# Transopse for orientation (determines keys in next step, we will conver this back to dictionary)
Filtered_Regs = Filtered_Regs.T
#convert to dict
Filtered_Regs = Filtered_Regs.to_dict(orient = 'index')
    
# We have now created one dictionary 'Filtered_Regs', which has format key=celltype, entry='TF:True/False'.
# If 'TF:True', that TF has activity >25th percentile and expression >25th percentile in that celltype. 

# Now for each TF in each celltype, if TF = True, append row of decoupler network table where source = TF to the dictionary. 
# This means we have a dictionary 'dc_dict' with format key = celltype, entry = network dataframe. 
dc_dict = {}
for cluster in Filtered_Regs.keys():
    keep_index=[]
    for num, source in enumerate(list(net_fixed.source),start=0):
        boolean_df = pd.DataFrame.from_dict(Filtered_Regs[cluster], orient='index', columns=['bool'])
        True_bool_index_list = boolean_df[boolean_df['bool']==True].index
        if source in True_bool_index_list:
            keep_index.append(num)
    net_fixed1 = net_fixed.iloc[keep_index]
    dc_dict[cluster] = net_fixed1

# The network dataframe contains interactions where the 'source' passed TF QC. However, the targets are not corrected in a celltype specifc way.
# Perform target raw expression cutoff
# To perform cutoff of gene expression in raw adata, subset these genes. 
subset_genes_tar = adata[:,adata.var.index]
clusters_tar = subset_genes_tar.obs['ClusterName'].cat.categories
genesarray_tar = subset_genes_tar[:,adata.var.index].X.toarray()
Raw_gene_exp_tar = pd.DataFrame(genesarray_tar,columns=adata.var.index,index=adata.obs['ClusterName'])
        
# Group the data by obs and calculate the mean for each gene in each cell type
Raw_gene_exp_tar = Raw_gene_exp_tar.groupby('ClusterName').mean()
        
# Calculate 25th percentile of gene expression for each celltype based on non-zero values
def percentile_nonzero(x):
    x_nonzero = x[x > 0]
    if len(x_nonzero) > 0:
        return x_nonzero.quantile(q=0.25)
    else:
        return 0
    
percentiles = Raw_gene_exp_tar.apply(percentile_nonzero, axis=1)

# Loop over each row (celltype) and set any expression values below the 25th percentile to 0
threshold = []
for celltype in Raw_gene_exp_tar.index:
    threshold = percentiles[celltype]
    Raw_gene_exp_tar.loc[celltype] = Raw_gene_exp_tar.loc[celltype].apply(lambda x: x if x >= threshold else 0)
        
# Remove any columns with a sum of 0
Raw_gene_exp_tar = Raw_gene_exp_tar.loc[:, (Raw_gene_exp_tar != 0).sum() > 0] 

# Create boolean table to create dictionary (same as source)
Raw_gene_exp_cutoff_tar = Raw_gene_exp_tar > 0

        
# Create target expression dict. keys = cell type, entry = gene:True/False
Target_expression_dict = {}
        
for cluster, entry in zip(list(Raw_gene_exp_cutoff_tar.index),Raw_gene_exp_cutoff_tar.iloc):
    temp_dict3 = dict(zip(entry.index,list(entry)))
    Target_expression_dict[cluster]=temp_dict3

#delete clutter
del subset_genes_tar
del clusters_tar
del Raw_gene_exp_cutoff_tar
del Raw_gene_exp_tar
del genesarray_tar
del adata

# Now we make use of our Filtered_Regs dictionary and Target_expression dict. 
# Reminder: Filtered_Regs dict contains key =celltype, entry = dataframe. 
# The 'source'(TFs) in each celltype have passed the celltype specifc TF QC (act/expr).
# However, perhaps some of the targets for these sources have not passed target expression QC celltype-wide.
# Therefore, we remove rows from each network in each key present in 'Filtered_Regs', if the 'target' gene failed Target expression QC. 
        
# Create the final dictionary, key = celltype, entry = network dataframe. 
result = {}
        
# Loop through each key in dc_dict
for key in dc_dict.keys():
  # Extract entry (a table in this case) based on the current key in dc_dict
    network = pd.DataFrame.from_dict(dc_dict[key])
            
    # Extract "bool" column from Target_expression_dict
    boolean1 = pd.DataFrame.from_dict(Target_expression_dict[key], orient='index', columns=['bool'])
    # Define a list of true values
    True_bool1 = boolean1[boolean1['bool']==True].index
    # Remove row from network table if target in row is not present in True_bool1
    networkcorr = network[network['target'].isin(True_bool1)]
            
    # Store the resulting DataFrame in the result dictionary with the key being the name of the obs attribute
    result[key] = networkcorr

# Remove source from each network if the TF has less than 5 interactions in that celltype.
# This correction is applied since DecoupleR calculates TF activity scores for TFs >5 targets. 
# Loop through the dictionary
for key, network_df in result.items():
    # Get the count of each source value in the dataframe
    source_counts = network_df['source'].value_counts()

    # Filter the dataframe to only keep rows where the source value appears at least 5 times
    filtered_network_df = network_df[network_df['source'].isin(source_counts[source_counts >= 5].index)]

    # Replace the original network dataframe with the filtered one
    result[key] = filtered_network_df
    
# Export network tables

# Specify the directory where you want to save the network tables
directory = os.path.join("")

if not os.path.exists(directory):
    os.makedirs(directory)
        
# Loop through each dataframe in the result dictionary
for key, df in result.items():
    # Construct the filename based on its key
    output_filename = key + '.csv'
    # Join the directory and filename to get the full path
    filepath = os.path.join(directory, output_filename)
    # export the dataframe to a CSV file
    df.to_csv(filepath, index=False)