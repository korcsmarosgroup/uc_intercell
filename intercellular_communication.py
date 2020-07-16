import pandas as pd
from itertools import permutations


# importing expression information from healthy and diseased conditions in different cells as pandas dataframes
df_imm = pd.read_csv('resources/filtered_imm_log2_TPM.csv', index_col = 0, header=0)
df_epi = pd.read_csv('resources/filtered_epi_log2_TPM.csv', index_col = 0, header=0)
df_fib = pd.read_csv('resources/filtered_fib_log2_TPM.csv', index_col = 0, header=0)


# select healthy/uninflamed average expression data in the 5 selected cell types
df_macrophage_healthy = df_imm['RNA.Macrophages_Healthy']
df_macrophage_uninflamed = df_imm['RNA.Macrophages_Uninflamed']

df_DC1_healthy = df_imm['RNA.DC1_Healthy']
df_DC1_uninflamed = df_imm['RNA.DC1_Uninflamed']

df_goblet_healthy = df_epi['RNA.Goblet_Healthy']
df_goblet_uninflamed = df_epi['RNA.Goblet_Uninflamed']

df_myofibroblast_healthy = df_fib['RNA.Myofibroblasts_Healthy']
df_myofibroblast_uninflamed = df_fib['RNA.Myofibroblasts_Uninflamed']

df_Treg_healthy = df_imm['RNA.Tregs_Healthy']
df_Treg_uninflamed = df_imm['RNA.Tregs_Uninflamed']


# dictionary to store dataframes and make it callable in a for cycle
healthy_cell_types = {"DC1": df_DC1_healthy,"Macrophage": df_macrophage_healthy, "Goblet": df_goblet_healthy,
                      "Myofibroblast": df_myofibroblast_healthy, 'Treg': df_Treg_healthy}

uninflamed_cell_types = {"DC1": df_DC1_uninflamed,"Macrophage": df_macrophage_uninflamed, "Goblet": df_goblet_uninflamed,
                         "Myofibroblast": df_myofibroblast_uninflamed, 'Treg': df_Treg_uninflamed}


# create dictionaries for source-target interactions and their annotations from OmniPath
interactions = {}
interaction_annotation = {}

with open ('resources/new_version/intercellular_interactions.tsv') as interaction:
    interaction.readline()
    for line in interaction:
        line = line.strip().split('\t')
        if line[7] not in interactions:
            interactions[line[7]] = []
        interactions[line[7]].append(line[6])

        if (line[7],line[0]) not in interaction_annotation:
            interaction_annotation[(line[7], line[0])] = []
        interaction_annotation[(line[7],line[0])].append((line[6],line[4]))


# create all of the possible interactions of cells (independently from the condition)
# important: directionality (A-B and B-A are different)
cell_cell_connections = permutations(healthy_cell_types.keys(), 2)


# going through each cell pair and select interactions between them where one part is represented in one condition only
for i in list(cell_cell_connections):

     # storing cell-type specific connections in sets
     healthy_interactions = set()
     UC_interactions = set()

     # iterating through source cell type expressed genes
     for source in healthy_cell_types[i[0]].keys():

        # selecting those ones which play role in intercellular communication as a transmitter
        if source in interactions.keys():

            # iterating through target cell type expressed genes
            for target in healthy_cell_types[i[1]].keys():

                # selecting those ones which play role in intercellular communication as a receiver
                if target in interactions[source]:
                    if str(healthy_cell_types[i[1]][target]) != 'nan' and str(healthy_cell_types[i[0]][source]) != 'nan':
                        healthy_interactions.add((source, target))

     # same method for uninflamed condition
     for source in uninflamed_cell_types[i[0]].keys():
         if source in interactions.keys():
            for target in uninflamed_cell_types[i[1]].keys():
                if target in interactions[source]:
                    if str(uninflamed_cell_types[i[1]][target]) != 'nan' and str(healthy_cell_types[i[0]][source]) != 'nan':
                        UC_interactions.add((source, target))


     # writing out the condition specific interactions
     with open(i[0] + "_" + i[1] + "_healthy_only.txt", 'w') as output_file_1:
         with open(i[0] + "_" + i[1] + "_UC_only.txt", 'w') as output_file_2:
                healthy_only = healthy_interactions.difference(UC_interactions)
                UC_only = UC_interactions.difference(healthy_interactions)
                print(UC_only)
                for i in healthy_only:
                    for source_annotation in interaction_annotation:
                        if i[0] == source_annotation[0]:
                            for target_annotation in interaction_annotation[source_annotation]:
                                if i[1] == target_annotation[0]:
                                    output_file_1.write(
                                        source_annotation[0] + "," + source_annotation[1] + "," + target_annotation[0]
                                        + "," + target_annotation[1] + "\n")

                for j in UC_only:
                    for source_annotation in interaction_annotation:
                        if j[0] == source_annotation[0]:
                            for target_annotation in interaction_annotation[source_annotation]:
                                if j[1] == target_annotation[0]:
                                    output_file_2.write(
                                        source_annotation[0] + "," + source_annotation[1] + "," + target_annotation[0]
                                        + "," + target_annotation[1] + "\n")
