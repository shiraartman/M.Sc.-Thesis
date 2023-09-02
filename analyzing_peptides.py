import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter


def get_best_peptides (file):
    """ This function loads a table from a csv file, calculates enrichment of each peptide, takes the 100 most enriched peptides
    and returns a dataframe containing the 100 most enriched peptides.
    - parameters:
             file: a string of the path to the table to be read.
    - return: 
             peptides: a dataframe of the 100 most enriched peptides. """
             
    data=pd.read_csv(file)
    o1=data["O1_reads"]
    o3=data["O3_reads"]
    data["Enrichment"]=(o3/o1)*100
    # to get the % of enrichment
    sorted_data=data.sort_values("Enrichment", ascending=False)
    # sort by % of enrichment
    top_100=sorted_data[0:101]
    #get the top 100 peptides
    peptides=top_100["peptide_AA"]
    # extract the peptides
    return peptides
    

def aa_positions (peptides):
    """ This function receives a dataframe of peptides and returns a dataframe with the amino acids of each peptide and
    their position in the peptide.
    - parameters:
             peptides: a dataframe containing peptides.
    - return: 
             aa_and_pos: a new dataframe containing amino acids from each peptide and their position in it. """

    aa=[]
    positions=[]
    aa_and_pos=pd.DataFrame(columns=["Amino Acid","Position"])
    for peptide in peptides:
        for aa_id in range(0,len(peptide)):
            positions.append(aa_id)
            aa.append(peptide[aa_id])
    positions=np.add(positions, 1).tolist()
    # so the index of the first position will be 1 instead of 0
    aa_and_pos["Amino Acid"]=aa
    aa_and_pos["Position"]=positions
    return aa_and_pos


def unique_pairs_counter (df, column_name1, column_name2):
    """ This function receives a dataframe and the names of two columns in it, counts unique pairs of the desired columns
    (e.g., amino acids and positions) and adds the counts to a new column in the given dataframe.
    - parameters:
             column_name1: a string of the name of the first desired column.
             column_name2: a string of the name of the second desired column.
             df: the dataframe in which the two columns are.
    - return: 
             unique_combos: the same given dataframe, but with a new column of the counts of each unique pair of both
             specified columns. """
    
    unique_combos=df.groupby([column_name1,column_name2]).size().reset_index().rename(columns={0:"Count"})
    return unique_combos
    

def get_aa_type (aa):
   """ This function gets an amino acid and returns its type.
   - parameters:
            aa: a string of the one letter coding of the amino acid.
   - return: 
            type_of_aa: a string of the type of the amino acid. """
            
   aa_dict={"Non-polar":("G","A","V","L","M","I"), "Aromatic":("F","Y","W"), "Polar":("P","S","N","T","C","Q"),
             "Positive":("K","R","H"), "Negative":("D","E"), "Stop":"*"}
   for type_of_aa, value in aa_dict.items():
         if aa in value:
             return type_of_aa
         

def type_and_position (aa_and_pos):
    """ This function gets a dataframe containing amino acids from each peptide and their position in it, and returns
    the same dataframe but with an added column of the type of each amino acid.
    - parameters:
             aa_and_pos: a dataframe containing amino acids from each peptide and their position in it.
    - return: 
             aa_and_pos: the same dataframe, with an added column of the type of each amino acid. """
             
    types=[]
    for a_acid in aa_and_pos["Amino Acid"]:
       types.append(get_aa_type(a_acid))
    aa_and_pos["Type of Amino Acid"]=types
    return aa_and_pos  
        

def aa_types_freq (peptides):
    """ This function receives a dataframe of peptides, and returns a dictionary of types of amino acids as keys and the count of each
    type among all the peptides as values.
    - parameters:
             peptides: a dataframe containing peptides.
    - return: 
            sorted_aa_types: a dictionary of types of amino acids as keys and their count among all peptides as values. """
    amino_acids=[]
    aa_types=[]
    for peptide in peptides:
        for aa in peptide:
            amino_acids.append(aa)
    for amino_acid in amino_acids:
        aa_types.append(get_aa_type(amino_acid))
    aa_types_counter=Counter(aa_types)
    return aa_types_counter
    

def main():
     peptides=get_best_peptides("PSI_peptides.csv")
     
     
     aa_vs_pos=aa_positions(peptides)
     aa_vs_pos_counter=unique_pairs_counter(aa_vs_pos, "Amino Acid", "Position")
    
     # heatmap 1 - amino acids as funtion of position
     plt.figure()
     pivoted1=aa_vs_pos_counter.pivot(index="Amino Acid", columns="Position", values="Count")
     sns.heatmap(pivoted1[:][::-1], xticklabels=True, yticklabels=True, cmap="Greens")
     plt.xlabel("Position", size=16)
     plt.ylabel("Amino Acid", rotation="vertical", size=16)
     plt.title("Amino Acids as Function of Position")
     plt.tight_layout()
     plt.savefig("Amino Acids as Function of Position.png", dpi=300)
     
    
     type_vs_pos=type_and_position(aa_vs_pos)
     type_vs_pos_counter=unique_pairs_counter(type_vs_pos, "Type of Amino Acid", "Position")
     
     # heatmap 2 - types of amino acids as function of position
     plt.figure()
     pivoted2=type_vs_pos_counter.pivot(index="Type of Amino Acid", columns="Position", values="Count")
     sns.heatmap(pivoted2[:][::-1], xticklabels=True, yticklabels=True, cmap="Blues", square=True, cbar_kws={'shrink': 0.45})
     #square=True and cbar_kws were added to have good proportions between the heatmap and the colorbar
     plt.xlabel("Position", size=12)
     plt.ylabel("Type of Amino Acid",rotation="vertical", size=12)
     plt.title("Types of Amino Acids as Function of Position")
     plt.tight_layout()
     plt.savefig("Types of Amino Acids as Function of Position.png", dpi=300)
     
     
     freq_of_aa_type=aa_types_freq(peptides)
     
     # bar plot - counts of amino acids types among all peptides
     plt.figure()
     types=freq_of_aa_type.keys()
     counts=freq_of_aa_type.values()
     plt.bar(types, counts, color="pink")
     plt.xlabel("Type of Amino Acid",size=12)
     plt.ylabel("Number of Occurrences",size=12)
     plt.title("Counts of Amino Acids Types")
     plt.tight_layout()
     sns.despine
     plt.savefig("Counts of Amino Acids Types.png", dpi=300)
     
     
     plt.show()

if __name__ == "__main__":
     main()







    
            
    



    