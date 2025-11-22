import pandas as pd
import numpy as np
import math

#
divrs = lambda seires : 0 if len(seires) == 0 else math.exp(seires.apply(lambda x : -(x*math.log(x,math.e))).sum())

#
serine_codons = {"AGT":'S"',
                 "AGC":'S"',
                 "TCT":"S'",
                 "TCC":"S'",
                 "TCA":"S'",
                 "TCG":"S'"}  

#
def ser_tpye(df:pd.DataFrame, #mut_df input
             target: str): #from_aa or to_aas
      
    df_input = df[["pos_aa","germline","top_seq","from_aa","to_aas"]]
    aa_from = df_input["from_aa"] 
    aa_to = df_input["to_aas"]
    aa_pos = int(df_input["pos_aa"])
    codon_germ = df["germline"][aa_pos*3-3:aa_pos*3]
    codon_top = df["top_seq"][aa_pos*3-3:aa_pos*3]

    if target == "from_aa": 
        try:        
            serine = serine_codons[codon_germ]
        except:
            serine = "S*"
     
    elif target == "to_aas":
        try:
            serine = serine_codons[codon_top]            
        except:
            serine = "S*"

    return serine

# looking for correct serine codon in the sequences dataframe
def get_ser_sequence(ser_df : pd.DataFrame,
                     seq_df : pd.DataFrame):

    serine_codons = {"AGT":'S"',
                     "AGC":'S"',
                     "TCT":"S'",
                     "TCC":"S'",
                     "TCA":"S'",
                     "TCG":"S'"}  
  
    clone_id = ser_df["clone_id"]
    pos_aa = int(ser_df["pos_aa"])
    seq_list = seq_df.loc[(seq_df["clone_id"] == clone_id),"sequence"].unique()

    for sq in seq_list:
        codon = sq[pos_aa*3-3:pos_aa*3]
        try:
            ser_final = serine_codons[codon]
            break
        except:
            ser_final = "S*"
     
    return ser_final

codon_dic_updated = {
                    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
                    'TCT': "S'", 'TCC': "S'", 'TCA': "S'", 'TCG': "S'",
                    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',  # * for STOP
                    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',

                    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
                    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

                    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
                    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                    'AGT': 'S"', 'AGC': 'S"', 'AGA': 'R', 'AGG': 'R',

                    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
                    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
                    }

def nt_transalte_104(nt_seq):
    import numpy as np
    translated = []
    for i in range(1,105):
        codon = nt_seq[i*3-3:i*3]
      
        if codon in list(codon_dic_updated.keys()):
            aa = codon_dic_updated[codon]  
        else:
            aa = np.nan
            
        translated.append(aa)
    return translated


def diversity(df):
    
    df_grouped = df.groupby("clone_id").agg({"germline":"unique"})
    df_grouped.reset_index(inplace=True)
    df_grouped["germline"]=df_grouped["germline"].apply(lambda x: x[0])
    df_grouped.rename({"germline":"germline_nt"},axis=1,inplace=True)
    df_grouped["germline_aa"] = df_grouped["germline_nt"].apply(nt_transalte_104)
    aas_list = df_grouped["germline_aa"].to_list()
    df_aas = pd.DataFrame(aas_list,columns=range(1,105))
    #df_aas.dropna(inplace=True)
    
    # for each position need:
    # 1. richness
    # 2. frequencies per aa
    # 3. diversity calculation
    
    div_list = []
    aa_range = list(range(1,105))
    
    import numpy as np
    temp_df = pd.DataFrame(0, index=np.arange(104), columns=["pos_aa","richness","diversity"])
    temp_df["pos_aa"] = aa_range
    
    for i in aa_range:
        aa_temp = df_aas.loc[:,i][df_aas.loc[:,i].notnull()]
        richness = len(aa_temp.unique())

        if richness == 0:         
            richness = 0
            diversity = 0
        
        else:
            freq = aa_temp.value_counts()/aa_temp.value_counts().sum()
            freq_table = freq.reset_index(name="freq").rename({"index":"aa"},axis=1)
            freq_table.index = freq_table.index + 1
            
            frequencies = freq_table["freq"].values
            
            import math
            t = 0
            
            for k in frequencies:
                
                Pi = k * math.log(k, math.e)
                t += Pi
            
            diversity = round(math.exp(-t))

        temp_df.loc[temp_df["pos_aa"]==i,["richness","diversity"]] = [richness,diversity]
    
    return temp_df

def getdiv(input_df:pd.DataFrame) -> int:
    return round(divrs(input_df.value_counts()/input_df.value_counts().sum()))