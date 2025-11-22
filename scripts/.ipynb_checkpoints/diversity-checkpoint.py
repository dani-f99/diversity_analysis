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
        aa_temp = df_aas[i]
        richness = len(aa_temp.unique())
        freq = aa_temp.value_counts()/aa_temp.value_counts().sum()
        freq_table = freq.reset_index(name="freq").rename({"index":"aa"},axis=1)
        freq_table.index = freq_table.index + 1
        
        frequencies = freq_table["freq"].values
        
        import math
        t = 0
        
        for k in frequencies:
            
            Pi = k * math.log(k, math.e)
            t += Pi
        
        diversity = round(math.exp(-t),3)
        
        temp_df.loc[temp_df["pos_aa"]==i,["richness","diversity"]] = [richness,diversity]
    
    return temp_df