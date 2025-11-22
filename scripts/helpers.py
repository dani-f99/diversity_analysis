import matplotlib.pyplot as plt
from scipy.stats import wilcoxon
from matplotlib import cm
import pandas as pd
import numpy as np
import math

# helper function used for diversity calculation
diversity = lambda seires : 0 if len(seires) == 0 else math.exp(seires.apply(lambda x : -(x*math.log(x,math.e))).sum())

# helper function to find unmathcing serine codons 1 distance away from codon
def find_close_ser(input: str,
                   codons: list):
    for i in codons.keys():
        len_m = sum(a==b for a, b in zip(input , i))
        codon = i
        
        if len_m == 2:
            return codons[codon]
    else:
        return "S*"
    
''' without "interpolation" '''
''' extracing serine value function '''

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

# aa order for graph x axis
aa_order = ["R","E","K","Q","D",
            "N",'S"',"G","W","C",
            "H","Y","A","T","P",
            "S'","V","M","I","L",
            "F"]

#serine variants
#ser_codon = {'S"':["AGT","AGC"],
#             "S'":["TCT","TCC","TCA","TCG"]}
serine_codons = {"AGT":'S"',
                 "AGC":'S"',
                 "TCT":"S'",
                 "TCC":"S'",
                 "TCA":"S'",
                 "TCG":"S'"}   

# custom function to round numbers upward
def round_up(number):
    num_dec = number
    num_round = round(number)
    
    if num_round < num_dec:
        value = num_round + 1
    else:
        value = num_round
    return value
    
# Creating matrix with maximum of 5 columns
def penta_matrix(n_cols):
    if round_up(n_cols/5) > 1:
        shape = (round_up(n_cols/5), 5)
    else:
        shape = (1,n_cols)

    shape_itr = [(i,j) for i in range(0,shape[0]) for j in range(0,shape[1])]

    return shape, shape_itr[:n_cols]

# frequnceis function
def analysis_information(input_df : pd.DataFrame, 
                         name_metric : str,
                         split_by : str = None,
                         pos_aa_range : range = None) -> [dict, pd.DataFrame, pd.DataFrame]:
    
    """
    A function that accepts mut_df format dataframe and returns dictionary of frequnceies (key=condition)
    where:
    1. input_df - mut_df dataframe.
    2. split_by - column name that will be used to filter the dataframe by its unique values. 
    3. pos_aa_range - amino acids position range (make sure the pos_aa is int and not float/str)
    """

    if name_metric not in ["div","rich"]:
        raise ValueError('{} not rich (richness) or div (diversity), enter correct value.'.format(name_metric))

    # Creating the template dataframe for the amino acid frequencies
    pos_aa_range = range(input_df.pos_aa.astype("int").min(),input_df.pos_aa.astype("int").max()+1)
    df = input_df[input_df["pos_aa"].isin(pos_aa_range)]
    freq_df_template = pd.DataFrame(data=0, index = pos_aa_range, columns=aa_order).astype(float)
    freq_dics = {}

    # Creating filtring condition for the output
    if split_by is None:
        cond_tempdf = ["from_aa.", "to_aas."]
    else:
        cond_tempdf = np.sort([".".join([by, col]) for col in ["from_aa","to_aas"] for by in df[split_by].unique()])

    # Creating diversity and richness dataframe to be filled later
    metric_df_template = pd.DataFrame(index = pos_aa_range, columns = cond_tempdf)
    metric_df = metric_df_template.copy()
    nclones_df = metric_df_template.copy()
    
    # Itirating over the conditions
    for c in cond_tempdf:
        c_split = list(filter(None,c.split(".")))
        temp_freq = freq_df_template.copy()
        
        # determining the input dataframe (according to the filtring) and output name.
        if len(c_split) == 1:
            temp_df = df
            name = c[:-1]
            split_i = 0

        elif len(c_split) == 2:
            temp_df = df[df[split_by]==c_split[0]]
            name = c
            split_i = 1

        # Getting the information of the frequncies, diversity and richness and assiging it to the output.
        for i in pos_aa_range:
            nclones_df.loc[nclones_df.index == i, c] = len(temp_df.loc[temp_df["pos_aa"] == i, "clone_id"].unique())

            if name_metric == "rich":
                frequencies = temp_df.loc[temp_df["pos_aa"] == i, c_split[split_i]].value_counts().to_frame()
                fractions = round(frequencies/frequencies.sum(),3)
                temp_freq.loc[temp_freq.index == i, fractions.T.columns] = fractions.T.values
                metric_df.loc[metric_df.index == i, c] = frequencies.shape[0]
            
            elif name_metric == "div":
                # finidng diversity value
                frequencies_og = temp_df.loc[temp_df["pos_aa"] == i, c_split[split_i]].value_counts().to_frame()
                diversity_val = round(diversity(frequencies_og.iloc[:,0]/frequencies_og.iloc[:,0].sum()))
                metric_df.loc[metric_df.index == i, c] = diversity_val
                
                # assiging the diversity values (rounded) to the datasets.
                frequencies_div = temp_df.loc[temp_df["pos_aa"] == i, c_split[split_i]].value_counts().to_frame().iloc[:diversity_val,]
                fractions_div = round(frequencies_div/frequencies_div.sum(),3)
                temp_freq.loc[temp_freq.index == i, fractions_div.T.columns] = fractions_div.T.values
           
        freq_dics[name] = temp_freq

    return freq_dics, metric_df, nclones_df

# Plotting bars which represents the frequencies of amino acids according to the
# analysis_information output.
def bar_freqs(input_freqs : dict,
              input_metric : pd.DataFrame,
              nclones_df: pd.DataFrame,
              name_metric : str,
              title : str,
              pos_plot : list = None,
              save_fig : bool = False,
              change_xticks : bool = False) -> plt.Figure:
    
    """
    input_freqs : dict -> 'freq_dics' output of the analysis_information function.
    input_metric : pd.DataFrame -> 'metric_df' output of the analysis_information function.
    nclones_df: pd.DataFrame -> 'nclones_df' output of the analysis_information function.
    name_metric : str -> Name of the metric used during the analysis (diversity usally).
    title : str -> title of the figure.
    pos_plot : list = None -> Determining the needed aa for the plot & creating parameters for figure
    save_fig : bool -> if True the figure will be saved into the output folder
    change_xticks : bool -> change the naming scheme of the x-ticks.
    """

    plt.rcParams['font.size'] = 20

    if name_metric not in ["div","rich"]:
        raise ValueError('{} not rich (richness) or div (diversity), enter correct value.'.format(name_metric))

    # Colors for the barplot legend
    n_colors = len(aa_order)
    colors = cm.tab20(np.linspace(0, 1, n_colors))
    colors = list(colors)
    colors[0] = "blue"
    colors[-1] = "lavender"

    # Determining the needed aa for the plot & creating parameters for figure
    if pos_plot is None:
        pos_plot = input_metric.index.values
    else:
        pos_plot = pos_plot

    n_pos = len(pos_plot)
    matrix_shape, matrix_itr = penta_matrix(n_pos)
    all_range = [(i,j) for i in range(0,matrix_shape[0]) for j in range(0,matrix_shape[1])]
    noaxis = [i for i in all_range if i not in matrix_itr]

    # Initiating the figure
    fig, axs = plt.subplots(matrix_shape[0],
                            matrix_shape[1], 
                            figsize=(8*matrix_shape[1],10*matrix_shape[0]),
                            constrained_layout=True)
    
    plt.subplots_adjust(hspace=0.25, wspace=0.2)
    
    # 1st loop - Itirating over the amino acid positions in the input list (pos_plot)
    # 2nd loop - Itirating over the datasets
    if n_pos <= 5:
        pos_iter_val = range(0,5)
    else:
        pos_iter_val = matrix_itr

    for pos_i, pos_itr in zip(pos_plot, pos_iter_val):
        if n_pos == 1:
            axi = axs
        else:
            axi = axs[pos_itr]

        # Frequncies template table
        pos_freq = pd.DataFrame(index =aa_order)

        for key in list(input_freqs.keys()):

            # Creating the needed format of the frequncies table for easy bar-plot
            temp_df = input_freqs[key]
            temp_df = temp_df[temp_df.index == pos_i].T
            temp_df.columns = [key]
            pos_freq = pd.concat([pos_freq, temp_df], axis=1)
        
        # The joined frequnceies at position i for all the dataset
        pos_freq = pos_freq.T.reset_index(names="dataset")

        # bar plot for position i
        pos_freq.plot(x="dataset", 
                      kind="bar", 
                      stacked=True,
                      ax=axi,
                      color = colors,
                      label="test")
        
        for ncol,xi in zip(input_metric.columns, range(0,len(pos_freq.columns))):
            metric_info = input_metric.loc[input_metric.index == pos_i, ncol].values[0]
            axi.text(x=xi-0.05, y=1+0.01, s= str(metric_info), fontsize=15)
            
        
        nclones_rmdupe = nclones_df.T.drop_duplicates().T
        nclones_posi = nclones_rmdupe[nclones_rmdupe.index == pos_i].values
        for val,i in zip(nclones_posi, range(1,len(nclones_posi)+1)):
            val_str = ", ".join([str(i) for i in val])
            type_str =  ", ".join([i.split(".")[0] for i in nclones_rmdupe.columns])
            axi.text(x=-0.05, y=1+0.06, s=f" Unique Clones of {type_str.upper()} : {val_str}", fontsize=18)
       
        # Configuring the bar axis
        axi.legend().set_visible(False)
        axi.set_ylabel("Amino Acid Fraction",  fontsize=30)
        axi.set_xlabel("")
        if change_xticks:
            axi.set_xticklabels(["SN.Germline", "SN.Sequenced", "SP.Germline", "SP.Sequenced"])
        axi.tick_params(axis='x', rotation=45, labelsize=30)
        axi.tick_params(axis='y', rotation=0, labelsize=30)
        axi.set_title("Amino Acid Position {}".format(pos_i), y=1.035 ,fontsize=30)

    # dropping empty axis
    for dax in noaxis:    
        axs[dax].set_axis_off()

    if matrix_shape[0] == 1:
        x_legend = (matrix_shape[1]-pos_itr)
        y_legend = 0
    else:
        x_legend = matrix_shape[1]-pos_itr[1]
        y_legend = matrix_shape[0]/2 - 0.25

    leg_fontsize = 12*3/2
    fig.legend(loc=7, bbox_to_anchor=(1.05, 0.5) ,labels=aa_order, fontsize=leg_fontsize, title="Amino \n Acid", title_fontsize=leg_fontsize, reverse=True)
    
    metric_dics = {"rich":"Richness", "div":"Diversity"}
    fig.suptitle(y=1.05, t="{} Stacked Bar Plot of Amino Acid Fraction Per Position ({})".format(metric_dics[name_metric],title), fontsize=24)

    if save_fig:
        plt.savefig("output\\{}_sbar_{}.png".format(name_metric, title), bbox_inches='tight')

    return plt.show()

###
class DiversityAnalysis():
    aa_order = ['R', 'E', 'K', 'Q', 'D', 'N', 'S"', 'G', 'W', 'C', 'H', 'Y', 'A', 'T', 'P', "S'", 'V', 'M', 'I', 'L', 'F']

    def __init__(self,
                 input_df : pd.DataFrame,
                 name_metric : str,
                 split_by : str = None,
                 pos_aa_range : range = None,
                 split_subjects : bool = False) -> [dict, pd.DataFrame, pd.DataFrame]:
        # frequnceis function
    
        """
        A function that accepts mut_df format dataframe and returns dictionary of frequnceies (key=condition)
        where:
        1. input_df - mut_df dataframe.
        2. split_by - column name that will be used to filter the dataframe by its unique values. 
        3. pos_aa_range - amino acids position range (make sure the pos_aa is int and not float/str)
        """

        self.subjects = input_df.subject_id.unique()

        if split_subjects:
            self.df_subjects = [input_df[input_df["subject_id"] == i] for i in self.subjects]
            self.dic_labels = self.subjects
        else:
            self.df_subjects = [input_df]
            self.dic_labels = "all_subjects"

        self.results = {}

        for sdf, lbl in zip(self.df_subjects, self.dic_labels):
        
            if name_metric not in ["div","rich"]:
                raise ValueError('{} not rich (richness) or div (diversity), enter correct value.'.format(name_metric))

            # Creating the template dataframe for the amino acid frequencies
            self.pos_aa_range = range(sdf.pos_aa.min(),sdf.pos_aa.max()+1)
            self.df = sdf[sdf["pos_aa"].isin(self.pos_aa_range)]
            self.freq_df_template = pd.DataFrame(data=0, index = self.pos_aa_range, columns=aa_order).astype(float)
            self.freq_dics = {}

            # Creating filtring condition for the output
            if split_by is None:
                self.cond_tempdf = ["from_aa.", "to_aas."]
            else:
                self.cond_tempdf = np.sort([".".join([by, col]) for col in ["from_aa","to_aas"] for by in self.df[split_by].unique()])

            # Creating diversity and richness dataframe to be filled later
            self.metric_df_template = pd.DataFrame(index = self.pos_aa_range, columns = self.cond_tempdf)
            self.metric_df = self.metric_df_template.copy()
            self.nclones_df = self.metric_df_template.copy()
            
            # Itirating over the conditions
            for c in self.cond_tempdf:
                self.c_split = list(filter(None,c.split(".")))
                self.temp_freq = self.freq_df_template.copy()
                
                # determining the input dataframe (according to the filtring) and output name.
                if len(self.c_split) == 1:
                    self.temp_df = self.df
                    self.name = c[:-1]
                    self.split_i = 0

                elif len(self.c_split) == 2:
                    self.temp_df = self.df[self.df[split_by]==self.c_split[0]]
                    self.name = c
                    self.split_i = 1

                # Getting the information of the frequncies, diversity and richness and assiging it to the output.
                for i in self.pos_aa_range:
                    self.nclones_df.loc[self.nclones_df.index == i, c] = len(self.temp_df.loc[self.temp_df["pos_aa"] == i, "clone_id"].unique())

                    if name_metric == "rich":
                        self.frequencies = self.temp_df.loc[self.temp_df["pos_aa"] == i, self.c_split[self.split_i]].value_counts().to_frame()
                        self.fractions = round(self.frequencies/self.frequencies.sum(),3)
                        self.temp_freq.loc[self.temp_freq.index == i, self.fractions.T.columns] = self.fractions.T.values
                        self.metric_df.loc[self.metric_df.index == i, c] = self.frequencies.shape[0]
                    
                    elif name_metric == "div":
                        # finidng diversity value
                        self.frequencies_og = self.temp_df.loc[self.temp_df["pos_aa"] == i, self.c_split[self.split_i]].value_counts().to_frame()
                        self.diversity_val = round(diversity(self.frequencies_og.iloc[:,0]/self.frequencies_og.iloc[:,0].sum()))
                        self.metric_df.loc[self.metric_df.index == i, c] = self.diversity_val
                        
                        # assiging the diversity values (rounded) to the datasets.
                        self.frequencies_div = self.temp_df.loc[self.temp_df["pos_aa"] == i, self.c_split[self.split_i]].value_counts().to_frame().iloc[:self.diversity_val,]
                        self.fractions_div = round(self.frequencies_div/self.frequencies_div.sum(),3)
                        self.temp_freq.loc[self.temp_freq.index == i, self.fractions_div.T.columns] = self.fractions_div.T.values
                
                self.freq_dics[self.name] = self.temp_freq
                
            self.results[lbl] = [self.freq_dics, self.metric_df, self.nclones_df]

    def get_data(self):
        return self.results
    
    @staticmethod
    def bar_freqs(input_freqs : dict,
              input_metric : pd.DataFrame,
              nclones_df: pd.DataFrame,
              name_metric : str,
              title : str,
              pos_plot : list = None,
              save_fig : bool = False) -> plt.Figure:
        
        """
        """
        if name_metric not in ["div","rich"]:
            raise ValueError('{} not rich (richness) or div (diversity), enter correct value.'.format(name_metric))

        # Colors for the barplot legend
        n_colors = len(aa_order)
        colors = cm.tab20(np.linspace(0, 1, n_colors))
        colors = list(colors)
        colors[0] = "blue"
        colors[-1] = "lavender"

        # Determining the needed aa for the plot & creating parameters for figure
        if pos_plot is None:
            pos_plot = input_div.index.values
        else:
            pos_plot = pos_plot

        n_pos = len(pos_plot)
        matrix_shape, matrix_itr = penta_matrix(n_pos)
        all_range = [(i,j) for i in range(0,matrix_shape[0]) for j in range(0,matrix_shape[1])]
        noaxis = [i for i in all_range if i not in matrix_itr]

        # Initiating the figure
        fig, axs = plt.subplots(matrix_shape[0],
                                matrix_shape[1], 
                                figsize=(6*matrix_shape[1],8*matrix_shape[0]))
        
        plt.subplots_adjust(hspace=0.25, wspace=0.2)

        # 1st loop - Itirating over the amino acid positions in the input list (pos_plot)
        # 2nd loop - Itirating over the datasets
        if n_pos <= 5:
            pos_iter_val = range(0,5)
        else:
            pos_iter_val = matrix_itr

        for pos_i, pos_itr in zip(pos_plot, pos_iter_val):
            if n_pos == 1:
                axi = axs
            else:
                axi = axs[pos_itr]

            # Frequncies template table
            pos_freq = pd.DataFrame(index =aa_order)

            for key in list(input_freqs.keys()):

                # Creating the needed format of the frequncies table for easy bar-plot
                temp_df = input_freqs[key]
                temp_df = temp_df[temp_df.index == pos_i].T
                temp_df.columns = [key]
                pos_freq = pd.concat([pos_freq, temp_df], axis=1)
            
            # The joined frequnceies at position i for all the dataset
            pos_freq = pos_freq.T.reset_index(names="dataset")

            # bar plot for position i
            pos_freq.plot(x="dataset", 
                        kind="bar", 
                        stacked=True,
                        ax=axi,
                        color = colors,
                        label="test")
            
            for ncol,xi in zip(input_metric.columns, range(0,len(pos_freq.columns))):
                metric_info = input_metric.loc[input_metric.index == pos_i, ncol].values[0]
                axi.text(x=xi-0.05, y=1+0.01, s= str(metric_info))
                
            
            nclones_rmdupe = nclones_df.T.drop_duplicates().T
            nclones_posi = nclones_rmdupe[nclones_rmdupe.index == pos_i].values
            for val,i in zip(nclones_posi, range(1,len(nclones_posi)+1)):
                val_str = ", ".join([str(i) for i in val])
                type_str =  ", ".join([i.split(".")[0] for i in nclones_rmdupe.columns])
                axi.text(x=0, y=1+0.06, s=f"-> Unique Clones of {type_str} : {val_str}", fontsize=12)
        
            # Configuring the bar axis
            axi.legend().set_visible(False)
            axi.set_ylabel("Amino Acid Fraction",  fontsize=13)
            axi.set_xlabel("")
            axi.tick_params(axis='x', rotation=45)
            axi.set_title("Amino Acid Position {}".format(pos_i),fontsize=17, y=1.035)

        # dropping empty axis
        for dax in noaxis:    
            axs[dax].set_axis_off()

        if matrix_shape[0] == 1:
            x_legend = (matrix_shape[1]-pos_itr)
            y_legend = 0
        else:
            x_legend = matrix_shape[1]-pos_itr[1]
            y_legend = matrix_shape[0]/2 - 0.25

        leg_fontsize = 12*matrix_shape[0]/2
        fig.legend(loc=7 ,labels=aa_order, fontsize=leg_fontsize, title="Amino \n Acid", title_fontsize=leg_fontsize, reverse=True)
        
        metric_dics = {"rich":"Richness", "div":"Diversity"}
        fig.suptitle(y=0.98-matrix_shape[0]*0.01-0.01, t="{} Stacked Bar Plot of Amino Acid Fraction Per Position ({})".format(metric_dics[name_metric],title), fontsize=24)

        if save_fig:
            plt.savefig("output\\{}_sbar_{}.png".format(name_metric, title), bbox_inches='tight')

        return plt.show()
    
    def plot_bars(self,
                  pos_aa : list,
                  save_figure : bool = False):
        
        for key in self.results.keys():
            title_str = f"subj_{str(key)}_head10"
            dict_frac, df_metric, nclones_df = self.results[key]
            bar_freqs(dict_frac, df_metric, nclones_df, "div", pos_plot=pos_aa, title=title_str, save_fig=save_figure)

    def plot_radar(self, pos_aa:int, save_figure:bool = False):
        subplots_n = len(self.results)
        subplots_labels = self.results.keys()

        n_rows = int(subplots_n/2  + (subplots_n%2)/2)
        axs_itr_og = [(nr, nc) for nr in range(0,n_rows) for nc in range(0,2)]
        axs_itr = axs_itr_og[0:subplots_n]

        fig, axs = plt.subplots(n_rows,2, subplot_kw={'projection': 'polar'}, figsize=(10*n_rows, 30))
        plt.subplots_adjust(hspace= 0.3, wspace= -0.5)
        theta = np.deg2rad(np.linspace(0,360,len(aa_order)+1))

        for axitr, key in zip (axs_itr, subplots_labels):
            temp_fractions, temp_metrics,  temp_nclones = self.results[key]

            fractions = pd.concat([temp_fractions[k][temp_fractions[k].index == pos_aa] for k in temp_fractions.keys()])
            fractions.index = temp_fractions.keys()
            fractions["R.D"] = fractions["R"]

            snval, spval= temp_nclones[temp_nclones.index == pos_aa].values[0][[0,2]]

            axs[axitr[0],axitr[1]].set_prop_cycle(color=['tab:blue', 'tab:green', '"tab:purple"', "tab:red"],
                                                  marker=['o', 'D', '^', "*"])

            for row, label in zip(fractions.iterrows(), fractions.index):
                vals = row[1].values
                axs[axitr[0],axitr[1]].plot(theta, 
                                            vals,
                                            #marker='o', 
                                            alpha=0.6, 
                                            markersize=15, 
                                            lw=2, 
                                            label=label)

                # Filling the shape
                axs[axitr[0],axitr[1]].fill(theta, vals, alpha=0.2, lw=0)

            # Setting the axis limit
            axs[axitr[0],axitr[1]].set_ylim(0,1)
            axs[axitr[0],axitr[1]].tick_params(axis="y", labelsize=15, color="grey")
            axs[axitr[0],axitr[1]].set_xticks(theta)
            axs[axitr[0],axitr[1]].set_xticklabels(aa_order+[""], size=20)
            axs[axitr[0],axitr[1]].set_title(f"Subject {key} \nUnique Clones:  Non-Spike = {str(snval)} , Spike = {str(spval)}\n", size= 20)

            if (len(axs_itr_og)-1) == subplots_n:
                axs[axs_itr_og[-1][0],axs_itr_og[-1][1]].set_axis_off()
        
        fig.suptitle(f'Radar Plot of Amino Acid Fractions at Position {str(pos_aa)}', fontsize=20, y=0.96)

        handles, labels = axs[axs_itr[0][0],axs_itr[0][1]].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center',bbox_to_anchor=(0.5, 0.94), fontsize=13, ncols=len(handles))
        
        if save_figure:
            plt.savefig(f"output\\radar_{str(pos_aa)}.png", bbox_inches='tight')

        plt.show()

    def do_wilcoxon(self, pos_aa:list = None, stats:str="to_aas"):
    
        if pos_aa is None:
            pos_aa = self.input_df.pos_aa.astype("int64").unique()

        tdiv_list = []
        raw_list = []
        div_datasets = self.results

        for key, item in div_datasets.items():
            temp_df = item[1].copy().astype("int64")
            temp_df = temp_df.loc[temp_df.index.isin(pos_aa),:]
            raw_df = temp_df.copy()
            raw_df["subject"] = key
            raw_list.append(raw_df)
            n_df = temp_df.shape[0] # degree of freedom n - groups
            statistic, pval = wilcoxon(x=temp_df[f"sn.{stats}"], y=temp_df[f"sp.{stats}"], alternative="greater") # non-parametric paired test
            tdiv_list.append([key, statistic, pval, n_df, n_df-2])
            wilcoxon_df = pd.DataFrame(data=tdiv_list, columns=[["subject","wilcoxon_w","p_value","n","dof"]])
    
        return wilcoxon_df, pd.concat(raw_list)


## amino acid to protein sequence translation function
def nt2aa(nt_seq : str) -> str:
    # Protein sequence output
    protein = []
    
    # Codon translation dictionary
    codon_dict =  {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
                   'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
                   'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
                   'AGC':'S"', 'AGT':'S"', 'AGA':'R', 'AGG':'R',                  
                   'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
                   'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
                   'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
                   'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
                   'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
                   'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
                   'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
                   'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
                   'TCA':"S'", 'TCC':"S'", 'TCG':"S'", 'TCT':"S'", 
                   'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
                   'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
                   'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    # Cheeking the raminder of length sequnce division by 3 for translation
    rem3 = len(nt_seq)%3
    if rem3 != 0:
        nt_seq = nt_seq[:-rem3]
    
    # Verifing that the sequence length is divided by 3, if not getting the highest number that does.
    for i in range(1,len(nt_seq),3):
        i_start = i-1
        i_end = i+2
        codon = nt_seq[i_start:i_end]

        try:
            protein.append(codon_dict[codon])

        except:
            protein.append("*")
    
    return protein


##
# filt_df[filt_df.ab_target == "sp"].groupby("pos_aa").agg({"from_aa":div_tuple, "to_aas":div_tuple}).reset_index(drop=False)
def get_setinter(df : pd.DataFrame,
                 split_by : str,
                 pos_aa : list = range(1,105),
                 upto : int = 104) -> pd.DataFrame:

    # Extracting nt values up to 104 (defualt value)
    extr_to104 = lambda x : x[0][:upto*3]
    # Extracint item i from string
    get_i = lambda str, i : str[i]

    # Creating unique germine sequence dataframe that will be translated later
    uniques_vals = df[split_by].unique()
    dfs_concat = []

    # Get filtred input df for want wanted aa_pos
    df = df[df["pos_aa"].isin(pos_aa)]

    # Iterating over the unique values of the split_by column
    for i in uniques_vals:
        # Creating germline sequence (unique) for each dataset and: 
        # 1) getting nucleotides up to pos 104 (defualt)
        # 2) Translating to amino acids 
        germline_raw = df.groupby(["clone_id",split_by]).agg({"germline":"unique"}).reset_index().drop_duplicates(subset=["germline"]).loc[:,[split_by,"germline"]].reset_index()
        germline_raw = germline_raw[germline_raw[split_by]==i]
        germline_raw[f"germline"] = germline_raw["germline"].apply(extr_to104)
        germline_raw["protein"] = germline_raw["germline"].apply(nt2aa)

        # Creating empty datadrame for the data
        temp_df = pd.DataFrame(index=pos_aa)

        # Getting the unique amino acid values for each position before (from_aa) and after (to_aas) mutations.
        # Merging the information to the temp df.
        vals_fromto = df.groupby("pos_aa").agg({"from_aa":"unique","to_aas":"unique"}).reset_index().set_index("pos_aa")
        temp_df = temp_df.merge(vals_fromto, left_index=True, right_index=True, how="left")

        # Getting unique values of the germline for each aa_position
        temp_df["germline"] = np.nan
        germ_list = []
        for iter in pos_aa:
            germ_list.append(germline_raw["protein"].apply(get_i, i=iter-1).unique())
    
        # appending the temp_df to the list for future concantination
        temp_df["germline"] = germ_list
        temp_df = temp_df[["germline","from_aa","to_aas"]].add_suffix(f"_{i}")
        dfs_concat.append(temp_df)
  
    # Merging, replacing nan with empty array and removing "_" (stop) and "*" (missing) from possible aa values.
    return_values = pd.concat(dfs_concat, axis=1)
    return_values = return_values.map(lambda x : np.array([]) if x is np.nan else x[(x != "_") & (x != "*")]).reset_index(drop=False, names="pos_aa")

    range_cdr = np.concatenate((np.arange(27,39,1), np.arange(56,66,1)))
    return_values.insert(loc=1, column="region", value=pd.Series(["cdr" if i in range_cdr else "fw" for i in return_values.pos_aa.values]))


    return return_values

##
def div_tuple(ser_input):
    aa_count = ser_input.value_counts()
    aa_fraction = aa_count / aa_count.sum()
    aa_div = int(round(diversity(aa_fraction),0))

    div_count = aa_count[:aa_div]
    div_fractions = div_count / div_count.sum()
    
    return [(aa, "{:.2f}".format(frac)) for aa,frac in zip(div_fractions.index, div_fractions.values)]