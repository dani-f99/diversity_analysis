-----------------------------
------Daniel-Fridman---------
--- System-Immunology-Lab --- 
------Haifa-University ------
-----------2025--------------
-----------------------------
-----------------------------
--- General-Information -----
This program's aim is to visualize the diversity patterns across the BCR heavy-chain variable region 
of BCR datasets originating from ImmuneDB tables.
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- Required-Python-Modules -
1.  configparser
2.  matplotlib
3.  scipy 
4.  pandas 
5.  numpy 
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- Guide ---
1. Configure the config.ini file as needed (see config section below).
2. provide the proper input tables (clones, clone_stats and sample_metadata) download from the ImmuneDB MySQL
   server or use 'lpa_preprocessing' / "substitution_survival_analysis" programs.
3. Run the notebook
4. Output figures will be saved into the output folder.
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- config.ini --------------
1. [data_path]
  a. path_clones -> path of the clones ImmuneDB tables, downloaded from the MySQL server.
  b. path_clone_stats -> path of the clone_stats ImmuneDB tables, downloaded from the MySQL server.
  c. path_sample_metadata -> path of the sample_metadata ImmuneDB tables, downloaded from the MySQL server.
-----------------------------
2. [figures]
  a. save_fig = boolean value, if True the figures will be saved into (folder will be creared via the code in cell no.2)
-----------------------------
-----------------------------
-----------------------------
-----------------------------
--- Subfolders --------------
1. output – output figures folder.
2. input - raw input tables derived from the MySQL server.
-----------------------------
-----------------------------
-----------------------------
-----------------------------
---Functions & Classes-------
1. helpers.analysis_information -> A function that accepts mut_df format dataframe and returns dictionary of frequnceies (key=condition)
2. helpers.bar_freqs -> Plotting bars which represents the frequencies of amino acids according to the analysis_information output.
3. helpers.DiversityAnalysis -> Class inititation, will modify and preprocess the mut_df information to analyze the clones diversity.
4. helpers.DiversityAnalysis.plot_radar -> plotting radar plot showing the amino acid fraction in specific position.
-----------------------------
-----------------------------

