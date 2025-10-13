Code used to derive PFNs from 2 timepoints of ABCD. 
First, stage_symlinks.py was used to creat symlinks of the path to scans that created parents directories of subject+session and also creates the txt file (file_scans) used in config.toml.
pNet is run using run_pnet_concat.py, a modified script of fmripnet.py: https://github.com/MLDataAnalytics/pNet/blob/main/fmripnet.py, in which Combine_Scan=true, file_subject_ID = None, and file_subject_folder = None.

To run run_pnet_concat.py, you need to activate the fmripnet conda environment and have the config.toml file in the same directory:
conda activate fmripnet
python run_pnet_concat.py config.toml


PFN_Similarity_ARI.R calculates the adjusted Rand index (ARI) of PFNs between and within subjects, across timepoints and pipelines, and visualizes these results.
