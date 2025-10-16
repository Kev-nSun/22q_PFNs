Code used to derive PFNs from 2 timepoints of ABCD. 
First, merge_dtseries_by_session.py is used to manually concatenate timeseries using wb_command -cifti-merge. This script sorts by subject and session, and then concatenates all scans within regardless of how many scans are in the session directory. The script also creates a txt file called 19_Scan_List_Concat.txt that is fed into config.toml as the variable file_scans.

To run merge_dtseries_by_session.py:

  python merge_dtseries_by_session.py \
      --root /cbica/projects/PFN_ABCD/abcd-hcp-pipeline_0.1.4_timeseries \
      --out  /cbica/projects/PFN_ABCD/long_PFN_scripts/pnet_inputs/merged_dtseries

Then, to run pnet, make sure config.toml is in the directory:

  conda activate fmripnet

  python /cbica/projects/PFN_ABCD/pNet/fmripnet.py -c config.toml --hpc

PFN_Similarity_ARI.R calculates the adjusted Rand index (ARI) of PFNs between and within subjects, across timepoints and pipelines, and visualizes these results.
