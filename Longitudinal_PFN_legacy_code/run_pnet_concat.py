# conda activate fmripnet
# python run_pnet_concat.py config.toml
import sys, tomli, pnet

cfg = tomli.load(open(sys.argv[1], "rb"))
ns   = cfg["necessary_settings"]
pfns  = cfg["pFN_settings"]
gfns  = cfg["gFN_settings"]
envs = cfg["hpc_settings"]["pnet_env"]
sub  = cfg["hpc_settings"]["submit"]
cres = cfg["hpc_settings"]["computation_resource"]

pnet.workflow_cluster(
    dir_pnet_result     = ns["dir_pnet_result"],
    dataType            = ns["dataType"],
    dataFormat          = ns["dataFormat"],
    file_Brain_Template = ns["file_Brain_Template"],
    file_scan           = ns["file_scans"],
    file_subject_ID     = None,
    file_subject_folder = None,
    file_gFN            = None if pfns["file_gFN"] == "None" else pfns["file_gFN"],
    K                   = ns["K"],
    Combine_Scan        = True, # <- force concatenation of scans
    method              = ns["method"],
    init                = "random",
    sampleSize          = gfns["sampleSize"],
    nBS                 = gfns["nBS"],
    nTPoints            = gfns["nTPoints"],
    Computation_Mode    = "CPU_Torch",
    # HPC bits:
    dir_env             = envs["dir_env"],
    dir_python          = envs["dir_python"],
    dir_pnet            = envs["dir_pnet"],
    submit_command      = sub["submit_command"],
    thread_command      = sub["thread_command"],
    memory_command      = sub["memory_command"],
    log_command         = sub["log_command"],
    computation_resource= cres,
)