GLOBAL_PREFIX = "honeycomb"
GLOBAL_OUTDIR = "results"

def set_global_prefix(prefix: str):
    """
    Set the output prefix for files
    
    Args:
        prefix (str): The prefix to be used
    """
    global GLOBAL_PREFIX
    GLOBAL_PREFIX = prefix


def set_global_outdir(outdir: str):
    """
    Set the output directory for files
    
    Args:
        prefix (str): The outdir to be used
    """
    global GLOBAL_OUTDIR
    GLOBAL_OUTDIR = outdir
