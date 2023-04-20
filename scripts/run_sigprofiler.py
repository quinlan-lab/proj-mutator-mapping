from SigProfilerExtractor import sigpro as sig
import argparse
from SigProfilerMatrixGenerator import install as genInstall

p = argparse.ArgumentParser()
p.add_argument(
    "--spectra",
    required=True,
    help=
    """Path to file containing mutation spectra counts in samples in SigProfiler format""",
)
p.add_argument(
    "--outdir",
    required=True,
    help="""Name of output directory to store results""",
)
args = p.parse_args()

if __name__ == "__main__":

    #genInstall.install('mm10')

    sig.sigProfilerExtractor(
        'matrix',
        args.outdir,
        args.spectra,
        maximum_signatures=10,
        nmf_replicates=100,
        opportunity_genome="mm10",
    )