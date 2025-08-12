#!/usr/bin/env python3
import os
import sys
import argparse
import logging


# Ensure this directory is on Python’s import path
HERE = os.path.dirname(__file__)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# Now we can import the step2 and step3 scripts as modules
import telocomp_Det_Ext as step1_2
import telocomp_trim as step1_3

def run_step1_2(ont_bam, hifi_bam, outdir, min_ratio):
    """
    Runs the second-stage script (motif discovery + overhang extraction)
    in‐process by calling its main().
    Produces subdirs outdir/ONT and outdir/HiFi.
    """
    os.makedirs(outdir, exist_ok=True)
    # build argv for step1_2
    old_argv = sys.argv
    sys.argv = [
        old_argv[0],
        "--ont", ont_bam,
        "--hifi", hifi_bam,
        "-o", outdir,
        "--min_ratio", str(min_ratio)
    ]
    try:
        step1_2.main()
    finally:
        sys.argv = old_argv
    return os.path.join(outdir, "ONT"), os.path.join(outdir, "HiFi")

def run_step1_3(ont_dir, hifi_dir, outdir, coverage, parallels):
    """
    Runs the third-stage script (trim & merge) in‐process by calling its main().
    """
    # build argv for step1_3
    old_argv = sys.argv
    sys.argv = [
        old_argv[0],
        "--ont_dir",  ont_dir,
        "--hifi_dir", hifi_dir,
        "-c", str(coverage),
        "-o", outdir,
        "-t", str(parallels)
    ]
    try:
        step1_3.main()
    finally:
        sys.argv = old_argv

def main():        
    # 自定义 Formatter，不显示 metavar
    class NoMetavarFormatter(argparse.HelpFormatter):
        def _format_action_invocation(self, action):
            if not action.option_strings:
                return super()._format_action_invocation(action)

            # 对于带参数的选项，只显示选项名，不显示 metavar
            return ', '.join(action.option_strings)
    
    parser = argparse.ArgumentParser(
        description="Your description here",
        formatter_class=NoMetavarFormatter
    )
    
    parser.add_argument("--ont_bam",  required=True, help="ONT BAM")
    parser.add_argument("--hifi_bam", required=True, help="HiFi BAM")
    parser.add_argument("-o","--out_dir", required=True, 
                        help="Directory where both step1_2 and step1_3 should write their outputs")
    parser.add_argument("-c","--coverage", type=float, default=100, 
                        help="Coverage parameter for step1_3 (passed to -c)")
    parser.add_argument("-p","--parallels",  type=int,   default=10,
                        help="Parallels for step1_3 (passed to -p)")
    parser.add_argument('--min_ratio', type=float, default=0.2, 
                        help="The proportion of the original genome sequence to the length of the reads, default=0.2")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s"
    )

    # Run step1_2
    logging.info(f"=== Running step1_2 → output in {args.out_dir}")
    ont_dir, hifi_dir = run_step1_2(
        args.ont_bam, args.hifi_bam, args.out_dir, args.min_ratio
    )

    # 同样解析成绝对路径，防止相对路径混淆
    ont_dir  = os.path.abspath(ont_dir)
    hifi_dir = os.path.abspath(hifi_dir)

    # Run step1_3
    logging.info(f"=== Running step1_3 on {ont_dir} & {hifi_dir}")
    run_step1_3(ont_dir, hifi_dir, args.out_dir, args.coverage, args.parallels)

    logging.info("All done.")

if __name__ == "__main__":
    main()
