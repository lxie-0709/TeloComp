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

def run_step1_2(ont_bam, hifi_bam, outdir):
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
        "-o", outdir
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
    parser = argparse.ArgumentParser(
        description="Wrapper: run step1_2 + step1_3 in sequence (imported modules)"
    )
    parser.add_argument("--ont_bam",  required=True, help="ONT BAM")
    parser.add_argument("--hifi_bam", required=True, help="HiFi BAM")
    parser.add_argument("-o","--out_dir", required=True,
                        help="Directory where both step1_2 and step1_3 should write their outputs")
    parser.add_argument("-c","--coverage", type=float, default=100,
                        help="Coverage parameter for step1_3 (passed to -c)")
    parser.add_argument("-p","--parallels",  type=int,   default=10,
                        help="Parallels for step1_3 (passed to -p)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s"
    )

    # Run step1_2
    logging.info(f"=== Running step1_2 → output in {args.out_dir}")
    ont_dir, hifi_dir = run_step1_2(
        args.ont_bam, args.hifi_bam, args.out_dir
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
