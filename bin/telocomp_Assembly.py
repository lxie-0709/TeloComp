#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
telo_comp_assemble.py

Depending on user flags, runs:
  - Flye-based assembly (with shortest-read fallback) when --flye is specified,
  - k-mer overlap assembly when --assemble/-a is specified,
on every FASTA in trim_L/ and trim_R/ in parallel.  After that, runs
the original polishing step unchanged.

Usage:
    python telo_comp_assemble.py \
      --dir_IN_L trim_L \
      --dir_IN_R trim_R \
      --flye \
      --assemble \
      -t 20 \
      -L long_reads.fq \
      -W read1.fq \
      -w read2.fq \
      -N /path/to/NextPolish
"""
import os
import re
import shutil
import argparse
import subprocess
import logging
import time
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
import itertools
from collections import deque


# module-level logger
log = logging.getLogger()


# --- Block #1: Flye + shortest reads fallback ---------------------------------
# Extract the shortest reads as contig
def extract_shortest_reads(fasta_path, output_path, seq_type, base_name, direction):
    log.info(f"Fallback: extracting shortest reads from {fasta_path}")
    if seq_type not in ["ONT", "HiFi", "merged"]:
        raise ValueError("seq_type must be one of 'ONT','HiFi','merged'")
    ont_out = os.path.join(output_path, f"{base_name}_shortest_ONT_{direction}.fasta")
    hifi_out = os.path.join(output_path, f"{base_name}_shortest_HiFi_{direction}.fasta")
    ont_h = open(ont_out, 'w') if seq_type=="merged" else None
    hifi_h = open(hifi_out,'w') if seq_type=="merged" else None
    out_h = None
    if seq_type!="merged":
        out = os.path.join(output_path, f"{base_name}_shortest_{seq_type}_{direction}.fasta")
        out_h = open(out,'w')
    seen_ont = seen_hifi = False
    for r in SeqIO.parse(fasta_path,'fasta'):
        if "min" in r.id:
            if seq_type=="merged":
                if "ONT" in r.id:
                    seen_ont=True
                    ont_h.write(f">{r.id}\n{r.seq}\n")
                elif "HiFi" in r.id:
                    seen_hifi=True
                    hifi_h.write(f">{r.id}\n{r.seq}\n")
            else:
                out_h.write(f">{r.id}\n{r.seq}\n")
    if seq_type=="merged":
        if ont_h: ont_h.close()
        if hifi_h: hifi_h.close()
        if not seen_ont and os.path.exists(ont_out): 
            os.remove(ont_out)
            log.debug(f"Removed empty ONT fallback file {ont_out}")
        if not seen_hifi and os.path.exists(hifi_out): 
            os.remove(hifi_out)
            log.debug(f"Removed empty HiFi fallback file {hifi_out}")        
    else:
        if out_h: out_h.close()

def safe_rmtree(path, max_retries=3, delay=1):
    """Safely remove directory tree with retries"""
    for i in range(max_retries):
        try:
            if os.path.exists(path):
                shutil.rmtree(path)
                log.debug(f"Successfully removed directory {path}")
                return True
        except (OSError, PermissionError) as e:
            if i < max_retries - 1:
                log.warning(f"Attempt {i+1} to remove {path} failed: {e}. Retrying...")
                time.sleep(delay)
            else:
                log.error(f"Failed to remove {path} after {max_retries} attempts: {e}")
                return False
    return True

def run_flye_one(fasta_path, threads):
    log.info(f"Running Flye on {fasta_path} (threads={threads})")
    fname = os.path.basename(fasta_path)
    base = fname.split('_')[0]
    seq_type = fname.split('_')[2]
    direction = 'L' if '_L' in fname else 'R'
    asm_dir = f"asm_{direction}"
    os.makedirs(asm_dir, exist_ok=True)
    out_pref = os.path.join(asm_dir, f"{base}_asm_{direction}")
    
    # First, remove any existing assembly or shortest files for this base+direction
    existing_files = glob.glob(os.path.join(asm_dir, f"{base}_*_{direction}.fasta"))
    for f in existing_files:
        try:
            os.remove(f)
            log.debug(f"Removed existing file: {f}")
        except OSError as e:
            log.warning(f"Could not remove existing file {f}: {e}")
    
    cmd = f"flye --meta --no-alt-contigs --nano-raw {fasta_path} --out-dir {out_pref} --threads {threads}"
    try:
        subprocess.run(cmd, shell=True, check=True)
        asm_f = os.path.join(out_pref, "assembly.fasta")
        newf = os.path.join(asm_dir, f"{base}_asm_{seq_type}_{direction}.fasta")
        if os.path.exists(asm_f):
            os.rename(asm_f, newf)
            log.info(f"Flye succeeded: wrote {newf}")
    except subprocess.CalledProcessError as e:
        log.warning(f"Flye failed on {fasta_path}, invoking shortest-read fallback. Error: {e}")
        # Remove any partial assembly files that might have been created
        partial_files = glob.glob(os.path.join(asm_dir, f"{base}_asm_*_{direction}.fasta"))
        for f in partial_files:
            try:
                os.remove(f)
                log.debug(f"Removed partial assembly file: {f}")
            except OSError:
                pass
        extract_shortest_reads(fasta_path, asm_dir, seq_type, base, direction)
    finally:
        # Use safe directory removal with retries
        safe_rmtree(out_pref)

# --- Block #2: k-mer overlap assembler ----------------------------------------

class SequenceAssembler:
    def __init__(self, min_overlap=50, error_rate=0.15, kmer_size=15):
        self.min_overlap, self.error_rate, self.kmer_size = min_overlap, error_rate, kmer_size
        log.debug(
            f"Initialized SequenceAssembler(min_ov={min_overlap}, "
            f"err={error_rate}, kmer={kmer_size})"
        )

        
    def read_fasta(self, fn):
        seqs, cur = [], []
        with open(fn) as f:
            for l in f:
                if l[0]==">":
                    if cur: seqs.append("".join(cur)); cur=[]
                else:
                    cur.append(l.strip())
            if cur: seqs.append("".join(cur))
        return seqs
    def _kmer_idx(self, s):
        return {s[i:i+self.kmer_size]:i for i in range(len(s)-self.kmer_size+1)}
    def _extend(self, s1,s2,i,j):
        m=t=0
        while i<len(s1) and j<len(s2) and t< len(s1)*2:
            if s1[i]==s2[j]: m+=1
            t+=1; i+=1; j+=1
        return (m/t) if t>=self.min_overlap else 0
    def find_best_overlap(self,a,b):
        best=(-1,None,None)
        idx=self._kmer_idx(b)
        for i in range(len(a)-self.kmer_size,-1,-1):
            seed=a[i:i+self.kmer_size]
            if seed in idx:
                sc=self._extend(a,b,i,idx[seed])
                if sc>best[0]: best=(sc,0,i)
        idx=self._kmer_idx(a)
        for i in range(len(b)-self.kmer_size,-1,-1):
            seed=b[i:i+self.kmer_size]
            if seed in idx:
                sc=self._extend(b,a,i,idx[seed])
                if sc>best[0]: best=(sc,1,i)
        if best[0]>= self.min_overlap*(1-self.error_rate):
            return best
        return None
    def assemble(self, seqs):
        dq=deque(sorted(seqs,key=len,reverse=True))
        log.info("Starting k-mer overlap assembly")
        while len(dq)>1:
            bs=-1; bp=None; merged=""
            for i,j in itertools.permutations(range(len(dq)),2):
                r=self.find_best_overlap(dq[i],dq[j])
                if r and r[0]>bs:
                    sc,dr,pos = r
                    merged = dq[i]+dq[j][pos+self.kmer_size:] if dr==0 else dq[j]+dq[i][pos+self.kmer_size:]
                    bs, bp = sc, (i,j)
            if not bp: break
            i,j=sorted(bp,reverse=True)
            dq.remove(dq[j]); dq.remove(dq[i])
            dq.appendleft(merged)
        log.info("Finished k-mer assembly")
        return dq[0] if dq else ""

def run_kmer_one(fasta_path, assembler):
    log.info(f"Running k-mer assembler on {fasta_path}")
    fname=os.path.basename(fasta_path)
    base=fname.split('_')[0]
    direction='L' if '_L' in fname else 'R'
    seqs=assembler.read_fasta(fasta_path)
    outdir=f"asm_{direction}"
    os.makedirs(outdir, exist_ok=True)
    
    # Remove any existing k-mer assembly for this base+direction
    existing_files = glob.glob(os.path.join(outdir, f"{base}_kmer_{direction}.fasta"))
    for f in existing_files:
        try:
            os.remove(f)
            log.debug(f"Removed existing k-mer file: {f}")
        except OSError as e:
            log.warning(f"Could not remove existing k-mer file {f}: {e}")
    
    if len(seqs)==1:
        dest = os.path.join(outdir, fname)
        shutil.copy(fasta_path, dest)
        log.info(f"Only one read; copied {dest}")
    else:
        contig=assembler.assemble(seqs)
        out_path = os.path.join(outdir, f"{base}_kmer_{direction}.fasta")
        with open(out_path, 'w') as w:
            w.write(f">{base}_kmer_{direction}\n{contig}\n")
        log.info(f"Wrote k-mer contig to {out_path}")

# --- Block #3: polishing (unchanged) ------------------------------------------

def clean_ansi_codes(text):
    """Remove ANSI escape codes from text"""
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    return ansi_escape.sub('', text)

def asm_contig_racon(input_file, lgsreads, wgs1, wgs2, threads, NextPolish):
    log.info(f"Starting polishing of {input_file}")
    tmp="tmp_dir"; os.makedirs(tmp,exist_ok=True)
    while True:
        cmds=[
            f"minimap2 -ax map-pb -t {threads} {input_file} {lgsreads} | samtools sort - -m 2g -o {tmp}/genome.lgs.bam",
            f"samtools index {tmp}/genome.lgs.bam",
            f"ls {os.path.abspath(tmp)}/genome.lgs.bam > {tmp}/pb.map.bam.fofn",
            f"python {NextPolish}/lib/nextpolish2.py -g {input_file} -l {tmp}/pb.map.bam.fofn -r clr -p {threads} -o {tmp}/merged_genome.lgspolish.fasta"
        ]
        try:
            for c in cmds:
                # Capture stderr to handle errors properly
                result = subprocess.run(c, shell=True, check=True, capture_output=True, text=True)
            break
        except subprocess.CalledProcessError as e:
            # Properly handle stderr that might be None
            err = e.stderr if e.stderr else str(e)
            if not err and e.stdout:
                err = e.stdout
            
            # Clean ANSI codes from error message
            clean_err = clean_ansi_codes(err)
            log.error(f"Command failed with error: {clean_err}")
            
            if "Failed to correct sequence" in clean_err:
                m=re.search(r'Failed to correct sequence:\s+(\S+)', clean_err)
                if not m: 
                    log.error("Polishing gave unknown failure")
                    return
                # Clean the sequence ID from ANSI codes
                fid = clean_ansi_codes(m.group(1))
                fid = re.sub(r'[^\w_\-]','', fid)
                log.warning(f"Polish failure on {fid}, removing and retrying")
                newf = remove_sequence_from_fasta(input_file, fid)
                if newf == input_file: 
                    log.error("Failed to remove sequence, aborting polishing")
                    return
                input_file = newf
                continue
            else:
                log.exception("Unrecoverable polish error")
                return
    outd="tmp_to_NP"; os.makedirs(outd,exist_ok=True)
    genome=f"{tmp}/merged_genome.lgspolish.fasta"
    if os.path.exists(genome):
        for rec in SeqIO.parse(genome,"fasta"):
            parts=rec.id.split('_')
            if len(parts) >= 4:
                chrn, et, d = parts[0], parts[1], parts[3]
                newf = f"{outd}/{chrn}_{et}_np_{d}.fasta"
                with open(newf,'w') as h:
                    SeqIO.write(rec,h,"fasta")
                log.debug(f"Wrote polish chunk {newf}")
            else:
                log.warning(f"Unexpected record ID format: {rec.id}")
    else:
        log.error(f"Polished genome file not found: {genome}")
        return
        
    safe_rmtree(tmp)
    asm_contig_polish(outd,lgsreads,wgs1,wgs2,threads,NextPolish)

def remove_sequence_from_fasta(fasta_file, seq_id):
    log.info(f"Removing failed seq {seq_id} from {fasta_file}")
    out_new="new_merged_chr_shortest.fasta"; od="tmp_to_NP"; os.makedirs(od,exist_ok=True)
    clean=re.sub(r'[^\w_\-]','',seq_id)
    removed=False
    with open(fasta_file) as inp, open(out_new,'w') as out:
        for rec in SeqIO.parse(inp,"fasta"):
            rid=re.sub(r'[^\w_\-]','',rec.id)
            if clean in rid:
                with open(f"{od}/{rid}.fasta",'w') as fh: 
                    SeqIO.write(rec,fh,"fasta")
                removed=True
                log.debug(f"Extracted failed record {rid}")
            else:
                SeqIO.write(rec,out,"fasta")
    return out_new if removed else fasta_file

def asm_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish):
    import glob
    log.info(f"Starting final merge/polish in {output_dir}")
    for f in glob.glob(f"{output_dir}/*.fasta"):
        name=os.path.splitext(os.path.basename(f))[0]
        FoFn=f"Dir_fofn/{name}"; os.makedirs(FoFn,exist_ok=True)
        shutil.copy(f, f"{FoFn}/{os.path.basename(f)}")
        sgs=f"{FoFn}/sgs.fofn"; lgs_file=f"{FoFn}/lgs.fofn"
        with open(sgs,'w') as h: h.write(wgs1+"\n"+wgs2+"\n")
        with open(lgs_file,'w') as h: h.write(lgsreads+"\n")
        t = threads
        work_dir = os.path.join(FoFn, "work")
        os.makedirs(work_dir, exist_ok=True)

        current_dir = os.getcwd()
        final_dir = os.path.join(current_dir, "files_NP")
        os.makedirs(final_dir, exist_ok=True)
        
        cfg=f"{FoFn}/run.cfg"
        with open(cfg,'w') as h:
            h.write(f"""[General]
   job_type=local
        job_prefix = nextPolish
        task = best
        rewrite = no
        rerun = 3
        parallel_jobs = 1
   multithread_jobs={threads}
   genome={os.path.basename(f)}
        genome_size = auto
        workdir = {work_dir}
        polish_options = -p {t}

        [sgs_option]
        sgs_fofn = sgs.fofn
        sgs_options = -max_depth 100 -bwa

        [lgs_option]
        lgs_fofn = lgs.fofn
        lgs_options = -min_read_len 5k -max_depth 100
        lgs_minimap2_options = -x map-pb
""")
        cmd=f"{NextPolish}/nextPolish {cfg}"
        try:
            subprocess.run(cmd,shell=True,check=True,capture_output=True,text=True)
            log.info(f"NextPolish succeeded in {FoFn}")
            for npf in os.listdir(f"{FoFn}/"):
                if npf.endswith(".fasta"):
                    shutil.move(f"{FoFn}/{npf}","files_NP/"+name+".fasta")
        except subprocess.CalledProcessError as e:
            log.warning(f"NextPolish failed for {name}, copying raw contig. Error: {e.stderr if e.stderr else str(e)}")
            shutil.copy(f,"files_NP/"+name+".fasta")
        safe_rmtree(FoFn)

def process_asm_merge_fasta_files(dirs, output_file, lgsreads, wgs1, wgs2, threads, NextPolish):
    uniq=set()
    with open(output_file,'w') as out:
        for d in dirs:
            if not os.path.exists(d):
                log.warning(f"Directory {d} does not exist, skipping")
                continue
            for f in os.listdir(d):
                if not f.endswith('.fasta'): continue
                path=os.path.join(d,f)
                # Only process files that don't have "shortest" in the name if there's a corresponding assembly file
                if 'shortest' in f:
                    # Check if there's a corresponding assembly file
                    base_name = f.split('_')[0]
                    direction = 'L' if '_L' in f else 'R'
                    asm_pattern = f"{base_name}_asm_*_{direction}.fasta"
                    asm_files = glob.glob(os.path.join(d, asm_pattern))
                    # If there's an assembly file, skip the shortest file
                    if asm_files:
                        log.debug(f"Skipping shortest file {f} because assembly file exists")
                        continue
                
                if 'asm_' in f:
                    for rec in SeqIO.parse(path,'fasta'):
                        id0=os.path.splitext(f)[0]
                        if id0 not in uniq:
                            uniq.add(id0)
                            out.write(f">{id0}\n{rec.seq}\n")
                else:
                    for rec in SeqIO.parse(path,'fasta'):
                        s=''.join(c for c in str(rec.seq) if c.islower())
                        id0=os.path.splitext(f)[0]
                        if s and id0 not in uniq:
                            uniq.add(id0)
                            out.write(f">{id0}\n{s}\n")
    asm_contig_racon(output_file,lgsreads,wgs1,wgs2,threads,NextPolish)

# --- Dispatcher -------------------------------------------------------------

def dispatch(dir_L, dir_R, threads, do_flye, do_asm, min_ov, err, ksize):
    log.info(f"Dispatching, L={dir_L}, R={dir_R}, flye={do_flye}, asm={do_asm}")
    
    # Check if input directories exist
    if not os.path.exists(dir_L):
        log.error(f"Input directory {dir_L} does not exist")
        return
    if not os.path.exists(dir_R):
        log.error(f"Input directory {dir_R} does not exist")
        return
        
    fasta_files = [
        os.path.join(dir_L, f) for f in os.listdir(dir_L) if f.endswith('.fasta')
    ] + [
        os.path.join(dir_R, f) for f in os.listdir(dir_R) if f.endswith('.fasta')
    ]
    
    if not fasta_files:
        log.warning("No FASTA files found in input directories")
        return
        
    with ThreadPoolExecutor(max_workers=min(threads, len(fasta_files))) as exe:
        futures=[]
        if do_flye:
            for f in fasta_files:
                futures.append(exe.submit(run_flye_one, f, threads))
        if do_asm:
            asm=SequenceAssembler(min_ov, err, ksize)
            for f in fasta_files:
                futures.append(exe.submit(run_kmer_one, f, asm))
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                log.error(f"Assembly task failed: {e}")
    log.info("Assembly dispatch complete")

# --- Main -------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog='FlyePipeline',
        usage='%(prog)s [options]',
        description='Run Flye/assembly pipeline',
        epilog='Text at the bottom of help'
    )
    
    # Directory arguments for left and right reads
    parser.add_argument('--dir_IN_L', required=True, metavar='', help='Directory containing left-aligned reads (FASTA format)')
    parser.add_argument('--dir_IN_R', required=True, metavar='', help='Directory containing right-aligned reads (FASTA format)')
    
    # Assembly method group
    assembly_group = parser.add_mutually_exclusive_group()
    assembly_group.add_argument('--flye', action='store_true', help='Flye assembly module (default: True)')
    assembly_group.add_argument('--assemble', '-a', action='store_true', help='Assemble using an alternative assembly module')
    parser.set_defaults(flye=True)  # Set default to enable Flye

    # General options
    parser.add_argument('-t', '--threads', type=int, default=20, metavar='', help='Number of threads to use (default: 20)')
    parser.add_argument('--min_overlap', type=int, default=50, metavar='', help='Minimum overlap length (default: 50)')
    parser.add_argument('--error_rate', type=float, default=0.15, metavar='', help='Error rate for assembly (default: 0.15)')
    parser.add_argument('--kmer_size', type=int, default=15, metavar='', help='K-mer size (default: 15)')

    # File paths
    parser.add_argument('-L', '--lgsreads', required=True, metavar='', help='Long-read sequencing data')
    parser.add_argument('-W', '--wgs1', required=True, metavar='', help='Path to WGS reads (read 1)')
    parser.add_argument('-w', '--wgs2', required=True, metavar='', help='Path to WGS reads (read 2)')
    parser.add_argument('-N', '--NextPolish', required=True, metavar='', help='Path to NextPolish tool')

    # Parse arguments
    args = parser.parse_args()

    # configure root logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s"
    )

    # -- Now decide what to run -- 
    # If the user gave neither --flye nor --assemble, default to flye:
    if not args.flye and not args.assemble:
        args.flye = True

    # If they asked for assemble, turn off flye:
    if args.assemble:
        args.flye = False

    dispatch(
      args.dir_IN_L, args.dir_IN_R,
      args.threads, args.flye, args.assemble,
      args.min_overlap, args.error_rate, args.kmer_size
    )

    polish_dirs=['asm_L','asm_R']
    process_asm_merge_fasta_files(
      polish_dirs, "asm_merged_All_chr.fasta",
      args.lgsreads, args.wgs1, args.wgs2,
      args.threads, args.NextPolish
    )

if __name__ == '__main__':
    main()