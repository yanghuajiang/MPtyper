import os, sys, click, glob, shutil, subprocess, collections, re
import pandas as pd
import numpy as np


def check_database(database):
    if not os.path.isdir(database) :
        database = os.path.join(os.path.dirname(__file__), database)
        if not os.path.isdir(database) :
            raise FileNotFoundError(f"Database dir not found: {database}")
        
    ref = os.path.join(database, 'reference.fa')
    if not os.path.isfile(ref) :
        raise FileNotFoundError(f"Reference genome not found in {database}")
    
    snvs = sorted(glob.glob(os.path.join(database, "*.def")))
    if not os.path.isfile(snvs[0]) :
        raise FileNotFoundError(f"SNVs file not found in {database}")
    
    return ref, snvs


def check_reads(reads) :
    reads_list = []
    for read in reads:
        if os.path.isfile(read) :
            reads_list.append(read)
        else :
            read = os.path.jion(os.path.dirname(__file__), read)
            if os.path.isfile(reads) :
                reads_list.append(read)
            else :
                raise FileNotFoundError(f"Reads not found: {read}")
    return " ".join(reads_list)


def result_df(gtt, pattern) :
    result = re.findall(pattern, gtt)
    df = pd.DataFrame(columns=["type","anc","der"])
    for r in result :
        t, a, d = r[0], int(r[1]), int(r[2])
        df = pd.concat(
            [df,
             pd.DataFrame([[t,a,d]], columns=df.columns)],
             ignore_index=True
        )
    df["anc"] = df["anc"].astype(int)
    df["der"] = df["der"].astype(int)
    return(df)

def p1_check(gtt, cutoff) : 
    pattern_P1 = r'(P1-2)\(([\w]+),([-]?[\w\.]+)\)'
    p1 = re.findall(pattern_P1,gtt)[0]
    _, a, d = p1[0], int(p1[1]), int(p1[2])
    if max(a, d) < cutoff :
        return("P1(N_D)")
    else :
        if a > d :
            return("P1-1({0:.2f})".format(a / (a + d)))
        elif d <= a :
            return("P1-2({0:.2f})".format(d / (a + d)))


def ec_check(gtt, cutoff) :
    pattern_EC = r'(EC[1-9])\(([\w]+),([-]?[\w\.]+)\)'  
    df_ec = result_df(gtt, pattern_EC)
    if df_ec[["anc","der"]].max().max() < cutoff :
        return("EC(N_D)")
    
    else :
        EC = []
        for ec in df_ec.itertuples() :
            if max(ec.anc, ec.der) < cutoff :
                continue
            elif ec.anc <= ec.der :
                EC.append("{0}({1})".format(ec.type, ec.der / (ec.der + ec.anc)))
        
        if not EC :
            return("")
        else :
            return(",".join(EC))
        

def barcode_check(gtt, cutoff) :
    pattern_barcode = r'(Barcode_[\d\.]+)\(([\w]+),([-]?[\w\.]+)\)'
    barcode = re.findall(pattern_barcode, gtt)
    bar = []
    for b in barcode :
        t, a, d = b[0], int(b[1]), int(b[2])
        if max(a, d) < cutoff :
            continue
        elif a <= d :
            bar.append("{0}({1})".format(t, d / (a + d)))

    if not bar :
        return("Barcode(N_D)")
    else:
        return(bar[-1])


def mr_check(gtt, cutoff) :
    pattern_MR = r'(MR[\w]+)\(([\w]+),([-]?[\w\.]+)\)'
    df = result_df(gtt, pattern_MR)

    if df[["anc","der"]].max().max() < cutoff :
        return("MR(N_D)")
        
    elif df["der"].sum() > 0 :
        row = df["der"].idxmax()
        return("{0}({1})".format(df["mr"][row], df["der"].max()/df["der"].sum()))
    
    elif df[df['mr'].str.contains('2063', na=False)]['anc'].sum() > 0 : # Position 2063 represents the most prevalent resistant mutation site (>95%).
        return("MR(non-MUT)")
    
    else:
        return("MR(N_D)")
    

class HelpfulCmd(click.Command):
    def format_help(self, ctx, formatter):
        dbs = glob.glob(os.path.join(os.path.dirname(__file__), "*"))
        dbs = [os.path.basename(fn) for fn in dbs if not fn.endswith('.py')]
        click.echo('''Usage: MPtyper.py [OPTIONS]

Options:                 
  Internally available: {0}
                   
  -db, --database TEXT   dirname for the database. [default: MP]
  -o, --prefix TEXT  prefix for the outputs.  [required]
  -r, --reads TEXT   files for short reads, can be specified at most twice. 
                     [required]
  -c, --consensus    flag to generate consensus sequences. (for phylogenetic
                     analysis)
  --min_depth FLOAT  minimum read depth for high quality bases. [only for
                     consensus sequences, default:3]
  --min_cons FLOAT   minimum proportion of consensus reads for high quality
                     bases. [only for consensus sequences, default:0.8]
  --cutoff FLOAT     cutoff of total depth in type-specific sites for validation of genotypes. 
                     [default:1, to assign genotypes for low-quality data]
  -b, --bam          flag to keep intermediate BAM file.
  --help             Show this message and exit.
'''.format(','.join(dbs)))


minimap2 = shutil.which('minimap2')
samtools = shutil.which('samtools')

@click.command(cls=HelpfulCmd)
@click.option('-o', '--prefix', help='prefix for the outputs.', required=True)
@click.option('-r', '--reads', help='files for short reads, can be specified at most twice.', required=True, multiple=True)
@click.option('-db', '--database', help='dirname for the database. [default: MP]', default='MP')
@click.option('-c', '--consensus', help='flag to generate consensus sequences (for phylogenetic analysis)', default=False, is_flag=True)
@click.option('--min_depth', help='minimum read depth for high quality bases. [only for consensus sequences, default:3]', default=3, type=float)
@click.option('--min_cons', help='minimum proportion of consensus reads for high quality bases. [only for consensus sequences, default:0.8]', default=0.8, type=float)
@click.option('--cutoff', help='cutoff of total depth in type-specific sites for validation of genotypes. [default:1, to assign genotypes for low-quality data]', default=1, type=int)
@click.option('-b', '--bam', help='flag to keep intermediate BAM file', default=False, is_flag=True)
def get_site_info(database, reads, prefix, consensus, bam, min_depth, min_cons, cutoff) :
    ref, snvs = check_database(database)
    read_list = check_reads(reads)

    map_cmd = f"{minimap2} -ax sr {ref} {read_list}|{samtools} sort -@8 -m 10G -O bam -l 0 -o {prefix}.bam"
    _ = subprocess.Popen(map_cmd, stderr=subprocess.PIPE, shell=True).communicate()

    map_info = subprocess.getoutput(f"samtools flagstat {prefix}.bam")
    mapped_reads = "0"
    mapped_percentage = "0"
    for line in map_info.splitlines() :
        if "mapped" in line :
            line = line.strip().split()
            if line[0] == "0" :
                if not bam :
                    os.unlink(f"{prefix}.bam")
                raise SystemExit(f"{prefix}\tERROR: No matched reads!")
            mapped_reads = line[0]
            mapped_percentage = line[4].replace("(","")
            break

    sam_cmd = f"{samtools} view -h -F 4 {prefix}.bam|{samtools} mpileup -AB -aaa -"
    p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, shell=True)
    depths = collections.defaultdict(list)

    for line in p.stdout :
        p = line.strip().split('\t')
        bases = re.sub(r'\$', '', re.sub(r'\^.', '', p[4])).upper()
        bases = collections.Counter(''.join([bb[int(n):] for n, bb in re.findall(r"[+-](\d+)([^+-]+)", "+0"+bases)]))
        depths[p[0]].append([bases.get('A', 0), bases.get('C', 0), bases.get('G', 0), bases.get('T', 0)])

    if not bam :
        os.unlink(f"{prefix}.bam")

    sequences = {}
    for n, c in depths.items() :
        depth = np.sum(c, 1)
        max_d = np.max(c, 1)
        q = np.ones(depth.size, dtype=bool)
        q[depth < np.max([np.mean(depth)*0., min_depth])] = 0
        q[max_d < min_cons * depth] = 0
        s = np.array(["A", "C", "G", "T"])[np.argmax(c, 1)]
        s[max_d == 0] = "-"
        s[q == 0] = np.vectorize(lambda x:x.lower())(s[q == 0])
        sequences[n] = ''.join(s)

    gtt = []
    for snv_file in snvs :
        genotypes = collections.OrderedDict()
        with open(snv_file, 'rt') as fin :
            for line in fin :
                p = line.strip().split()
                anc, der = p[3][0].upper(), p[3][-1].upper()
                cont, site = p[1], int(p[2])-1
                anc_i, der_i = {"A":0, "C":1, "G":2, "T":3}[anc], {"A":0, "C":1, "G":2, "T":3}[der]
                anc_d, der_d = depths[cont][site][anc_i], depths[cont][site][der_i]
                if p[0] not in genotypes :
                    genotypes[p[0]] = [anc_d, der_d]
                else :
                    genotypes[p[0]][0], genotypes[p[0]][1] = genotypes[p[0]][0] + anc_d, genotypes[p[0]][1] + der_d

        gt = []
        for t, (a, d) in genotypes.items() :
            gt.append("{0}({1},{2})".format(t,a,d))
        gtt.append(';'.join(gt))

    gtt = "\t".join(gtt)

    cutoff = cutoff 
    p1 = p1_check(gtt, cutoff)
    ec = ec_check(gtt, cutoff)
    barcode = barcode_check(gtt, cutoff)
    mr =  mr_check(gtt, 1) # Resistant type id comfirmed even if there is only one read support.
    

    gtt_r = p1 + ";" + ec + "\t" + barcode + "\t" + mr # simplified format
    with open(f"{prefix}.genotypes", "wt") as fout :
        fout.write("#Prefix\tInput\tmapped_count\tmapped_percentage\t{0}\n".format('\t'.join([ os.path.basename(fn).rsplit('.', 1)[0] for fn in snvs ])))
        fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(prefix, ','.join(reads[:2]), mapped_reads, mapped_percentage, gtt_r))
        sys.stdout.write("{0}\t{1}\n".format(prefix, gtt_r))

    if consensus :
        with open(f"{prefix}.consensus.fa", "wt") as fout :
            for n, s in sequences.items() :
                fout.write(f'>{prefix}__{n} {gtt} lowercase_for_uncertain_bases=(depth<{min_depth} or proportion<{min_cons})\n{s}\n')


if __name__ == '__main__' :
    get_site_info()


