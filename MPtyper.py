import os, sys, click, glob, shutil, subprocess, collections, re, numpy as np


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

def p1_check(gtt, cutoff) : 
    pattern_P1 = r'(P1-2)\(([\w]+),([-]?[\w\.]+)\)'
    p1_type = []
    p1 = re.findall(pattern_P1,gtt)[0]
    if float(p1[2]) == -1 :
        p1_type.append("N_D")
    else :
        if float(p1[2]) >= cutoff :
            p1_type.append("P1-2({0})".format(p1[2]))
        elif float(p1[2]) <= 1 - cutoff :
            p1_type.append("P1-1({0:.2f})".format(1 - float(p1[2])))
        else :
            p1_type.append("P1-1({0:.2f}),P1-2({1})".format(1 - float(p1[2]),p1[2]))

    return(p1_type[0])

def ec_check(gtt, cutoff) : 
    pattern_EC1 = r'(EC1)\(([\w]+),([-]?[\w\.]+)\)'
    pattern_EC2 = r'(EC2)\(([\w]+),([-]?[\w\.]+)\)'
    pattern_EC3 = r'(EC3)\(([\w]+),([-]?[\w\.]+)\)'
    ec1 = re.findall(pattern_EC1,gtt)[0]
    ec2 = re.findall(pattern_EC2,gtt)[0]
    ec = []
    if float(ec1[2]) >= 1 - cutoff :
        ec.append("EC1({0})".format(ec1[2]))
    if float(ec2[2]) >= 1 - cutoff :
        ec.append("EC2({0})".format(ec2[2]))

    if not ec :
        ec.append("N_D")
        
    return(",".join(ec))

def barcode_check(gtt, cutoff) :
    pattern_barcode = r'(Barcode_[\d\.]+)\(([\w]+),([-]?[\w\.]+)\)'
    barcode = re.findall(pattern_barcode, gtt)
    bar = []
    bar_mix = []
    for b in barcode :
        if float(b[2]) >= cutoff :
            bar.append("{0}({1})".format(b[0],b[2]))
        elif float(b[2]) > 1 - cutoff and float(b[2]) < cutoff :
            bar_mix.append("{0}({1})".format(b[0],b[2]))

    if bar_mix :
        bar_mix.extend(bar)
        return(",".join(bar_mix))
    
    elif bar :
        return(bar[-1])
    
    elif not bar :
        bar.append("N_D")
        return(",".join(bar))

def mr_check(gtt) :
    pattern_MR = r'(MR[\w]+)\(([\w]+),([-]?[\w\.]+)\)'
    mr = re.findall(pattern_MR,gtt)
    mut = []
    non_mut = []
    for m in mr :
        if float(m[2]) > 0 :
            mut.append(m[0])
        elif float(m[2]) == 0 :
            non_mut.append(m[0])

    if not mut :
        if len(non_mut) == len(mr) :
            mut.append("non-MUT")
        elif any("2063" in i for i in non_mut) :
            mut.append("Close to non-MUT")
        else :
            mut.append("Close to N_D")
    
    return(",".join(mut))

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
  --cutoff FLOAT     cutoff prability for checks of mixed infection. 
                     [default:0.5, output only single genotype result]
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
@click.option('--cutoff', help='cutoff prability for checks of mixed infection. [default:0.5, output only single genotype results]', default=0.5, type=float)
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
            if d > 0 :
                gt.append("{0}({1},{2:.2f})".format(t, d, d/(a+d)))
            elif a > 0 :
                gt.append("{0}({1},{2:.2f})".format(t, a, 0))
            else :
                gt.append("{0}({1},{2})".format(t, 0, -1))
        gtt.append(';'.join(gt))

    gtt = "\t".join(gtt)

    cutoff = cutoff # Set a cutoff for checks of mixed infection:
    p1 = p1_check(gtt, cutoff)
    ec = ec_check(gtt, cutoff)
    barcode = barcode_check(gtt, cutoff)
    mr =  mr_check(gtt)
    

    gtt_r = p1 + ";" + ec + "\t" + barcode + "\t" + mr # simplified output format
    with open(f"{prefix}.genotypes", "wt") as fout :
        fout.write("#Prefix\tInput\tmapped_count\tmapped_percentage\t{0}\n".format('\t'.join([ os.path.basename(fn).rsplit('.', 1)[0] for fn in snvs ])))
        fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(prefix, ','.join(reads[:2]), mapped_reads, mapped_percentage, gtt))
        sys.stdout.write("{0}\t{1}\n".format(prefix, gtt_r))

    if consensus :
        with open(f"{prefix}.consensus.fa", "wt") as fout :
            for n, s in sequences.items() :
                fout.write(f'>{prefix}__{n} {gtt} lowercase_for_uncertain_bases=(depth<{min_depth} or proportion<{min_cons})\n{s}\n')


if __name__ == '__main__' :
    get_site_info()
