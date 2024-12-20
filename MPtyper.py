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
@click.option('-b', '--bam', help='flag to keep intermediate BAM file', default=False, is_flag=True)
def get_site_info(database, reads, prefix, consensus, bam, min_depth, min_cons) :
    ref, snvs = check_database(database)
    read_list = check_reads(reads)

    map_cmd = f"{minimap2} -ax sr {ref} {read_list}|{samtools} view -F 4 -h|{samtools} sort -@8 -m 10G -O bam -l 0 -o {prefix}.bam"
    _ = subprocess.Popen(map_cmd, stderr=subprocess.PIPE, shell=True).communicate()

    sam_cmd = f"{samtools} mpileup -AB -aaa {prefix}.bam"
    p = subprocess.Popen(sam_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    
    depths = collections.defaultdict(list)
    for line in p.stdout :
        if not line :
            raise SystemExit(f"{prefix}\tERROR: No matched reads!")
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
            if d > a :
                gt.append("{0}({1:.2f})".format(t, d/(a+d)))
            elif d < a :
                gt.append("{0}({1:.2f})".format(t, -a/(a+d)))
            elif d == a and a != 0:
                gt.append("{0}({1:.2f})".format(t, 0.5))
            else:
                gt.append("{0}({1})".format(t, "N. D."))
        gtt.append(','.join(gt))
    gtt = "\t".join(gtt)
    with open(f"{prefix}.genotypes", "wt") as fout :
        fout.write("#Prefix\tReads\t{0}\n".format('\t'.join([ os.path.basename(fn).rsplit('.', 1)[0] for fn in snvs ])))
        fout.write("{0}\t{1}\t{2}\n".format(prefix, ','.join(reads[:2]), gtt))
        sys.stdout.write("{0}\t{1}\t{2}\n".format(prefix, ','.join(reads[:2]), gtt))

    if consensus :
        with open(f"{prefix}.consensus.fa", "wt") as fout :
            for n, s in sequences.items() :
                fout.write(f'>{prefix}__{n} {gtt} lowercase_for_uncertain_bases=(depth<{min_depth} or proportion<{min_cons})\n{s}\n')


if __name__ == '__main__' :
    get_site_info()
