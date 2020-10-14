import sys
import logging
import argparse
import pysam
from collections import defaultdict
import numpy as np
import tempfile
import subprocess
from shlex import split as shsplit
from pprint import pformat
import os
import shutil
import itertools
from cyvcf2 import VCF, Writer

# Create a custom logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class PhaseBlock(object):
    def __init__(self, id, chrom, start, end, phase, status):
        self.id = id  # identifier of the phase block in the BAM, coordinate of the first SNV
        self.chrom = chrom
        self.start = start  # first coordinate of a read in the phase block
        self.end = end  # last coordinate of a read in the phase block
        assert isinstance(phase[0], str), "Phase should be a string!"
        assert isinstance(phase[0], str), "Phase should be a string!"
        self.phase = phase  # List with 1 and/or 2
        self.status = status  # biphasic, monophasic, unphased
        self.bams = []
        self.vcfs = {k: None for k in ['1', '2', 'u']}

    def __repr__(self):
        return f"{str(self.chrom)}:{str(self.start)}-{str(self.end)}_H{'-'.join([str(h) for h in self.phase])}"

    def make_bams(self, bam):
        """
        This function will take
        - the bam file

        And produces one or two temporary bam file(s) with the reads of this locus
        if the locus is phased (phase_block.phase is not 'u'/unphased)
        then only take reads assigned to this phase

        If a bam file ends up empty, just return None for that block

        [x] implementation done
        [ ] test done
        """
        tmp_bam_paths = []
        for phase in self.phase:
            reads_in_block = 0
            handle, tmppath = tempfile.mkstemp(suffix=".bam")
            tmpbam = pysam.AlignmentFile(tmppath, mode='wb', template=bam)
            if phase == 'u':
                for read in bam.fetch(contig=self.chrom, start=self.start, end=self.end):
                    tmpbam.write(read)
                    reads_in_block += 1
            else:
                phase = int(phase)
                for read in bam.fetch(contig=self.chrom, start=self.start, end=self.end):
                    if read.has_tag('HP') and read.get_tag('HP') == phase:
                        tmpbam.write(read)
                        reads_in_block += 1
            os.close(handle)
            tmpbam.close()
            if reads_in_block > 0:
                # TODO: IDE COMPLAINS ABOUT INDEX NOT SET IN INIT
                pysam.index(tmppath)
                tmp_bam_paths.append(tmppath)
            else:
                tmp_bam_paths.append(None)
        self.bams = tmp_bam_paths

    def sniffles(self, sample="foo", support=5):
        """
        [x] implementation done
        [ ] test done

        support: minimal number of supporting reads.
        Needs to be evaluated.
        For unphased regions, this number doubles.
        """
        for tmpbam, phase in zip(self.bams, self.phase):
            if tmpbam:
                # cov = get_coverage(tmpbam, block)
                handle_1, tmppath = tempfile.mkstemp(suffix=".vcf")
                # Used default values in sniffles to filter SVs based on homozygous or heterozygous allelic frequency (AF).
                # Will not attempt to remove calls based on the FILTER field in VCF, which only shows unresovled insertion length other than PASS.
                FNULL = open(os.devnull, 'w')
                if self.status == "unphased":
                    support *= 2
                subprocess.call(shsplit(
                    f"sniffles --genotype --min_homo_af 0.8 --min_het_af 0.3 -s {support} -m {tmpbam} -v {tmppath}"),
                    stdout=FNULL,
                    stderr=subprocess.STDOUT)
                handle_2, tmpsamp = tempfile.mkstemp()
                with open(tmpsamp, 'w') as outf:
                    outf.write(f"{sample}")
                c1 = subprocess.Popen(
                    shsplit(f"bcftools reheader -s {tmpsamp} {tmppath}"), stdout=subprocess.PIPE)
                c2 = subprocess.Popen(shsplit(f"bcftools sort"),
                                      stdin=c1.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
                handle_3, compressed_vcf = tempfile.mkstemp(prefix=self.__repr__(),
                                                            suffix=".vcf.gz")
                subprocess.call(shsplit("bgzip -c"), stdin=c2.stdout, stdout=handle_3)
                subprocess.call(shsplit(f"tabix {compressed_vcf}"))
                os.close(handle_1)
                os.close(handle_2)
                os.close(handle_3)
                os.remove(tmppath)
                os.remove(tmpsamp)
                os.remove(tmpbam)
                os.remove(tmpbam + '.bai')
                self.vcfs[phase] = compressed_vcf


def main():
    """
    [ ] implementation done
    [ ] test done
    """
    args = get_args()

    # Setting logging file
    f_handler = logging.FileHandler(os.path.join(os.getcwd(), args.log_file))
    f_handler.setLevel(logging.DEBUG)
    f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    f_handler.setFormatter(f_format)
    logger.addHandler(f_handler)
    logger.info("Used params\n{}".format(pformat(vars(args))))

    bam = pysam.AlignmentFile(args.bam, "rb")
    vcfs_per_chromosome = []
    for chrom_info in bam.get_index_statistics():  # Iterate over all chromosomes separately
        if chrom_info.mapped > 0:  # chrom_info is a namedtuple
            chrom = chrom_info.contig
            eprint(f"Working on chromosome {chrom}")
            phase_blocks = check_phase_blocks(bam, chrom)
            # Adding unphased blocks by complementing
            phase_blocks.extend(get_unphased_blocks(
                phase_blocks, bam.get_reference_length(chrom), chrom))
            # LOG THE NUMBER OF BIPHASIC, MONOPHASIC AND UNPHASED BLOCKS

            for block in phase_blocks:
                eprint(f"Working on block {block}")
                block.make_bams(bam)
                block.sniffles(sample=args.name, support=args.minimum_suport_read)

            eprint(f"Concatenating VCFs of {chrom}, haplotype 1")
            h1_vcf = concat_vcf([ph.vcfs['1'] for ph in phase_blocks if ph.vcfs['1']])
            eprint(f"Concatenating VCFs of {chrom}, haplotype 2")
            h2_vcf = concat_vcf([ph.vcfs['2'] for ph in phase_blocks if ph.vcfs['2']])
            eprint(f"Concatenating VCFs of {chrom}, unphased")
            unph_vcf = concat_vcf([ph.vcfs['u'] for ph in phase_blocks if ph.vcfs['u']])
            eprint(f"Merging haplotypes for {chrom}")
            hbams = make_hap_bams(bam, chrom)
            chrom_vcf = merge_haplotypes(hbams, h1_vcf, h2_vcf, unph_vcf)
            vcfs_per_chromosome.append(chrom_vcf)
    concat_vcf(vcfs_per_chromosome, output=args.vcf)


def get_args():
    """
    [x] implementation done
    [ ] test done
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Use Sniffles on a phased bam to get phased SV calls",
        add_help=True)
    parser.add_argument("-b", "--bam",
                        help="Phased bam to perform phased SV calling on",
                        required=True)
    parser.add_argument("-v", "--vcf",
                        help="output VCF file",
                        required=True)
    parser.add_argument("-l", "--log",
                        help="Log file", dest='log_file', type=str,
                        default="sniphles.log")
    parser.add_argument("-s", "--minimum_suport_read",
                        help="Minimum support read to call SV equals to -s in sniffles",
                        dest='minimum_suport_read',
                        type=int,
                        default=4)
    parser.add_argument("-n", "--name",
                        help="Sample name of output VCF",
                        required=True)
    return parser.parse_args()


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_phase_blocks(bam, chromosome):
    """
    [x] implementation done
    [ ] test done

    TODO: check if regions with >2 blocks exist
    """
    phase_dict = defaultdict(list)
    coordinate_dict = defaultdict(list)
    for read in bam.fetch(contig=chromosome):
        if read.has_tag('HP'):
            phase_dict[read.get_tag('PS')].append(str(read.get_tag('HP')))
            coordinate_dict[read.get_tag('PS')].extend([read.reference_start, read.reference_end])
    phase_blocks = []
    for block_identifier in phase_dict.keys():
        if '1' in phase_dict[block_identifier] and '2' in phase_dict[block_identifier]:
            phase_blocks.append(
                PhaseBlock(
                    id=block_identifier,
                    chrom=chromosome,
                    start=np.amin(coordinate_dict[block_identifier]),
                    end=np.amax(coordinate_dict[block_identifier]),
                    phase=['1', '2'],
                    status='biphasic')
            )
        else:
            phase_blocks.append(
                PhaseBlock(
                    id=block_identifier,
                    chrom=chromosome,
                    start=np.amin(coordinate_dict[block_identifier]),
                    end=np.amax(coordinate_dict[block_identifier]),
                    phase=[phase_dict[block_identifier][0]],
                    status='monophasic')
            )
    return sorted(phase_blocks, key=lambda x: x.start)


def get_unphased_blocks(phase_blocks, chromosome_end_position, chromosome_id):
    """
    Returns intervals per chromosome where no phasing information is available.

    [x] implementation done
    [x] test done

    Parameters
    ----------
        phase_blocks : PhaseBlock[]
            List of known phase block instances.

        chromosome_end_position : int
            Index for the end position of a chromosome

        chromosome_id : str
            ID of the chromosome

    Returns
    -------
        unphased_blocks : PhaseBlock[]
            Intervals in a chromosome where the phase is not known as list of PhaseBlock instances.
    """

    def compare_and_update_phased_blocks(previous, current):
        if previous.start <= current.start <= previous.end:
            """As the input is sorted we check if the start of the current block lies within the previous block.
            If that is the case we extend to the start of the previous interval to span the whole region."""
            previous.end = max(previous.end,
                               current.end)  # merge these and extend the end of the phased block, disregards phases.
            return previous
        else:
            return current

    phase_blocks = list(itertools.accumulate(phase_blocks, compare_and_update_phased_blocks))

    start_positions = sorted(list(set([block.start for block in phase_blocks])))

    unphased_intervals_starts = [0]
    unphased_intervals_ends = []

    for start in start_positions:
        max_end_position = max([block.end for block in phase_blocks if block.start == start])

        # The end positions of a known interval are the start of an unphased region
        unphased_intervals_starts.append(max_end_position)

        # The start positions of known intervals are the end of an unphased region
        unphased_intervals_ends.append(start)

    unphased_intervals_ends.append(chromosome_end_position)

    unphased_intervals = zip(unphased_intervals_starts, unphased_intervals_ends)

    unphased_blocks = [PhaseBlock(
        id='NOID',
        chrom=chromosome_id,
        start=interval[0],
        end=interval[1],
        phase=['u'],
        status='unphased'
    ) for interval in unphased_intervals if interval[0] != interval[1]]

    return sorted(unphased_blocks, key=lambda x: x.start)


def get_coverage(tmpbam, block):
    """
    [x] implementation done
    [ ] test done
    """
    tmpdir = tempfile.mkdtemp(prefix="mosdepth")
    handle, tmpbed = tempfile.mkstemp(suffix=".bed")
    with open(tmpbed, 'w') as outf:
        outf.write(f"{block.chrom}\t{block.start}\t{block.end}\n")
    subprocess.call(
        shsplit(f"mosdepth -n -x -b {tmpbed} {tmpdir}/{block.chrom}.{block.start} {tmpbam}"))
    cov = np.loadtxt(f"{tmpdir}/{block.chrom}.{block.start}.regions.bed.gz", usecols=3, dtype=float)
    os.close(handle)
    os.remove(tmpbed)
    shutil.rmtree(tmpdir)
    return cov


def concat_vcf(vcfs, output=tempfile.mkstemp(suffix=".vcf")[1]):
    """
    [X] implementation done
    [ ] test done
    """
    if vcfs:
        c = subprocess.Popen(
            shsplit(f"bcftools concat -a {' '.join(vcfs)}"), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        subprocess.call(shsplit(f"bcftools sort -o {output}"), stdin=c.stdout, stdout=FNULL, stderr=subprocess.DEVNULL)
        # remove temp vcf files
        # for vcf in vcfs:
        #    os.remove(vcf)
        return output
    else:
        return None


def make_hap_bams(bam, chrom):
    """
    [x] implementation done
    [ ] test done
    Partition reads according to phase; discard unphased
    """
    hap_bams = []
    outs = []
    for h in [0, 1]:
        hap_bams.append(tempfile.mkstemp(suffix=".bam")[1])
        outs.append(pysam.AlignmentFile(hap_bams[h], mode='wb', template=bam))
    for read in bam.fetch(contig=chrom):
        if read.has_tag('HP'):
            if read.get_tag('HP') == '1':
                outs[0].write(read)
            else:
                outs[1].write(read)
    for h in [0, 1]:
        outs[h].close()
        # TODO: IDE COMPLAINS ON UNKNOWN INDEX
        pysam.index(hap_bams[h])
    return hap_bams


def make_header(vcf, name="SAMPLE"):
    """
    [x] implementation done
    [ ] test done
    """
    header = ['##fileformat=VCFv4.1', '##source=sniphles']
    for line in vcf.header_iter():
        if line["HeaderType"] == 'CONTIG':
            header.append('##contig=<ID={}>'.format(line['ID']))
    header.extend([
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##ALT=<ID=BND,Description="Translocation">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End of the structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">',
        '##INFO=<ID=HAPSUPPORT,Number=1,Type=String,Description="Support per haplotype.">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'.format(name)])
    return header


def merge_haplotypes(hbams, h1_vcf, h2_vcf, unph_vcf):
    """
    [x] implementation done
    [ ] test done
    """
    # Also think about removing the VCFs, hbams

    # merging h1/2
    handle1, tmptxt = tempfile.mkstemp(suffix=".txt")
    handle2, mvcf = tempfile.mkstemp(suffix=".vcf")
    np.savetxt(tmptxt, [h1_vcf + h2_vcf], fmt='%s')
    # Parameters explained:
    # maximum allowed distance btwn SVs < 1000 bp
    # do not remove SV even if types are different e.g. INS in h1, DEL in h2
    # do not remove SV based on strand
    # do not estimate distance based on the size of SV
    # minimal size of SV >= 0 bp
    subprocess.call(shsplit(f"SURVIVOR merge {tmptxt} 1000 1 0 0 0 0 {mvcf}"))
    os.close(handle1); os.close(handle2)

    # force calling on h1/2
    hvcfs = []
    for h in [1, 2]:
        hvcfs.append(tempfile.mkstemp(suffix=f".{h}.vcf")[1])
        subprocess.call(shsplit(f"sniffles -m {hbams[h-1]} -v {hvcfs[h-1]} --Ivcf {mvcf}"))
        os.remove(hbams[h-1])
    os.remove(mvcf)

    # merging w/ unph_vcf
    handle1, rawvcf = tempfile.mkstemp(suffix=".vcf")
    handle2, chromvcf = tempfile.mkstemp(suffix=".vcf")
    np.savetxt(tmptxt, hvcfs + [unph_vcf], fmt='%s')
    subprocess.call(shsplit(f"SURVIVOR merge {tmptxt} 10 1 0 0 0 0 {rawvcf}"))
    # XXX 2,3 swapped compared to @wouter haplomerge.py
    allele_dict = {0: 'HOM_REF', 1: 'HET', 2: 'HOM_ALT', 3: 'UNKNOWN'}
    unph_gt_to_str = {1: "1/0", 2: "1/1"}
    ph_gt_to_str = {0: {2: "0|1"}, 2: {0: "1|0", 2: "1|1"}}
    vcf = VCF(rawvcf, gts012=True)
    with open(chromvcf, 'w') as f:
        f.write("\n".join(make_header(vcf)) + "\n")
        for v in vcf:
            if v.gt_types[2] == 3:  # gt ordered by h1_vcf/h2_vcf/unphased_vcf
                if v.gt_types[0] == 3 or v.gt_types[1] == 3:
                    assert False, "at least one SV call in phased region is missing"
                else:
                    # if phased, each genotype has to be HOM not HET
                    if v.gt_types[0] == 1 or v.gt_types[1] == 1:
                        eprint(f"{v.ID} w/ gt {v.gt_types} removed due to HET called on a haplotype")
                        continue
                    elif v.gt_types[0] == 0 and v.gt_types[1] == 0:  # not a variant
                        eprint(f"{v.ID} w/ gt {v.gt_types} removed due to no variant")
                        continue
                    else:
                        gt = ph_gt_to_str[v.gt_types[0]][v.gt_types[1]]
            else:
                if v.gt_types[0] != 3 or v.gt_types[1] != 3:
                    assert False, "duplicate calls from phased and unphased regions"
                else:  # if unphased, HOM and HET are both valid
                    if v.gt_types[2] == 0:  # not a variant
                        eprint(f"{v.ID} w/ gt {v.gt_types} removed due to no variant")
                        continue
                    else:
                        gt = unph_gt_to_str[v.gt_types[2]]

            info = {'SVLEN': v.INFO.get('SVLEN'),
                    'END': v.end,
                    'SVTYPE': v.INFO.get('SVTYPE'),
                    'HAPSUPPORT': '-'.join([allele_dict[gt] for gt in v.gt_types])}
            f.write("{chrom}\t{pos}\t{idf}\t{ref}\t{alt}\t{q}\t{filt}\t{info}\t{form}\t{sam}\n"
                    .format(
                        chrom=v.CHROM,
                        pos=v.start,
                        idf=v.ID,
                        ref=v.REF,
                        alt=','.join(v.ALT),
                        q=v.QUAL or '.',
                        filt='.',
                        info=';'.join(['{}={}'.format(k, v) for k, v in info.items()]),
                        form='GT',
                        sam=gt))
    vcf.close()
    os.close(handle1); os.close(handle2)
    os.remove(rawvcf); os.remove(tmptxt)
    for f in h1_vcf + h2_vcf + unph_vcf + hvcfs + hbams:
        os.remove(f)
    return chromvcf


if __name__ == '__main__':
    main()
