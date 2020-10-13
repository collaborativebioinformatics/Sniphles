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
        self.phase = phase  # List with 1 and/or 2
        self.status = status  # biphasic, monophasic, unphased

    def __repr__(self):
        return f"{str(self.chrom)}:{str(self.start)}-{str(self.end)} => {'-'.join([str(h) for h in self.phase])}"


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
            tmpdvcf = tempfile.mkdtemp(prefix="sniffles")
            phase_blocks = check_phase_blocks(bam, chrom)
            # Adding unphased blocks by complementing
            phase_blocks.extend(get_unphased_blocks(phase_blocks, bam.get_reference_length(chrom)))
            # LOG THE NUMBER OF BIPHASIC, MONOPHASIC AND UNPHASED BLOCKS
            variant_files = defaultdict(list)

            for block in phase_blocks:
                eprint(f"Working on block {block}")
                tmpbams = make_bams(bam, block)
                for tmpbam, phase in zip(tmpbams, block.phase):
                    if tmpbam:
                        cov = get_coverage(tmpbam, block)
                        # XXX (Evaluation needed) Do not attempt to call SVs if coverage of phased block < 10
                        if cov >= 10:
                            tmpvcf = sniffles(tmpdvcf, tmpbam, block.status)
                            variant_files[phase].append(tmpvcf)
                        os.remove(tmpbam)
                        os.remove(tmpbam + '.bai')
            h1_vcf = concat_vcf(variant_files['1'])
            h2_vcf = concat_vcf(variant_files['2'])
            unph_vcf = concat_vcf(variant_files['u'])
            chrom_vcf = merge_haplotypes(h1_vcf, h2_vcf, unph_vcf)
            vcfs_per_chromosome.append(chrom_vcf)
            shutil.rmtree(tmpdvcf)
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
            phase_dict[read.get_tag('PS')].append(read.get_tag('HP'))
            coordinate_dict[read.get_tag('PS')].extend([read.reference_start, read.reference_end])
    phase_blocks = []
    for block_identifier in phase_dict.keys():
        if 1 in phase_dict[block_identifier] and 2 in phase_dict[block_identifier]:
            phase_blocks.append(
                PhaseBlock(
                    id=block_identifier,
                    chrom=chromosome,
                    start=np.amin(coordinate_dict[block_identifier]),
                    end=np.amax(coordinate_dict[block_identifier]),
                    phase=[1, 2],
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


def get_unphased_blocks(phase_blocks, chromosome_end_position):
    """
    [x] implementation done
    [ ] test done

    Returns intervals per chromosome where no phasing information is available.

    Parameters
    ----------
        phase_blocks : PhaseBlock[]
            List of known phase block instances.

        chromosome_end_position : Int
            Index for the end position of a chromosome

    Returns
    -------
        unphased_blocks : PhaseBlock[]
            Intervals in a chromosome where the phase is not known as list of PhaseBlock instances.
    """

    start_positions = sorted([block.start for block in phase_blocks])

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

    chromosome = phase_blocks[0].chrom
    unphased_blocks = [PhaseBlock(
        id='NOID',
        chrom=chromosome,
        start=interval[0],
        end=interval[1],
        phase=['u'],
        status='unphased'
    ) for interval in unphased_intervals if interval[0] != interval[1]]

    return sorted(unphased_blocks, key=lambda x: x.start)


def make_bams(bam, block):
    """
    This function will take
    - the bam file
    - the chromosome we're working on
    - the phase block to isolate

    And produces one or two temporary bam file(s) with the reads of this locus
    if the locus is phased (phase_block.phase is not 'u'/unphased)
    then only take reads assigned to this phase

    If a bam file ends up empty, just return None for that block

    [x] implementation done
    [ ] test done
    """
    tmp_bam_paths = []
    for phase in block.phase:
        reads_in_block = 0
        _, tmppath = tempfile.mkstemp(suffix=".bam")
        tmpbam = pysam.AlignmentFile(tmppath, mode='wb', template=bam)
        if phase == 'u':
            for read in bam.fetch(contig=block.chrom, start=block.start, end=block.end):
                tmpbam.write(read)
                reads_in_block += 1
        else:
            for read in bam.fetch(contig=block.chrom, start=block.start, end=block.end):
                if read.has_tag('HP') and read.get_tag('HP') == phase:
                    tmpbam.write(read)
                    reads_in_block += 1
        tmpbam.close()
        if reads_in_block > 0:
            try:
                pysam.index(tmppath)
            except:
                eprint(
                    f"Problem indexing {tmppath} for {block} when phase is {phase} and the file has {reads_in_block} reads")
                raise()
            tmp_bam_paths.append(tmppath)
        else:
            tmp_bam_paths.append(None)
    return tmp_bam_paths


def get_coverage(tmpbam, block):
    """
    [x] implementation done
    [ ] test done
    """
    tmpdir = tempfile.mkdtemp(prefix="mosdepth")
    _, tmpbed = tempfile.mkstemp(suffix=".bed")
    with open(tmpbed, 'w') as outf:
        outf.write(f"{block.chrom}\t{block.start}\t{block.end}\n")
    subprocess.call(
        shsplit(f"mosdepth -n -x -b {tmpbed} {tmpdir}/{block.chrom}.{block.start} {tmpbam}"))
    cov = np.loadtxt(f"{tmpdir}/{block.chrom}.{block.start}.regions.bed.gz", usecols=3, dtype=float)
    os.remove(tmpbed)
    shutil.rmtree(tmpdir)
    return cov


def sniffles(tmpdvcf, tmpbam, status, support=5):
    """
    [x] implementation done
    [ ] test done

    support: minimal number of supporting reads. Needs to be evaluated. For unphased regions, this number doubles.
    """
    handle, tmppath = tempfile.mkstemp(prefix=tmpdvcf, suffix=".vcf")
    # Used default values in sniffles to filter SVs based on homozygous or heterozygous allelic frequency (AF).
    # Will not attempt to remove calls based on the FILTER field in VCF, which only shows unresovled insertion length other than PASS.
    FNULL = open(os.devnull, 'w')
    if status == "unphased": support *= 2
    subprocess.call(shsplit(
        f"sniffles --genotype --min_homo_af 0.8 --min_het_af 0.3 -s {support} -m {tmpbam} -v {tmppath}"),
        stdout=FNULL,
        stderr=subprocess.STDOUT)
    c = subprocess.Popen(shsplit(f"bcftools sort {tmppath}"), stdout=subprocess.PIPE)
    handle, compressed_vcf = tempfile.mkstemp(suffix=".vcf.gz")
    subprocess.call(shsplit("bgzip -c"), stdin=c.stdout, stdout=handle)
    subprocess.call(shsplit(f"tabix {compressed_vcf}"))
    os.remove(tmppath)
    return compressed_vcf


def concat_vcf(vcfs, output=tempfile.mkstemp(suffix=".vcf")[1]):
    """
    [X] implementation done
    [ ] test done
    """
    if vcfs:
        cmd = f"bcftools concat -a {' '.join(vcfs)} | bcftools sort -o {output}"
        subprocess.check_output(shsplit(cmd))
        # remove temp vcf files
        for vcf in vcfs:
            os.remove(vcf)
        return output
    else:
        return None


def merge_haplotypes(H1, H2):
    """
    [ ] implementation done
    [ ] test done
    """
    # Also think about removing the VCFs


if __name__ == '__main__':
    main()
