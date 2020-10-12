from argparse import ArgumentParser
import pysam
from collections import defaultdict
import numpy as np
import tempfile
import subprocess
import shlex
import os


class PhaseBlock(object):
    def __init__(self, id, start, end, phase, status):
        self.id = id
        self.start = start
        self.end = end
        self.phase = phase
        self.status = status  # biphasic, monophasic, unphased


def main():
    args = get_args()
    bam = pysam.AlignmentFile(args.bam, "rb")
    vcfs_per_chromosome = []
    for chrom in bam.references:  # Iterate over all chromosomes separately
        print(f"Working on chromosome {chrom}")
        phase_blocks = check_phase_blocks(bam, chrom)
        # Need to add in the unphased blocks too by complementing

        variant_files = defaultdict(list)
        for block in phase_blocks:
            tmpbams = make_bams(bam, chrom=chrom, phase_block=block)
            for tmpbam, phase in zip(tmpbams, block.phase):
                tmpvcf = sniffles(tmpbam, status=block.status)
                variant_files[phase].append(tmpvcf)
        H1 = concat_vcf(variant_files['1'])
        H2 = concat_vcf(variant_files['2'])
        chrom_vcf = merge_haplotypes(H1, H2)
        vcfs_per_chromosome.append(chrom_vcf)


def get_args():
    parser = ArgumentParser(description="Use Sniffles on a phased bam to get phased SV calls")
    parser.add_argument("-b", "--bam", "phased bam to perform phased SV calling on")
    parser.add_argument("-v", "--vcf", "output VCF file")
    return parser.parse_args()


def check_phase_blocks(bam, chromosome):
    phase_dict = defaultdict(list)
    coordinate_dict = defaultdict(list)
    for read in bam.fetch(contig=chromosome):
        if read.has_tag('HP'):
            phase_dict[read.get_tag('PS')].append(read.get_tag('HP'))
            coordinate_dict[read.get_tag('PS')].extend([read.reference_start, read.reference_end])
    phase_blocks = []
    for block_identifier in phase_dict.keys():
        if '1' in phase_dict[block_identifier] and '2' in phase_dict[block_identifier]:
            phase_blocks.append(
                PhaseBlock(
                    id=block_identifier,
                    start=np.amin(coordinate_dict[block_identifier]),
                    end=np.amax(coordinate_dict[block_identifier]),
                    phase=[1, 2],
                    status='biphasic')
            )
        else:
            phase_blocks.append(
                PhaseBlock(
                    id=block_identifier,
                    start=np.amin(coordinate_dict[block_identifier]),
                    end=np.amax(coordinate_dict[block_identifier]),
                    phase=[phase_dict[block_identifier][0]],
                    status='monophasic')
            )
        return sorted(phase_blocks, key=lambda x: x.start)


def make_bams(bam, chrom, phase_block):
    tmp_bam_paths = []
    for phase in phase_block.phase:
        handle, tmppath = tempfile.mkstemp(suffix=".bam")
        tmpbam = pysam.AlignmentFile(tmppath, mode='wb', template=bam)
        for read in bam.fetch(contig=chrom, start=phase_block.start, end=phase_block.end):
            if read.get_tag('HP') == phase:
                tmpbam.write(read)
        tmp_bam_paths.append(tmppath)


def sniffles(tmpbam, status):
    handle, tmppath = tempfile.mkstemp(suffix=".vcf")
    subprocess.call(shlex.split(f"sniffles -m {tmpbam} -v {tmppath} --genotype"))
    os.remove(tmpbam)
    if status == 'monophasic':
        return tmppath
    else:
        return filter_vcf(tmppath)


def filter_vcf(tmpvcf):
    pass
    # Also think about removing the older vcf


def concat_vcf(vcfs):
    pass
    # Also think about removing the VCFs


def merge_haplotypes(H1, H2):
    pass
    # Also think about removing the VCFs


if __name__ == '__main__':
    main()
