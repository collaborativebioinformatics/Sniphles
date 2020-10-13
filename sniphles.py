import sys
import argparse
import pysam
from collections import defaultdict
import numpy as np
import tempfile
import subprocess
import shlex
import os
from cyvcf2 import VCF, Writer


class PhaseBlock(object):
    def __init__(self, id, start, end, phase, status):
        self.id = id  # identifier of the phase block in the BAM
        self.start = start  # first coordinate of a read in the phase block
        self.end = end  # last coordinate of a read in the phase block
        self.phase = phase  # '1' or '2'
        self.status = status  # biphasic, monophasic, unphased


def main():
    """
    [ ] implementation done
    [ ] test done
    """
    args = get_args()
    bam = pysam.AlignmentFile(args.bam, "rb")
    vcfs_per_chromosome = []
    for chrom in bam.references:  # Iterate over all chromosomes separately
        eprint(f"Working on chromosome {chrom}")
        phase_blocks = check_phase_blocks(bam, chrom)
        phase_blocks.extend(get_unphased_blocks(phase_blocks, 0, bam.get_reference_length(
            chrom)))  # Adding unphased blocks by complementing
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
    """
    [x] implementation done
    [ ] test done
    """
    parser = argparse.ArgumentParser(
                formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                description="Use Sniffles on a phased bam to get phased SV calls",
                add_help=True)
    parser.add_argument("-b", "--bam", help="Phased bam to perform phased SV calling on")
    parser.add_argument("-v", "--vcf", help="output VCF file")

    if len(sys.argv) == 1:
         parser.print_help(sys.stderr)
         sys.exit(1)

    return parser.parse_args()


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_phase_blocks(bam, chromosome):
    """
    [x] implementation done
    [ ] test done
    """
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


def get_unphased_blocks(phase_blocks, chromosome_start_position, chromosome_end_position):
    """
    [x] implementation done
    [ ] test done

    Returns intervals per chromosome where no phasing information is available.

    Parameters
    ----------
        phase_blocks : PhaseBlock[]
            List of known phase block instances.

        chromosome_start_position : Int
            Index for the start position of a chromosome

        chromosome_end_position : Int
            Index for the end position of a chromosome

    Returns
    -------
        unphased_blocks : PhaseBlock[]
            Intervals in a chromosome where the phase is not known as list of PhaseBlock instances.
    """

    start_positions = sorted([block.start for block in phase_blocks])

    unphased_intervals_starts = [chromosome_start_position]
    unphased_intervals_ends = []

    for start in start_positions:
        max_end_position = max([block.end for block in phase_blocks if block.start == start])

        unphased_intervals_starts.append(max_end_position)  # The end positions of a known interval are the start of
        # an unphased region

        unphased_intervals_ends.append(start)  # The start positions of known intervals are the end of an unphased
        # region

    unphased_intervals_ends.append(chromosome_end_position)

    unphased_intervals = zip(unphased_intervals_starts, unphased_intervals_ends)

    unphased_blocks = [PhaseBlock(
        id='NOID',
        start=interval[0],
        end=interval[1],
        phase=[],
        status='unphased'
    ) for interval in unphased_intervals if interval[0] != interval[1]]

    return sorted(unphased_blocks, key=lambda x: x.start)


def make_bams(bam, chrom, phase_block):
    """
    [x] implementation done
    [ ] test done
    """
    tmp_bam_paths = []
    for phase in phase_block.phase:
        handle, tmppath = tempfile.mkstemp(suffix=".bam")
        tmpbam = pysam.AlignmentFile(tmppath, mode='wb', template=bam)
        for read in bam.fetch(contig=chrom, start=phase_block.start, end=phase_block.end):
            if read.get_tag('HP') == phase:
                tmpbam.write(read)
        tmp_bam_paths.append(tmppath)
    return tmp_bam_paths


def sniffles(tmpbam, status):
    """
    [x] implementation done
    [ ] test done
    """
    handle, tmppath = tempfile.mkstemp(suffix=".vcf")
    subprocess.call(shlex.split(
        f"sniffles -m {tmpbam} -v {tmppath} --genotype"))  # I would set minumum coverage to 2 then we can filter later -s 2 or get the coverage per bam
    if status == 'monophasic':
        return tmppath
    else:
        return filter_vcf(tmppath)


def filter_vcf(tmpvcf):
    """
    [x] implementation done
    [ ] test done
    """
    vcf = VCF(tmpvcf)
    handle, tmppath = tempfile.mkstemp(suffix=".vcf")
    w = Writer(tmppath, vcf)
    for variant in vcf:
        pass  # If the variant is not supported by almost all of the reads, then remove it
    os.remove(tmpvcf)
    return tmppath


def concat_vcf(vcfs):
    """
    [ ] implementation done
    [ ] test done
    """
    pass
    # Also think about removing the VCFs


def merge_haplotypes(H1, H2):
    """
    [ ] implementation done
    [ ] test done
    """
    pass
    # Also think about removing the VCFs


if __name__ == '__main__':
    main()
