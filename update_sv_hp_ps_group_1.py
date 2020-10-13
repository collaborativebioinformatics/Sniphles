
#!/usr/bin/env python3

"""
This script update vcf file to add both HP haplotag and PS phasing block info fields, It takes as input vcf file, hp, ps.
"""
import argparse
import sys, os
from operator import itemgetter
from collections import Counter

# Python program to print
# green text with red background
#
# from colorama import init
# from termcolor import colored
#
# init()



# Part that processing input arguments
def get_args():
    parser = argparse.ArgumentParser(epilog="%(prog)s version 0.01. use command -h for info.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Phase SVs Using haplotyped reads in tab format',
                                     add_help=True, )
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.01')
    # parser.add_argument('input', help='Input file', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    # parser.add_argument('output', help='Output file', nargs="?", type=argparse.FileType('w'), default=sys.stdout)

    parser.add_argument('input', nargs='?', help="Structural variant vcf file",
                             type=argparse.FileType('r'),
                             default=sys.stdin)
    parser.add_argument('hp', nargs='?', help="tab delimeted read\thp\tps file",
                             type=argparse.FileType('r'))
    parser.add_argument('output', nargs='?', help="Output file, PS and HP will be added.",
                                 type=argparse.FileType('w+'),
                                 default=sys.stdout)

    parser.set_defaults(func=update_vcf)

    # if not argument print help.
    if len(sys.argv) == 1 and  sys.stdin.isatty():  # sys.stdin.isatty() returns false if there's something in stdin
         parser.print_help(sys.stderr)
         sys.exit(1)

    args = parser.parse_args()


    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()

def update_vcf(args):
    # check if the input from stdin
    if not sys.stdin.isatty(): # there is nothing in the stdin
        if args.input.name.endswith("gz"):
            import gzip
            myfile = gzip.open(args.input.name, 'rt') # t is not a must normally it is default.
        else:
            myfile = args.input
    else:
        myfile = args.input

    # read the Haplotyped reads file as dictionary
    hp_dic = {}
    with args.hp as hp_in:
        for line in hp_in:
            id, hp, ps = line.split()
            hp_dic[id] = [hp.rsplit(":", 1)[-1], ps.rsplit(":", 1)[-1]] # read -> [hp, ps]


    with myfile as data_in, args.output as data_out:
        for line in data_in:
            reads = []
            if line.startswith('##'):
                data_out.write(line)
            elif line.startswith("#"):
                # data_out.write("##INFO=<ID=HP,Number=1,Type=Integer,Description=\"Haplotype identifier\">\n")
                data_out.write("##INFO=<ID=CONFLICT,Number=.,Type=Integer,Description=\"The Phase is conflict or not\">\n")
                data_out.write("##INFO=<ID=HP_RATIO,Number=.,Type=String,Description=\"Phase Ratio of 1\">\n")
                data_out.write("##FORMAT=<ID=PS,Number=.,Type=Integer,Description=\"Phase set identifier\">\n")
                data_out.write(line)
            else:
                line_split = line.split()
                #1)Adding PS field to FORMAT column for this cases 1|0,0|1,1|1,0|0,.|. for all -1.YES
                #2)Lets try to analyze homozygous variant to cover uniparental disomy cases
                #3)REF_NO_CONFLICT,HMZ_NO_CONFLICT,HET_NO_CONFLICT,HET_SNP_ALLELE_CONFLICT,HET_SNP_MISSING????
                #Binary number,00,01,02
                #4)Calculate ratio for HET.Count of 1's/Total=HP_SV_READ_RATIO_1.Count of 2's/Total=HP_SV_READ_RATIO_2.YES
                if line_split[-1].split(":", 1)[0] == "1/1" or line_split[-1].split(":", 1)[0] == "0/0" or line_split[-1].split(":", 1)[0] == "./.":  # no gt to phase
                    line_split[-2] = "{}:{}".format(line_split[-2], "PS")
                    line_split[-1] = "{}:{}".format(line_split[-1], "-1")
                    data_out.write("{}\n".format("\t".join(line_split)))
                elif line_split[-1].split(":", 1)[0] in ["0/1", "1/0"]:
                    reads = [i for i in line_split[7].split(";")  if i.startswith("RNAMES")][0].split("=",1)[-1].split(",")
                    myvalues = list(map(hp_dic.get, reads))  # list of lists first element id hp second is ps or None in case there are no reads with hp and ps to support this sv
                    # If any value not None
                    total_count = len(myvalues)
                    count_of_1s = 0
                    count_of_2s = 0
                    count_of_None = 0
                    if any(myvalues): # any value is not none
                        ps_dict = categorize_ps(myvalues)
                        if 0 in list(ps_dict.values()): # means that the hp is conflicting do not update anything and add flag that is is conflicting
                            for i in myvalues:
                                if i is not None:
                                    if i[0] == '1':
                                        count_of_1s+=1
                                    elif i[0] == '2':
                                        count_of_2s+=1
                                else:
                                    count_of_None=+1
                            ratios = count_of_1s,count_of_2s,count_of_None
                            line_split[7] = "{info};CONFLICT={conflict};HP_RATIO={hp_ratio}".format(info=line_split[7], conflict=1, hp_ratio=ratios)

                            line_split[-2] = "{}:{}".format(line_split[-2], "PS")
                            line_split[-1] = "{}:{}".format(line_split[-1], ",".join(ps_dict.keys()))
                            data_out.write("{}\n".format("\t".join(line_split)))
                        else: # update the gt field and ps to sv

                            line_split[7] = "{info};CONFLICT={conflict}".format(info=line_split[7], conflict=0)
                            line_split[-2] = "{}:{}".format(line_split[-2], "PS")
                            # if values are negative then it is hp=1 1|0 else it is hp2 0|1
                            # line_split[-1] = line_split[-1].replace("/", "|")
                            hp_new_value = line_split[-1].split(':')
                            if list(ps_dict.values())[0] < 1: # haplotype 1
                                hp_new_value[0] = "1|0"
                            else:
                                hp_new_value[0] = "0|1"
                            hp_new_value = ":".join(hp_new_value)
                            line_split[-1] = "{}:{}".format(hp_new_value, ",".join(ps_dict.keys()))
                            data_out.write("{}\n".format("\t".join(line_split)))
                    else: # all are none
                        line_split[7] = "{info};CONFLICT=2".format(info=line_split[7])
                        data_out.write("{}\n".format("\t".join(line_split)))


# Test case [['1', '23200'], ['2', '23200'], ['2', '23200'], ['1', '23200'], ['2', '23200'], ['2', '23200'], ['1', '23200'], ['2', '23200'], ['1', '23200'], ['2', '23200'], ['1', '23200'], ['1', '23200'], ['2', '23200'], ['2', '23200'], ['2', '23200']]
def categorize_ps(myvalues):
    ps_dict = {}
    for i in myvalues:
        if not i:
            continue
        ps = i[1]
        hp = int(i[0])
        if ps in ps_dict:
            if hp == 1:
                if ps_dict[ps] < 0 :  # if we put =  it will use the paralment decision
                    ps_dict[ps] = ps_dict[ps] - 1
                else: #conflict
                    ps_dict[ps] = 0

            else: # means that it is haplotype 2 hp=2
                if ps_dict[ps] > 0: # if we put =  it will use the paralment decision
                    ps_dict[ps] = ps_dict[ps] + 1
                else: #conflict
                    ps_dict[ps] = 0
        else:
            if hp == 1:
                ps_dict[ps] = -1
            else:
                ps_dict[ps] = 1
    return ps_dict


def main():
    args = get_args()



if __name__ == "__main__":
    main()
