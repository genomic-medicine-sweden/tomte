#!/usr/bin/env python

import argparse


def get_rrna_gtf(gtf, oupt):
    lines_rrna = []

    with open(gtf, "r") as gtf_file:
        for line in gtf_file:
            if line.startswith("#"):
                lines_rrna.append(line)

            if "rRNA" in line:
                lines_rrna.append(line)
    with open(oupt, "w") as outfile:
        outfile.write("".join(lines_rrna))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description="""Generate gtf with only rRNA.""",
    )
    parser.add_argument("--gtf", type=str, help="Transcript annotation file in gtf format", required=True)
    parser.add_argument("--output", type=str, help="output gtf file name", required=True)
    args = parser.parse_args()

get_rrna_gtf(args.gtf, args.output)
