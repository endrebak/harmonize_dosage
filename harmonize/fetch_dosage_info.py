
from subprocess import check_output
from io import BytesIO

import pandas as pd

def get_header_row(output):

    for i, line in enumerate(output.decode().split("\n")):
        if line.startswith("#CHROM"):
            return i


def fetch_dosage_info(dosage_path):

    command = "bcftools view {}".format(dosage_path)

    output = check_output(command, shell=True)

    header_row = get_header_row(output)

    df = pd.read_table(BytesIO(output), sep="\t", header=header_row, usecols=range(8))

    alt_allele_freq = df.INFO.str.split("AF=|;", expand=True)[1].astype(float)
    alt_allele_freq.name = "AF"

    subset = df["ID REF ALT".split()]

    result = pd.concat([subset, alt_allele_freq], 1)

    return result
