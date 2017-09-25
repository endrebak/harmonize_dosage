from io import StringIO

import pytest

import pandas as pd

from numpy import isclose


@pytest.fixture
def palindromic_table():

    c = """rsid   A1 A2 REF ALT       AF      MAF  FreqA1
rs1   C  G   C   G  0.45  0.55 0.01
rs2   T  A   T   A  0.69802  0.30198  0.2850"""

    return pd.read_table(StringIO(c), sep="\s+", index_col="rsid")


@pytest.fixture
def expected_result():
    s = StringIO("""rsid Inferable
rs1 True
rs2 False""")
    return pd.read_table(s, sep=" ", header=0, index_col=0)


def non_inferable_palindromes(pt):

    switch_alleles = lambda a: {"A": "T",
                                "T": "A",
                                "C": "G",
                                "G": "C",
                                "N": "N"}[a]

    palindromic = (pt.A1.apply(switch_alleles) == pt.A2)
    palindromic_but_inferrable = isclose(pt.FreqA1, pt.MAF, .08)

    df = pd.DataFrame(~palindromic_but_inferrable, index=pt.index, columns=["Inferable"])

    return df


def test_palindrome(palindromic_table, expected_result):

    palindromic_but_not_inferrable = non_inferable_palindromes(palindromic_table)

    print(palindromic_but_not_inferrable)
    print(expected_result)

    assert palindromic_but_not_inferrable.equals(expected_result)


@pytest.fixture
def ambiguous_table():

# exposure effect = 0.5
# effect allele = A
# other allele = G

# outcome effect = -0.05
# effect allele = A
# other allele = C



    c = """rsid   A1 A2 REF ALT       AF      MAF  FreqA1
rs1   C  G   A   C  0.45  0.55 0.01
rs2   C  G   A   C  0.45  0.55 0.01
rs3   T  A   T   G  0.69802  0.30198  0.2850"""

    return pd.read_table(StringIO(c), sep="\s+", index_col="rsid")


def test_ambiguous(pt):

    """
    """
