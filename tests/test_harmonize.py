from io import StringIO

TODO: create functions that return non-overlapping indexes describing the type of SNP within the dataframe

import pytest

import pandas as pd

from numpy import isclose

@pytest.fixture
def final_fixture():

    """Thanks to Gibran Hemani for creating this test set."""

    c = """rsid   E1 E2 O1 O2       OAF   EAF note
a  A  G  A  G        0.5  0.5                 fine
b  A  G  G  A        0.5  0.5              to_swap
c  A  G  T  A        0.5  0.5         incompatible
d  A  T  A  T        0.2  0.8     to_flip_and_swap
e  A  T  A  T        0.2  0.2                 fine
f  A  T  A  T        0.5  0.5          palindromic
g  A  T  T  A        0.5  0.5  palindromic_to_swap
h  C  G  G  C        0.2  0.2              to_flip
i  G  T  A  C        0.5  0.5     to_flip_and_swap
j  G  T  C  A        0.5  0.5              to_flip
k  A  G  A  G        0.2  0.2                 fine
l  T  G  G  T        0.2  0.8              to_swap"""

    return pd.read_table(StringIO(c), sep="\s+", index_col="rsid")


@pytest.fixture
def pt(): # palindromic_table

    c = """rsid   E1 E2 O1 O2       OAF  EAF
rs1   C  G   C   G  0.45  0.01
rs2   T  A   T   A  0.69802  0.2850"""

    return pd.read_table(StringIO(c), sep="\s+", index_col="rsid")


@pytest.fixture
def expected_result_palindromes():
    s = StringIO("""rsid NonInferable
rs1 True
rs2 False""")
    return pd.read_table(s, sep=" ", header=0, index_col=0)

switch_alleles = lambda a: {"A": "T",
                            "T": "A",
                            "C": "G",
                            "G": "C",
                            "N": "N"}[a]

def is_palindromic(pt):

    return (pt.E1.apply(switch_alleles) == pt.E2)


def non_inferable_palindromes(pt, threshold):

    palindromic = is_palindromic(pt)
    palindromic_df = pt[palindromic]

    # Check that the EAF is close to OAF
    inferrable_and_right = isclose(palindromic_df.EAF, palindromic_df.OAF, threshold)
    # Check that the EAF is close to OAF
    inferrable_but_reverse = isclose(palindromic_df.EAF, 1 - palindromic_df.OAF, threshold)

    print("inferrable_and_right", inferrable_and_right)
    print("inferrable_but_reverse", inferrable_but_reverse)

    palindromic_but_not_inferrable = palindromic & ~(inferrable_and_right | inferrable_but_reverse)

    df = pd.DataFrame(palindromic_but_not_inferrable, index=palindromic_df.index, columns=["NonInferable"])

    return df


def inferable_palindromes(pt, threshold):

    palindromic = is_palindromic(pt)
    inferrable = isclose(pt.EAF, pt.OMAF, threshold)

    palindromic_and_inferrable = palindromic & inferrable

    df = pd.DataFrame(palindromic_and_inferrable, index=pt.index, columns=["NonInferable"])

    return df


def test_palindrome(pt, expected_result_palindromes):

    palindromic_but_not_inferrable = non_inferable_palindromes(pt, 0.08)

    print(palindromic_but_not_inferrable)
    print(expected_result_palindromes)

    assert palindromic_but_not_inferrable.equals(expected_result_palindromes)


@pytest.fixture
def at(): # ambiguous_table

    # rsid   A1 A2 REF ALT       AF      MAF  FreqA1

    c = """rsid   E1 E2 O1 O2       OAF      EAF
rs4   A  G   A   C  0.20 0.50
rs5   A  G   T   C  0.20 0.50"""

    return pd.read_table(StringIO(c), sep="\s+", index_col="rsid")


def ambiguous_snps(at):

    equal1, switched1 = (at.E1 == at.O1), (at.E1 == at.O1.apply(switch_alleles))
    equal2, switched2 = (at.E2 == at.O2), (at.E2 == at.O2.apply(switch_alleles))

    snp_ambiguous = ~((equal1 == equal2) | (switched1 == switched2)).to_frame()
    snp_ambiguous.columns = ["Ambiguous"]

    return snp_ambiguous


@pytest.fixture
def expected_result_ambiguous():
    s = StringIO("""rsid Ambiguous
rs4 True
rs5 False""")
    return pd.read_table(s, sep=" ", header=0, index_col=0)


def test_ambiguous(at, expected_result_ambiguous):

    result = ambiguous_snps(at)

    print(result)
    print(expected_result_ambiguous)

    assert result.equals(expected_result_ambiguous)


def switch_reverse_palindromes():

    palindromic_but_inferable
