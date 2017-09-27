from io import StringIO

# TODO: create functions that return non-overlapping indexes describing the type of SNP within the dataframe

import pytest

import pandas as pd

from numpy import isclose

switch_alleles = lambda a: {"A": "T",
                            "T": "A",
                            "C": "G",
                            "G": "C",}[a] # "N": "N"

@pytest.fixture
def df():

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

    return pd.read_table(StringIO(c), sep="\s+")

def _palindrome(m):

    return (m.E1 == m.E2.apply(switch_alleles)) & (m.O1 == m.O2.apply(switch_alleles)) & (m.E1 == m.O1) & (m.E2 == m.O2)


def _non_inferrable_palindrome(m, tolerance):

    return _palindrome(m) & ~_exposure_and_outcome_allele_frequencies_differ_from_fifty_percent(m, tolerance)


def _inferrable_palindrome(m, tolerance):

   return _palindrome(m) & _exposure_and_outcome_allele_frequencies_differ_from_fifty_percent(m, tolerance)



def _exposure_and_outcome_allele_frequencies_differ_from_fifty_percent(m, tolerance):

   outcome_differs = (m.OAF > 0.5 + tolerance) | (m.OAF < 0.5 - tolerance)
   exposure_differs = (m.EAF > 0.5 + tolerance) | (m.EAF < 0.5 - tolerance)

   return outcome_differs & exposure_differs


def _allele_frequencies_similar_within_tolerance(m, tolerance):

    return pd.Series(isclose(m.OAF, m.EAF, tolerance), index=m.index)



def inferrable_palindromic_snp(m, tolerance):

    x = _inferrable_palindrome(m, tolerance)
    return m.loc[x].rsid


def non_inferrable_palindromic_snp(m, tolerance):

    x = _non_inferrable_palindrome(m, tolerance)
    return m.loc[x].rsid


def _fine_snps(m, tolerance):

    similar_af = _allele_frequencies_similar_within_tolerance(m, tolerance)
    not_palindromic = ~_non_inferrable_palindrome(m, tolerance)

    return (m.E1 == m.O1) & (m.E2 == m.O2) & similar_af & not_palindromic


def fine_snps(m, tolerance):

    fine = _fine_snps(m, tolerance)

    return m.loc[fine].rsid


def test_is_fine(df):

    fine = fine_snps(df, 0.08)

    assert fine.reset_index(drop=True).equals(pd.Series(["a", "e", "k"]))


def test_is_non_inferrable_palindromic_snp(df):

    non_inferrable_palindromic = non_inferrable_palindromic_snp(df, 0.08)

    print(non_inferrable_palindromic)

    assert non_inferrable_palindromic.reset_index(drop=True).equals(pd.Series(["f"]))


def is_inferrable_palindromic_snp(m, tolerance):

    x = _inferrable_palindrome(m, tolerance) & ~_fine_snps(m, tolerance)

    return m.loc[x].rsid


def test_is_inferrable_palindromic_snp(df):

    inferrable_palindromic = is_inferrable_palindromic_snp(df, 0.08)

    print(inferrable_palindromic)

    assert inferrable_palindromic.reset_index(drop=True).equals(pd.Series(["d"]))


def _to_swap(m, tolerance):

    return (m.E1 == m.O2) & (m.E2 == m.O1) & ~_palindromic_to_swap(m, tolerance) & ~_to_flip(m, tolerance)


def to_swap(m, tolerance):

    return m.loc[_to_swap(m, tolerance)].rsid


def palindromic_to_swap(m, tolerance):

    return m.loc[_palindromic_to_swap(m, tolerance)].rsid


def _palindromic_to_swap(m, tolerance):

    return (m.E1 == m.O2) & (m.E2 == m.O1) & (m.E1.apply(switch_alleles) == m.O1) & (m.E2.apply(switch_alleles) == m.O2) & ~_exposure_and_outcome_allele_frequencies_differ_from_fifty_percent(m, tolerance)


def test_to_swap(df):

    swap = to_swap(df, 0.08)

    print(swap)

    # swap.reset_index(drop=True)

    assert swap.reset_index(drop=True).equals(pd.Series(["b", "l"]))


def test_palindromic_to_swap(df):

    palindromic_swap = palindromic_to_swap(df, 0.08)

    print(palindromic_swap)
    assert palindromic_swap.reset_index(drop=True).equals(pd.Series(["g"]))


def to_flip(m, tolerance):

    return m.loc[_to_flip].rsid

def _to_flip(m, tolerance):

    return _allele_frequencies_similar_within_tolerance(m, tolerance) & (m.E1.apply(switch_alleles) == m.O1) & (m.E2.apply(switch_alleles) == m.O2)
