from io import StringIO

# TODO: create functions that return non-overlapping indexes describing the type of SNP within the dataframe

import pytest

import pandas as pd


from harmonize.harmonize import *

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


def test_incompatible(df):

    incom = incompatible(df, 0.08)

    print(incom)

    assert incom.reset_index(drop=True).equals(pd.Series(["c"]))


def test_to_swap_and_flip(df):

    sf = to_swap_and_flip(df, 0.08)

    print(sf)
    assert sf.reset_index(drop=True).equals(pd.Series(["i"]))


def test_to_flip(df):

    flip = to_flip(df, 0.08)

    print(flip)
    assert flip.reset_index(drop=True).equals(pd.Series(["h", "j"]))

def test_to_swap(df):

    swap = to_swap(df, 0.08)

    print(swap)

    # swap.reset_index(drop=True)

    assert swap.reset_index(drop=True).equals(pd.Series(["b", "l"]))


def test_palindromic_to_swap(df):

    palindromic_swap = palindromic_to_swap(df, 0.08)

    print(palindromic_swap)
    assert palindromic_swap.reset_index(drop=True).equals(pd.Series(["g"]))

def test_is_fine(df):

    fine = fine_snps(df, 0.08)

    assert fine.reset_index(drop=True).equals(pd.Series(["a", "e", "k"]))


def test_is_non_inferrable_palindromic_snp(df):

    non_inferrable_palindromic = non_inferrable_palindromic_snp(df, 0.08)

    print(non_inferrable_palindromic)

    assert non_inferrable_palindromic.reset_index(drop=True).equals(pd.Series(["f"]))

def test_is_inferrable_palindromic_snp(df):

    inferrable_palindromic = is_inferrable_palindromic_snp(df, 0.08)

    print(inferrable_palindromic)

    assert inferrable_palindromic.reset_index(drop=True).equals(pd.Series(["d"]))


def test_assign_table(df):

    assign_table(df, 0.08)

    assert 0
