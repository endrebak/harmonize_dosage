from io import StringIO

import pytest

import pandas as pd


from harmonize.harmonize import *


@pytest.fixture
def df2():

    """Thanks to Gibran Hemani for creating this test set."""

    c = """rsid   E1 E2 O1 O2       OAF   EAF
a  A  G  G  G        0.5  0.5
b  G  G  G  G        0.5  0.5
c  T  C  A  G        0.5  0.5
d  C  T  A  G        0.5  0.5
e  C  T  A  G        0.2  0.4
g  A  T  A  T        0.2  0.8
h  A  T  A  T        0.2  0.87
i  A  T  A  T        0.2  0.88
j  A  T  A  T        0.2  0.89
k  A  T  A  T        0.43  0.57
l  A  T  A  T        0.42  0.58
m  A  T  A  T        0.41  0.59"""
    return pd.read_table(StringIO(c), sep="\s+")


@pytest.fixture
def expected_result_df2():

    c = """rsid   E1 E2 O1 O2       OAF   EAF  SNPType
a  A  G  G  G        0.5  0.5  incompatible
b  G  G  G  G        0.5  0.5  incompatible
c  T  C  A  G        0.5  0.5  flipped
d  C  T  A  G        0.5  0.5  flipped_and_swapped
e  C  T  A  G        0.2  0.4  flipped_and_swapped
g  A  T  A  T        0.2  0.8  inferrable_palindromic
h  A  T  A  T        0.2  0.87 inferrable_palindromic
i  A  T  A  T        0.2  0.88 inferrable_palindromic
j  A  T  A  T        0.2  0.89 incompatible
k  A  T  A  T        0.43  0.57 noninferrable_palindromic
l  A  T  A  T        0.42  0.58 noninferrable_palindromic
m  A  T  A  T        0.41  0.59 inferrable_palindromic"""

    return pd.read_table(StringIO(c), sep="\s+")


@pytest.fixture
def df1():

    """Thanks to Gibran Hemani for creating this test set."""

    c = """rsid   E1 E2 O1 O2       OAF   EAF note
a  A  G  A  G        0.5  0.5                 fine
b  A  G  G  A        0.5  0.5              swapped
c  A  G  T  A        0.5  0.5         incompatible
d  A  T  A  T        0.2  0.8     flipped_and_swapped
e  A  T  A  T        0.2  0.2                 fine
f  A  T  A  T        0.5  0.5          palindromic
g  A  T  T  A        0.5  0.5  palindromic_swapped
h  C  G  G  C        0.2  0.2              flipped
i  G  T  A  C        0.5  0.5     flipped_and_swapped
j  G  T  C  A        0.5  0.5              flipped
k  A  G  A  G        0.2  0.2                 fine
l  T  G  G  T        0.2  0.8              swapped"""

    return pd.read_table(StringIO(c), sep="\s+")

@pytest.fixture
def expected_result_df1():

    c = """rsid E1 E2 O1 O2  OAF  EAF                 note                    SNPType
a  A  G  A  G  0.5  0.5                 fine                       fine
b  A  G  G  A  0.5  0.5              swapped                    swapped
c  A  G  T  A  0.5  0.5         incompatible               incompatible
d  A  T  A  T  0.2  0.8     flipped_and_swapped     inferrable_palindromic
e  A  T  A  T  0.2  0.2                 fine                       fine
f  A  T  A  T  0.5  0.5          palindromic  noninferrable_palindromic
g  A  T  T  A  0.5  0.5  palindromic_swapped        palindromic_swapped
h  C  G  G  C  0.2  0.2              flipped                    flipped
i  G  T  A  C  0.5  0.5     flipped_and_swapped           flipped_and_swapped
j  G  T  C  A  0.5  0.5              flipped                    flipped
k  A  G  A  G  0.2  0.2                 fine                       fine
l  T  G  G  T  0.2  0.8              swapped                    swapped"""

    return pd.read_table(StringIO(c), sep="\s+")

def test_bed_allele():

    with pytest.raises(KeyError) as exec_info:
        switch_alleles("B")

    assert "Non-valid allele in dataset: B" in str(exec_info)





def test_incompatible(df1):

    incom = incompatible(df1, 0.08)

    print(incom)

    assert incom.reset_index(drop=True).equals(pd.Series(["c"]))


def test_to_flip_and_swap(df1):

    sf = to_flip_and_swap(df1, 0.08)

    print(sf)
    assert sf.reset_index(drop=True).equals(pd.Series(["i"]))


def test_to_flip(df1):

    flip = to_flip(df1, 0.08)

    print(flip)
    assert flip.reset_index(drop=True).equals(pd.Series(["h", "j"]))


def test_to_swap(df1):

    swap = to_swap(df1, 0.08)

    print(swap)

    assert swap.reset_index(drop=True).equals(pd.Series(["b", "l"]))


def test_palindromic_to_swap(df1):

    palindromic_swap = palindromic_to_swap(df1, 0.08)

    print(palindromic_swap)
    assert palindromic_swap.reset_index(drop=True).equals(pd.Series(["g"]))


def test_is_fine1(df1):

    fine = fine_snps(df1, 0.08)

    assert fine.reset_index(drop=True).equals(pd.Series(["a", "e", "k"]))


def test_is_fine2(df2):

    fine = fine_snps(df2, 1.09)

    print(fine)

    assert fine.empty

def test_is_non_inferrable_palindromic_snp1(df1):

    non_inferrable_palindromic = non_inferrable_palindromic_snp(df1, 0.08)

    print(non_inferrable_palindromic)

    assert non_inferrable_palindromic.reset_index(drop=True).equals(pd.Series(["f"]))


def test_is_non_inferrable_palindromic_snp2(df2):

    non_inferrable_palindromic = non_inferrable_palindromic_snp(df2, 0.08)

    print(non_inferrable_palindromic)

    assert non_inferrable_palindromic.reset_index(drop=True).equals(pd.Series(["k", "l"]))


def test_inferrable_palindromic_snp1(df1):

    inferrable_palindromic = inferrable_palindromic_snp(df1, 0.08)

    print(inferrable_palindromic)

    assert inferrable_palindromic.reset_index(drop=True).equals(pd.Series(["d"]))

def test_is_inferrable_palindromic_snp2(df2):

    inferrable_palindromic = inferrable_palindromic_snp(df2, 0.08)

    print(inferrable_palindromic)

    assert inferrable_palindromic.reset_index(drop=True).equals(pd.Series(["g", "h", "i", "m"]))


def test_assign_table1(df1, expected_result_df1):

    result = assign_table(df1, 0.08)

    print(result)
    print(expected_result_df1)

    # print(result.dtypes)
    # print(expected_result_df1.dtypes)

    # print(result.index)
    # print(expected_result_df1.index)


    assert result.SNPType.equals(expected_result_df1.SNPType)



def test_assign_table2(df2, expected_result_df2):

    result = assign_table(df2, 0.08)

    print(result)
    print(expected_result_df2)

    # print(result.dtypes)
    # print(expected_result_df2.dtypes)

    # print(result.index)
    # print(expected_result_df2.index)


    assert result.equals(expected_result_df2)
