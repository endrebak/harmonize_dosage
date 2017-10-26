from io import StringIO

import pytest

import pandas as pd
from numpy import allclose

from harmonize.flip_swap import (flip_swap_remove, flip_swap_remove_exposure_alleles, swap_dosages)


@pytest.fixture
def exposure_table():

    c = """rs ALT REF
a  A  G
b  A  G
c  A  G
d  A  T
e  A  T
f  A  T
g  A  T
h  C  G
i  G  T
j  G  T
k  A  G
l  T  G"""

    return pd.read_table(StringIO(c), sep="\s+", index_col=0)


@pytest.fixture
def expected_result_exposure_table():

    c = """rs ALT REF
a  A  G
b  G  A
d  A  T
e  A  T
g  T  A
h  G  C
i  A  C
j  C  A
k  A  G
l  G  T"""

    return pd.read_table(StringIO(c), sep="\s+", index_col=0)

@pytest.fixture
def table():

    c = """rsid D1 D2 G1 G2  GAF  DAF                 note                    SNPType
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
    return pd.read_table(StringIO(c), sep="\s+", index_col=0)


@pytest.fixture()
def dosages():
    return pd.read_table(StringIO('''a_dose b_dose c_dose d_dose e_dose f_dose g_dose h_dose i_dose j_dose k_dose l_dose
0.5075130671930126 0.6646550304735508 0.45777424139957534 1.071761126886858 0.436854404298715 0.594036819968637 1.9349088247231299 1.8886705366489525 0.8162360277277216 1.0076002714788566 0.4520110894154601 0.2736366193037105
0.5004172840913197 0.21095987657153747 1.4870948694391923 1.8576468410105191 0.01851747237110657 1.7390990679354936 0.07719353296994136 0.4521553701206604 0.3487909379081455 1.782841803381271 0.371937965043585 1.4629459877760314
1.6616343437709837 1.8564140717822688 1.6097753873063108 0.4205175729986621 0.11166296561587719 0.27231866851836095 0.106642537805677 0.026029450444908786 1.5809031660710282 1.714374677015734 1.1555326731456808 1.1364352755316034
0.758296744995087 1.8732386417230384 1.0147184281761823 1.3730437486724474 0.02787779046727845 0.03312954383124178 0.0321372707362948 1.8856745302493647 1.2857408874034968 1.5140479339584874 1.516972296755377 0.2124103568336495
1.4005022406061318 0.39092144068589874 1.8659338378524475 1.6833101574797324 1.973712363287082 1.1483468391652196 0.17765884344199012 0.8203629222345827 1.0520400339857368 1.1891185263426203 0.40744333944562117 1.7210069030470356\n'''), sep=" ")

@pytest.fixture()
def expected_dosages():
    c = """a_dose b_dose d_dose e_dose g_dose h_dose i_dose j_dose k_dose l_dose
0.5075130671930126 1.3353449695264492 0.9282388731131421 0.43685440429871497 0.06509117527687036 1.8886705366489525 1.1837639722722784 1.0076002714788566 0.4520110894154601 1.7263633806962895
0.5004172840913197 1.7890401234284625 0.14235315898948087 0.01851747237110657 1.9228064670300586 0.4521553701206604 1.6512090620918545 1.782841803381271 0.371937965043585 0.5370540122239686
1.6616343437709835 0.1435859282177312 1.579482427001338 0.1116629656158772 1.893357462194323 0.026029450444908786 0.4190968339289718 1.714374677015734 1.1555326731456808 0.8635647244683966
0.758296744995087 0.12676135827696156 0.6269562513275526 0.02787779046727845 1.9678627292637052 1.8856745302493647 0.7142591125965032 1.5140479339584874 1.5169722967553771 1.7875896431663505
1.4005022406061318 1.6090785593141013 0.31668984252026755 1.9737123632870817 1.8223411565580099 0.8203629222345827 0.9479599660142632 1.1891185263426205 0.4074433394456212 0.27899309695296437"""
    return pd.read_table(StringIO(c), sep="\s+")


def test_flip_swap(dosages, table):

    to_swap, to_flip, to_remove = flip_swap_remove(table, 0.08)

    print(to_swap, to_flip, to_remove)

    assert to_swap == ['b', 'd', 'g', 'i', 'l'] and to_flip == ['d', 'h', 'i', 'j'] and to_remove == ['c', 'f']


def test_correct_dosages(dosages, expected_dosages):

    swapped = swap_dosages(dosages, ['b', 'd', 'g', 'i', 'l'], ['c', 'f'])

    print(swapped)
    print(expected_dosages)

    print(swapped.dtypes)
    print(expected_dosages.dtypes)

    print(swapped.index)
    print(expected_dosages.index)

    assert allclose(swapped, expected_dosages)


def test_flip_swap_exposure_alleles(exposure_table, expected_result_exposure_table):

    actual_result = flip_swap_remove_exposure_alleles(exposure_table, ['d', 'h', 'i', 'j'], ['b', 'd', 'g', 'i', 'l'], ['c', 'f'])

    print(actual_result, "actual")
    print(expected_result_exposure_table, "expected")

    print(actual_result.dtypes)
    print(expected_result_exposure_table.dtypes)

    print(actual_result.index)
    print(expected_result_exposure_table.index)

    assert actual_result.equals(expected_result_exposure_table)
