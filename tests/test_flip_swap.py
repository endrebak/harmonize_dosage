from io import StringIO

import pytest

import pandas as pd

from harmonize.flip_swap import flip_swap_remove


@pytest.fixture
def table():

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
    return pd.read_table(StringIO(c), sep="\s+", index_col=0)

@pytest.fixture()
def dosages():
    return pd.read_table(StringIO('''a b c d e f g h i j k l
0.5075130671930126 0.6646550304735508 0.45777424139957534 1.071761126886858 0.436854404298715 0.594036819968637 1.9349088247231299 1.8886705366489525 0.8162360277277216 1.0076002714788566 0.4520110894154601 0.2736366193037105
0.5004172840913197 0.21095987657153747 1.4870948694391923 1.8576468410105191 0.01851747237110657 1.7390990679354936 0.07719353296994136 0.4521553701206604 0.3487909379081455 1.782841803381271 0.371937965043585 1.4629459877760314
1.6616343437709837 1.8564140717822688 1.6097753873063108 0.4205175729986621 0.11166296561587719 0.27231866851836095 0.106642537805677 0.026029450444908786 1.5809031660710282 1.714374677015734 1.1555326731456808 1.1364352755316034
0.758296744995087 1.8732386417230384 1.0147184281761823 1.3730437486724474 0.02787779046727845 0.03312954383124178 0.0321372707362948 1.8856745302493647 1.2857408874034968 1.5140479339584874 1.516972296755377 0.2124103568336495
1.4005022406061318 0.39092144068589874 1.8659338378524475 1.6833101574797324 1.973712363287082 1.1483468391652196 0.17765884344199012 0.8203629222345827 1.0520400339857368 1.1891185263426203 0.40744333944562117 1.7210069030470356\n'''), sep=" ")


def test_flip_swap(dosages, table):

    to_swap, to_flip, to_remove = flip_swap_remove(table, 0.08)

    print(to_swap, to_flip, to_remove)

    assert to_swap == ['b', 'g', 'i', 'l'] and to_flip == ['h', 'i', 'j'] and to_remove == ['c', 'f']
