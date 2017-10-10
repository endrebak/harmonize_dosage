import pytest

import pandas as pd
from io import StringIO

from harmonize.create_exposure_outcome_table import create_exposure_outcome_table


@pytest.fixture
def exposure_table():
    c = """REF ALT rsid AF
A G RS1 0.63
C T RS2 0.14"""

    return pd.read_table(StringIO(c), sep="\s+", header=0, index_col=2)


@pytest.fixture
def outcome_table():
    c = """rsid A1 A2 FreqA1
RS1 A g 0.53
RS2 C t 0.42"""

    return pd.read_table(StringIO(c), sep="\s+").set_index("rsid")


@pytest.fixture
def expected_result():
    c = """rsid E1 E2 O1 O2 OAF EAF
RS1 A G A G 0.53 0.63
RS2 C T C T 0.42 0.14"""

    return pd.read_table(StringIO(c), sep="\s+")


def test_create_exposure_outcome_table(exposure_table, outcome_table, expected_result):

    result = create_exposure_outcome_table(exposure_table, outcome_table)

    print(result)
    print(expected_result)

    print(result.index)
    print(expected_result.index)

    print(result.dtypes)
    print(expected_result.dtypes)

    assert result.equals(expected_result)
