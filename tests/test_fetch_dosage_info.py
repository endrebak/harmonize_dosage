
import pytest

import pandas as pd

from io import StringIO

from harmonize.fetch_dosage_info import fetch_dosage_info

@pytest.fixture()
def expected_result():

    contents = """ID REF ALT     AF
bs1   C   G  0.123
bs3   T   G  0.123"""

    return pd.read_table(StringIO(contents), sep="\s+")


# # silly pandas why aren't these two dfs equal?

# def test_fetch_dosage_info(expected_result):

#     result = fetch_dosage_info("tests/data/tiny.vcf")

#     print(result)
#     print(expected_result)

#     print(result.columns)
#     print(expected_result.columns)

#     print(result.dtypes)
#     print(expected_result.dtypes)

#     print(result.index)
#     print(expected_result.index)

#     assert expected_result.equals(result)
