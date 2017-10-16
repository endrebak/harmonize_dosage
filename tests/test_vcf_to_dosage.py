import pytest
import pandas as pd

from numpy import float64

from io import StringIO

from harmonize.vcf_to_dosage import vcf_to_dosage


@pytest.fixture
def vcf():

    path = "tests/data/tiny.vcf"

    return path


@pytest.fixture
def expected_result():
    contents = """ID    bs1_dose  bs3_dose
NA1   1.3    2
NA2  0.01    0
NA3     2    1"""

    return pd.read_table(StringIO(contents), sep="\s+", header=0, index_col=0,
                         dtype={"bs1_dose": float64, "bs3_dose": float64})


def test_vcf_to_dosage(vcf, expected_result):

    result = vcf_to_dosage(vcf)

    print(result)

    print(expected_result)

    assert result.equals(expected_result)
