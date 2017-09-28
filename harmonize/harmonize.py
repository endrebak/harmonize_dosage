from itertools import product

import pandas as pd

from numpy import isclose



switch_alleles = lambda a: {"A": "T",
                            "T": "A",
                            "C": "G",
                            "G": "C",}[a] # "N": "N"


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


def is_inferrable_palindromic_snp(m, tolerance):

    x = _inferrable_palindrome(m, tolerance) & ~_fine_snps(m, tolerance)

    return m.loc[x].rsid


def _to_swap(m, tolerance):

    return (m.E1 == m.O2) & (m.E2 == m.O1) & ~_palindromic_to_swap(m, tolerance) & ~_to_flip(m, tolerance)


def to_swap(m, tolerance):

    return m.loc[_to_swap(m, tolerance)].rsid


def palindromic_to_swap(m, tolerance):

    return m.loc[_palindromic_to_swap(m, tolerance)].rsid


def _palindromic_to_swap(m, tolerance):

    return (m.E1 == m.O2) & (m.E2 == m.O1) & (m.E1.apply(switch_alleles) == m.O1) & (m.E2.apply(switch_alleles) == m.O2) & ~_exposure_and_outcome_allele_frequencies_differ_from_fifty_percent(m, tolerance)


def to_flip(m, tolerance):

    return m.loc[_to_flip(m, tolerance)].rsid

def _to_flip(m, tolerance):

    return _allele_frequencies_similar_within_tolerance(m, tolerance) & (m.E1.apply(switch_alleles) == m.O1) & (m.E2.apply(switch_alleles) == m.O2) & ~_palindromic_to_swap(m, tolerance)



def _to_swap_and_flip(m, tolerance):
    # to swap(m.E1 == m.O2) & (m.E2 == m.O1) & ~_palindromic_to_swap(m, tolerance)

    return (m.E1.apply(switch_alleles) == m.O2) & (m.E2.apply(switch_alleles) == m.O1) & ~_fine_snps(m, tolerance) & ~_palindrome(m)


def to_swap_and_flip(m, tolerance):

    return m.loc[_to_swap_and_flip(m, tolerance)].rsid

def _incompatible(m, tolerance):

    return ~_to_flip(m, tolerance) & ~_palindromic_to_swap(m, tolerance) & ~_to_swap(m, tolerance) & ~_inferrable_palindrome(m, tolerance) & ~_non_inferrable_palindrome(m, tolerance) & ~_fine_snps(m, tolerance) & ~_to_swap_and_flip(m, tolerance)


def incompatible(m, tolerance):

    return m.loc[_incompatible(m, tolerance)].rsid


def assign_table(m, tolerance):

    d = {
         "fine": fine_snps(m, tolerance),
         "inferrable_palindromic": is_inferrable_palindromic_snp(m, tolerance),
         "noninferrable_palindromic": non_inferrable_palindromic_snp(m, tolerance),
         "to_swap": to_swap(m, tolerance),
         "to_flip": to_flip(m, tolerance),
         "to_swap_and_flip": to_swap_and_flip(m, tolerance),
         "palindromic_to_swap": palindromic_to_swap(m, tolerance),
         "incompatible": incompatible(m, tolerance)
    }

    for (li, i), (lj, j) in product(d.items(), repeat=2):
        i, j = set(i), set(j)
        if not li == lj:
            assert not i.intersection(j), "{}: {}\n{} {}\n{}".format(li, str(i), lj, str(j), set(i).intersection(j))

    m.insert(len(m.columns), "Harmonize", "")
    for l, v in d.items():
        print(l)
        print(v.index)
        m.loc[v.index, "Harmonize"] = l

    print(m)
