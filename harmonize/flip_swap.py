swap = -1
flip_and_swap = -1
palindromic_to_swap = -1

to_flip = 1


def flip_swap_remove(assigned_table, tol=0.08):

    snptype = assigned_table.SNPType.str

    swap_idx = snptype.contains("swap")
    flip_idx = snptype.contains("flip")
    remove_idx = snptype.contains("incompatible|noninferrable_palindromic")

    to_swap = list(assigned_table.ix[swap_idx].index.to_series())
    to_flip = list(assigned_table.ix[flip_idx].index.to_series())
    to_remove = list(assigned_table.ix[remove_idx].index.to_series())

    return to_swap, to_flip, to_remove


def swap_dosages(dosages, to_swap, to_remove):

    dosages = dosages.drop(to_remove, axis=1)
    dosages.loc[:, to_swap] = 2 - dosages.loc[:, to_swap]

    return dosages


def flip_exposure_alleles(exposure):
    pass
