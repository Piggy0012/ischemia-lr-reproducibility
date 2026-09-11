# Derived executable subset of scverse/liana merge d4211373692e7b9c10210488ccb1efe06452b097.
# Function bodies are AST-identical to the archived upstream source. No new LIANA package is loaded.
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import beta, rankdata
def _logg(*args, **kwargs):
    raise RuntimeError('Logging branch was not expected for this three-column aggregation audit')
def _rank_aggregate(lr_res: pd.DataFrame, specs: dict[str, tuple[str, bool | None]], aggregate_method: Literal['rra', 'mean'], verbose: bool=False) -> NDArray[np.floating]:
    """
    Aggregate method ranks

    Parameters
    ----------
    lr_res
        joined results from all methods
    specs
        specs dictionary where method_name:(score_name, score_desc)
    aggregate_method
        method by which to aggregate the ranks
    verbose
        Verbosity flag

    Returns
    -------
    An array of values /w length of lr_res.shape[0]
    """
    if aggregate_method not in ('rra', 'mean'):
        raise ValueError(f"`aggregate_method` must be 'rra' or 'mean', got {aggregate_method!r}.")
    specs = {method: spec for method, spec in specs.items() if spec[0] in lr_res.columns}
    columns: dict[str, bool | None] = {}
    for col, asc in specs.values():
        if columns.setdefault(col, asc) != asc:
            raise ValueError(f'Column `{col}` is ranked in opposite directions by different methods.')
    if len(columns) < 2:
        _logg(f'Aggregating over {len(columns)} score(s) only: {sorted(columns)}.', level='warn', verbose=verbose)
    rmat: NDArray[np.floating] = np.column_stack([rankdata(lr_res[col] if asc else -lr_res[col], method='average') for col, asc in columns.items()])
    if aggregate_method == 'rra':
        return _robust_rank_aggregate(rmat)
    return np.asarray(np.mean(rmat, axis=1) / rmat.shape[0], dtype=np.float64)

def _corr_beta_pvals(p: NDArray[np.floating], k: int) -> NDArray[np.floating]:
    """
    Correct beta p-values

    Parameters
    ----------
    p
        (min) p-value
    k
        total number of rows

    Returns
    -------
    An array with corrected p-values
    """
    p = np.clip(p * k, a_min=0, a_max=1)
    return p

def _rho_scores(rmat: NDArray[np.floating], dist_a: NDArray[np.integer], dist_b: NDArray[np.integer]) -> NDArray[np.floating]:
    """
    Calculate Beta Distribution Rho Scores

    Parameters
    ----------
    rmat
        a matrix where rows are the ranks/n for each interaction, while
        columns correspond to each method
    dist_a
        non-negative shape param a
    dist_b
        non-negative shape param b

    Returns
    -------
    A vector of pvals as implemented in the RRA method
    """
    rmat = np.sort(rmat, axis=1)
    p = beta.cdf(rmat, dist_a, dist_b)
    p = np.min(p, axis=1)
    rho = _corr_beta_pvals(p, k=rmat.shape[1])
    return rho

def _robust_rank_aggregate(rmat: NDArray[np.floating]) -> NDArray[np.floating]:
    """
    Calculate Robust Rank Aggregate as in Kolde et al., 2012

    Parameters
    ----------
    rmat
        Matrix with interaction ranks (rows) for each method (columns)

    Returns
    -------
    An array with p-values for each row
    """
    rmat = rmat / np.max(rmat, axis=0)
    dist_a = np.repeat([np.arange(rmat.shape[1])], rmat.shape[0], axis=0) + 1
    dist_b = rmat.shape[1] - dist_a + 1
    return _rho_scores(rmat, dist_a, dist_b)
