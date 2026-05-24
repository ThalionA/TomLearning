"""Spatial CCA pipeline for inter-areal communication across learning -- Tom edition.

A Tom-dataset port of the striatum CCA pipeline (see
``StriatumACC/Striatum project/cca/``). The same scientific question (does
inter-areal communication strengthen / reorient / change direction across
learning?) and the same analyses (lagged CCA + IFI, held-out-CC permutation
nulls, membership, principal angles, partial CCA, cross-pair subspace
similarity), applied to the hippocampal/cortical cohort with the eight pair
list Tom hard-codes in ``HC_V1_Code/HC_V1_temporal.m``.
"""

from . import (
    analysis,
    config,
    core,
    crosspair,
    dataio,
    lagged,
    membership,
    partial,
    pipeline,
    stage3,
    subspace,
    surrogate,
    sweep,
)

__all__ = [
    "analysis",
    "config",
    "core",
    "crosspair",
    "dataio",
    "lagged",
    "membership",
    "partial",
    "pipeline",
    "stage3",
    "subspace",
    "surrogate",
    "sweep",
]
