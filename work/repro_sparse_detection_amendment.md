# Sparse detection during corrected-count inference

Recorded 2026-09-11 before the first corrected-count LIANA fit, while DecontX
models were still being completed. Raw/native full-network outputs remain
unchanged and retain their original script snapshot.

The installed LIANA 1.10.0 `_get_props` uses `X_mask.getnnz(axis=0)` divided by
the cell count. Sparse matrices can store explicit numerical zeros, which
`getnnz` includes. Fractional corrected estimates can underflow when converted
to the float32 normalized matrix used by the existing LIANA pipeline.

The current corrected-count pipeline therefore records positive-entry counts
before and after float32 normalization and the count of explicit zeros, then
calls `eliminate_zeros()` before LIANA. This enforces the intended >0 detection
rule without rounding estimated counts to integers. The independent custom
score gate uses actual >0 comparisons in float64 and separately reports >=1
estimated-count sensitivity. Code and count-source differences remain labeled.

The inspected upstream helper is copied to
`work/repro_literature/liana_1_10_0_common.py`, under its existing BSD license.
This numerical representation check does not change model fitting, identity
selection, expression-proportion threshold, resource, or animal labels.
