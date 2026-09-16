# Detector-Based Oracles Are Empty Where QEC Compilers Must Be Validated

Source and built PDF of the paper describing the measurements in
[`../research/qec_oracle_power.py`](../research/qec_oracle_power.py).

Everything it claims is reproducible from this repository:

| claim in the paper | where it comes from |
|---|---|
| 335 mutants, detectors and DEM detect none | `research/qec_oracle_power.py` |
| two negative controls, 200/200 equivalent | same file, `MUTATIONS` |
| Stim transformations clean, with declared power | `tools/fuzz_stim_transforms.py` |

Build with `pdflatex paper_oracle_power_qec.tex` twice.
