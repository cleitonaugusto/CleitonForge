# Detector-Based Oracles Are Empty Where QEC Compilers Must Be Validated

**DOI: [10.5281/zenodo.22802524](https://doi.org/10.5281/zenodo.22802524)**

Source and built PDF of the paper describing the measurements in
[`../research/qec_oracle_power.py`](../research/qec_oracle_power.py).

Everything it claims is reproducible from this repository:

| claim in the paper | where it comes from |
|---|---|
| 335 mutants, detectors and DEM detect none | `research/qec_oracle_power.py` |
| two negative controls, 200/200 equivalent | same file, `MUTATIONS` |
| Stim transformations clean, with declared power | `tools/fuzz_stim_transforms.py` |

Build with `pdflatex paper_oracle_power_qec.tex` twice.
