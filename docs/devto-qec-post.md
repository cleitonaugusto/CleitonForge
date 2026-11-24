Quantum error correction compilers rewrite circuits all the time. They swap gates for native ones, drop gates nobody reads, merge layers. Every one of those rewrites (a "pass") has to leave the circuit doing the same thing. So how do you check that it did?

A common answer: compare the detector error model (DEM) before and after the pass. If the DEM is the same, the decoder sees the same thing, so the pass is fine. A preprint from June gives a fast, exact procedure for comparing DEMs and proposes it for exactly this job, verifying compiler passes ([arXiv:2606.14677](https://arxiv.org/abs/2606.14677)).

Here is a pass that gets through that check.

```
before the pass               after the pass
R 0                           R 0
X_ERROR(0.01) 0               X_ERROR(0.01) 0
M 0                           M 0
                              X 0
M 0                           M 0
DETECTOR rec[-1] rec[-2]      DETECTOR rec[-1] rec[-2]
```

The detector compares two measurements of the same qubit, so it should never fire. After the pass there is an X between them, and it fires on every single shot.

![Running the check: DEM identical, detector sampler 0.0, real detector 1.0, check_pass reports changed](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/ucwe2g4dytcxo0b76rcb.png)

The DEM is identical. Stim's detector sampler says 0.0. Neither one is buggy. Both are defined relative to the circuit's own noiseless run, and a deterministic X moves the measured value and that reference together, so there is nothing left for them to see. If you decode the new circuit's shots with the old circuit's detectors, you get 1.0, every time.

So comparing DEMs can't verify a pass on its own, however exact the comparison is. That isn't a mistake in anyone's math. It's a limit of what a DEM records.

## It's not just a toy

I took Stim's rotated surface code (memory in Z and in X, distance 3 and 5) and applied four kinds of change, 50 times each per circuit. Before counting anything, I confirmed every mutant really changed the circuit, without trusting my own tool: I decoded its noiseless shots with the original circuit's detectors. That agreed in 1200 out of 1200 cases.

![Bar chart: DEM comparison and detector sampling catch 0% of "X before a measurement" mutants, stabiliser flows catch 100%](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/0npjs2k2zb8moqku5yw6.png)

Three of the four changes are caught by everything. The fourth, one X right before a measurement, is caught by neither the DEM comparison nor detector sampling. Zero out of 200.

What does catch it is checking the circuit's **stabiliser flows** in both directions: every flow of the old circuit has to be a flow of the new one, and the other way round.

## Then my own check fell into the opposite hole

I was ready to write "use flows, not DEMs". Then the repetition code broke that sentence.

In a repetition code memory every qubit starts in |0>. A CX between qubits in |0> does nothing to the noiseless state. So if a pass drops one of the parity-check CXs, the circuit still computes exactly the same thing, and the flows say "equivalent". They're right. But under noise, an error that CX used to carry to an ancilla now reaches fewer detectors. The DEM catches it 100 times out of 100. The flows catch it zero times.

![Table: X before a measurement is blind to the DEM and caught by flows; dropping a parity-check CX is caught by the DEM and blind to flows; check_pass catches both](https://dev-to-uploads.s3.us-east-2.amazonaws.com/uploads/articles/qehp88u74cc9s9wguxaf.png)

Each check is blind exactly where the other one sees. A pass is validated by both or it isn't validated.

## Try it

I put both checks in one function:

```bash
pip install "cleitonforge[qec]"
```

```python
from cleitonforge import qec

v = qec.check_pass(before, after)   # stim.Circuit, or a deltakit_circuit.Circuit
print(v)
v.equivalent       # flows match both ways
v.dem_equal        # what a DEM comparison would say
v.dem_blind        # DEM equal, circuit changed
v.faults_changed   # circuit equal, DEM changed: errors are detected differently
```

And if you already have a test for your passes, you can measure what it can't see:

```python
qec.oracle_power(my_check, circuits)   # my_check(before, after) -> bool
```

It applies those same mutations and tells you, per kind, how many real changes your check let through.

## I ran it on a real compiler

[Deltakit](https://github.com/Deltakit), from Riverlane, has two compilers. I ran every circuit pass I could run in both of them: native gate compilation, dead gate removal, layer merging, the MLIR passes between Stim import and export. That came to 480 before-and-after pairs, on surface, repetition and planar codes, with and without noise. None of them changed a circuit.

The only bug the run found was in my own tool. The very first run, with no pass at all, said the DEM had changed everywhere, because the import and export rewrite detector coordinates and I was comparing DEMs as text. A decoder doesn't care about coordinates. Fixed in 0.6.1, and the comparison now works on error mechanisms.

## Limits, so nobody has to dig for them

Distances 3 and 5, memory experiments only. No lattice surgery, nothing adaptive. Four kinds of mutation. "No problem found in Deltakit" means no problem on these inputs, not that the passes are proven correct.

## Links

- The note with all the numbers: [doi:10.5281/zenodo.24517526](https://doi.org/10.5281/zenodo.24517526)
- Why detector-based checks are blind here, in more depth: [doi:10.5281/zenodo.22802524](https://doi.org/10.5281/zenodo.22802524)
- Code and the script behind every number: [github.com/cleitonaugusto/CleitonForge](https://github.com/cleitonaugusto/CleitonForge)
- Package: [pypi.org/project/cleitonforge](https://pypi.org/project/cleitonforge/)

If you maintain a QEC compiler and want it run against your passes, open an issue on the repo. I'd like to see what it finds.
