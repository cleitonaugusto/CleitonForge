"""The smallest error-correction circuit where the usual oracles are blind.

Four instructions, one qubit, one detector. Two circuits that differ by a
single inserted X:

    R 0                     R 0
    M 0                     M 0
                            X 0        <- the fault
    M 0                     M 0
    DETECTOR rec[-1] rec[-2]    DETECTOR rec[-1] rec[-2]

The detector is the parity of two measurements of the same qubit, so it is
deterministic. Inserting the X flips the second measurement, which flips that
parity from 0 to 1.

And yet:

  detector error models          identical
  detection-event rate, each     0.0000 and 0.0000
  circuit sampled on its own

  detection-event rate, the      0.0000 and 1.0000
  mutant judged by the
  original's reference

The reason both self-referential numbers are zero is that stim builds the
detector reference sample from whichever circuit it is given. The fault moves
the outcome and the reference together, so nothing is left over to see. The
DEM is identical for the same reason: a deterministic gate contributes no error
mechanism, it relocates the frame the mechanisms are expressed in.

This matters because "identical DEM, identical detector samples" is a common
way to argue that a compiler pass is safe. It is a sound argument only for
faults that are not deterministic Paulis, and nothing about the method says so.

    .venv-qec/bin/python qec_minimal_witness.py
"""
import stim

SHOTS = 10000

GOOD = stim.Circuit("""
    R 0
    M 0
    M 0
    DETECTOR rec[-1] rec[-2]
""")

BAD = stim.Circuit("""
    R 0
    M 0
    X 0
    M 0
    DETECTOR rec[-1] rec[-2]
""")


def main():
    print(__doc__.split("\n\n")[0])
    print()

    print("self-referential (each circuit sampled with its own sampler):")
    for name, circuit in (("correct", GOOD), ("X inserted", BAD)):
        rate = circuit.compile_detector_sampler(seed=1).sample(SHOTS).mean()
        print(f"  {name:<12} detection-event rate {rate:.4f}")

    print("\ncross-referential (mutant judged by the original's definitions):")
    converter = GOOD.compile_m2d_converter()
    for name, circuit in (("correct", GOOD), ("X inserted", BAD)):
        records = circuit.compile_sampler(seed=1).sample(SHOTS)
        rate = converter.convert(measurements=records,
                                 separate_observables=False).mean()
        print(f"  {name:<12} detection-event rate {rate:.4f}")

    same_dem = (str(GOOD.detector_error_model()) ==
                str(BAD.detector_error_model()))
    same_records = (GOOD.compile_sampler(seed=1).sample(SHOTS) ==
                    BAD.compile_sampler(seed=1).sample(SHOTS)).all()

    print(f"\ndetector error models identical: {same_dem}")
    print(f"measurement records identical:   {same_records}")

    assert same_dem, "the DEMs must agree for this witness to mean anything"
    assert not same_records, "the circuits must actually differ"
    print("\nThe circuits differ, the DEM cannot tell, and neither can a "
          "detector\nsampler that takes its reference from the circuit under "
          "test.")


if __name__ == "__main__":
    main()
