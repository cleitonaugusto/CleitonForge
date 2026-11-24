"""Tests for cleitonforge.qec. Run with: pytest cforge-py/tests"""

import pytest

stim = pytest.importorskip("stim")

from cleitonforge import qec  # noqa: E402


def surface_code(distance=3, rounds=3):
    return stim.Circuit.generated(
        "surface_code:rotated_memory_z",
        distance=distance,
        rounds=rounds,
        after_clifford_depolarization=0.001,
        before_measure_flip_probability=0.001,
    )


def test_witness_is_missed_by_the_dem_and_caught_by_the_flows():
    good, bad = qec.witness()
    verdict = qec.check_pass(good, bad)
    assert verdict.dem_equal
    assert not verdict.equivalent
    assert verdict.dem_blind
    assert verdict.witness is not None


def test_a_circuit_is_equivalent_to_itself():
    c = surface_code()
    verdict = qec.check_pass(c, c.copy())
    assert verdict.equivalent and verdict.dem_equal and not verdict.dem_blind


def test_removing_a_gate_nobody_reads_is_equivalent():
    # A basis change left before a reset, the case of Deltakit/deltakit#286.
    before = stim.Circuit("R 0\nH 0\nM 0\nH 0\nR 0\nM 0")
    after = stim.Circuit("R 0\nH 0\nM 0\nR 0\nM 0")
    assert qec.flows_equivalent(before, after)[0]


def test_a_trailing_gate_matters_only_if_the_final_state_does():
    before = stim.Circuit("R 0\nM 0\nH 0")
    after = stim.Circuit("R 0\nM 0")
    assert qec.flows_equivalent(before, after)[0]
    assert not qec.flows_equivalent(before, after, discard_final_state=False)[0]


def test_different_measurement_counts_differ():
    assert not qec.flows_equivalent(stim.Circuit("M 0"), stim.Circuit("M 0 0"))[0]


def test_deltakit_style_objects_are_converted():
    class Wrapped:
        def as_stim_circuit(self):
            return stim.Circuit("R 0\nM 0")

    assert qec.flows_equivalent(Wrapped(), stim.Circuit("R 0\nM 0"))[0]
    with pytest.raises(TypeError):
        qec.to_stim(object())


def test_the_dem_cannot_see_a_pauli_before_a_measurement():
    power = qec.oracle_power(qec.dem_differs, [surface_code()], trials=10)
    blind = power["X before a measurement"]
    assert blind.changed == 10
    assert blind.caught == 0
    assert power["drop a two-qubit gate"].rate == 1.0


def test_the_pauli_mutants_change_the_circuit_by_an_independent_measure():
    # The flows decide which mutants changed the circuit, so check them against
    # something that does not use flows: decode the mutant's noiseless
    # measurements with the original's detectors and observable. A fired
    # detector or a flipped observable means the mutant computes something else.
    import random

    c = surface_code().flattened()
    convert = c.without_noise().compile_m2d_converter()
    rng = random.Random(0)
    for _ in range(10):
        mutant = qec.mut_pauli_before_measurement(c, rng).without_noise()
        shots = mutant.compile_sampler().sample(64)
        det, obs = convert.convert(measurements=shots, separate_observables=True)
        assert det.any() or obs.any()
        assert qec.flows_differ(c, mutant)


def test_every_mutation_applies_to_a_surface_code():
    power = qec.oracle_power(qec.flows_differ, [surface_code()], trials=5)
    for name, p in power.items():
        assert p.changed + p.equivalent == 5, name


def test_mutations_do_not_touch_their_input():
    c = surface_code()
    text = str(c)
    qec.oracle_power(qec.dem_differs, [c], trials=3)
    assert str(c) == text
