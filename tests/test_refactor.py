"""Current state snapshots must remain independent of mutable live arrays."""

from evo_sim.models.plant_insect import simulate
from evo_sim.models.plant_insect.helpers.agents import CommunityState
from tests.test_acceptance import tiny


def test_snapshot_independence():
    result = simulate(tiny())
    state = CommunityState.from_arrays(result["snapshots"][0])
    saved = state.snapshot(0)
    state.plants.defense[0, 0] = 0.987
    assert saved["D"][0, 0] != 0.987
    assert len(state.plants) == len(saved["plant_ids"])
    assert len(state.insects) == len(saved["insect_ids"])


if __name__ == "__main__":
    test_snapshot_independence()
    print("PASS snapshot independence")
