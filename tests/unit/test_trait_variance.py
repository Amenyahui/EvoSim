"""Count-weighted descriptive variance; type differences are not additive genetic variance."""

import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.summaries.qtl import partition_trait_variance
from tests.support import raises


def test_trait_variance_partition_weights_types_and_preserves_empty_states():
    traits = pd.DataFrame(
        dict(
            generation=[0, 0, 0, 1],
            kind=["plant"] * 4,
            trait=["x"] * 4,
            species=["P1", "P2", "P3", "P1"],
            N=[1, 3, 0, 0],
            mean=[0.0, 1.0, np.nan, np.nan],
            variance=[0.04, 0.01, np.nan, np.nan],
        )
    )
    pooled, absent = partition_trait_variance(traits).iloc
    assert pooled.N == 4 and pooled.type_count == 2 and pooled["mean"] == 0.75
    np.testing.assert_allclose(
        [pooled.within_type_variance, pooled.between_type_variance, pooled.total_variance],
        [0.0175, 0.1875, 0.205],
    )
    assert absent.N == 0 and absent.type_count == 0 and np.isnan(absent.total_variance)
    assert partition_trait_variance(traits.iloc[:0]).empty
    raises(lambda: partition_trait_variance(pd.concat([traits, traits.iloc[:1]])), ValueError)
    for column, value in (("N", -1), ("N", 0.5), ("mean", np.nan), ("variance", -1)):
        invalid = traits.astype({"N": float})
        invalid.loc[0, column] = value
        raises(lambda: partition_trait_variance(invalid), ValueError)
