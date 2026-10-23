import numpy as np
import pandas as pd

from evo_sim.models.plant_insect.analysis.plot_diagnostics import cell_abundance, centered_xy


def test_plot_coordinates_and_extinction():
    p1, p2, h = cell_abundance(
        dict(
            plant_cells=np.array([[2, 7]]),
            plant_species=np.array([0]),
            insect_cells=np.array([[3, 8]]),
        )
    )
    assert p1[2, 7] == 1 and h[3, 8] == 1 and p2.sum() == 0
    data = pd.DataFrame(
        [
            dict(species="H1", case="x", replicate=1, generation=g, trait=t, mean=v)
            for g, v in [(0, 0.2), (1, 0.3), (2, np.nan)]
            for t in ["I_1", "I_2"]
        ]
    )
    g, x, y = centered_xy(data, "H1", "x", 1, ("I_1", "I_2"))
    assert list(g) == [0, 1]
    np.testing.assert_allclose(x, [0, 0.1])


def test_summary_batch_has_no_figure_side_effect():
    import json
    from pathlib import Path
    from tempfile import TemporaryDirectory

    from evo_sim.models.plant_insect.analysis.common import Batch, batch_parser

    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "data").mkdir()
        (root / "audits").mkdir()
        (root / "audits/integrity_audit.json").write_text(json.dumps({"runs": 1}))
        (root / "data/run_status.csv").write_text("case,replicate,status\nx,1,completed\n")
        batch = Batch.open(root)
        assert batch.summary_dir.is_dir()
        assert not (root / "figures").exists()
        assert not hasattr(batch_parser("summary").parse_args(["--input", tmp]), "out")
