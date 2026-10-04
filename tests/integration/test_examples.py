"""Smoke tests for the repository-only numerical examples."""

import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).parents[2]


@pytest.mark.parametrize("module,arguments", [
    ("examples.baker_formula.baker_formula_demo", ()),
    ("examples.bernatska.bernatska_trigonal_demo", ()),
    (
        "examples.bobenko_reyman_semenov_tian_shansky."
        "kowalewski_genus_three",
        ("--steps", "40", "--samples", "5", "--stop", "0.2"),
    ),
    (
        "examples.bobenko_reyman_semenov_tian_shansky.curve_data",
        (),
    ),
    (
        "examples.kovalevskaya_original.kovalevskaya_original_demo",
        ("--stage", "demo", "--tol", "1e-5"),
    ),
    (
        "examples.manakov.manakov_kleinian_demo",
        ("--dps", "20", "--steps", "40", "--samples", "5",
         "--start=-0.1", "--stop=0.1", "--phases", "0.13"),
    ),
    (
        "examples.neumann_moser.neumann_moser_kleinian_demo",
        ("--dps", "20", "--steps", "40", "--samples", "5",
         "--start=-0.1", "--stop=0.1", "--phases", "0.13"),
    ),
    (
        "examples.neumann_moser.neumann_moser_kleinian_demo",
        ("--dps", "25", "--start=0", "--stop=0.1",
         "--steps", "128", "--samples", "5", "--tol", "1e-8",
         "--initial-state", "-1", "-2", "0.1", "0.2", "1", "-13.75", "-7"),
    ),
    (
        "examples.neumann_moser.neumann_moser_genus_three_demo",
        ("--dps", "20", "--steps", "40", "--samples", "5",
         "--start=-0.1", "--stop=0.1", "--phases", "0.13"),
    ),
    (
        "examples.onishi.onishi_determinant_demo",
        ("--genus", "2", "--n", "2", "--kiepert-n", "2",
         "--coordinate", "1", "--dps", "20"),
    ),
    (
        "examples.onishi_trigonal_genus_3.lemma_5_1_demo",
        ("--dps", "25", "--curves", "purely-trigonal"),
    ),
    ("examples.rn_desitter_9d.rn_desitter_geodesic_demo", ("--help",)),
])
def test_example(module, arguments):
    """Each curated entry point runs without optional dependencies."""
    environment = os.environ.copy()
    environment["PYTHONWARNINGS"] = "error::DeprecationWarning"
    subprocess.run(
        [sys.executable, "-m", module, *arguments],
        cwd=ROOT,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
    )
