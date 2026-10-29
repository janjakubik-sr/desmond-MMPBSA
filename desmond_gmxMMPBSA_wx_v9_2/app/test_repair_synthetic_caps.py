#!/usr/bin/env python3
"""Focused regression tests for synthetic-cap nonbonded exclusions."""

from __future__ import annotations

import unittest

import numpy as np

from app.repair_synthetic_caps import _score_candidates, graph_neighborhood


class PairSpecificExclusionTests(unittest.TestCase):
    def test_one_five_pair_is_not_lost_to_cap_wide_neighborhood(self) -> None:
        # Model CH3-C-N-CA-HA.  CH3/HA is a 1-5 pair and must remain in the
        # nonbonded search even though HA is close to another atom in the cap.
        adjacency = [set() for _ in range(5)]
        for left, right in zip(range(4), range(1, 5)):
            adjacency[left].add(right)
            adjacency[right].add(left)

        old_cap_wide = graph_neighborhood((0, 1), adjacency, depth=3)
        methyl_pair_exclusions = graph_neighborhood((0,), adjacency, depth=2)

        self.assertIn(4, old_cap_wide)
        self.assertNotIn(4, methyl_pair_exclusions)

    def test_candidate_score_ignores_only_explicit_pair_exclusions(self) -> None:
        candidates = np.asarray(
            [
                [[[0.0, 0.0, 0.0]]],
                [[[0.20, 0.0, 0.0]]],
            ],
            dtype=float,
        ).reshape(2, 1, 3)
        environment_xyz = np.asarray([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0]])
        best, _, _, minimum = _score_candidates(
            candidates,
            np.asarray([0]),
            environment_xyz,
            np.asarray([1, 2]),
            np.asarray([0.17, 0.17, 0.17]),
            0.10,
            np.asarray([[True, False]]),
        )

        self.assertEqual(best, 1)
        self.assertAlmostEqual(float(minimum[0]), 0.05)
        self.assertAlmostEqual(float(minimum[1]), 0.15)


if __name__ == "__main__":
    unittest.main()
