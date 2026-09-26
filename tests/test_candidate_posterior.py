"""Tests for the candidate-posterior endpoint, with a brute-force cross-check."""

from __future__ import annotations

import itertools
import random
from fractions import Fraction

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def parse_fraction(text: str) -> Fraction:
    n, d = text.split("/")
    return Fraction(int(n), int(d))


def fmt(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def brute_force_peaks(atom_groups: list) -> dict[int, Fraction]:
    """Enumerate every isotope-choice path explicitly (no merging shortcuts)."""
    dist: dict[int, Fraction] = {}
    per_group_outcomes = []
    for group in atom_groups:
        outcomes = []
        for combo in itertools.product(group["isotopes"], repeat=group["count"]):
            shift = sum(iso["shift_mda"] for iso in combo)
            prob = Fraction(1)
            for iso in combo:
                prob *= parse_fraction(iso["abundance"])
            outcomes.append((shift, prob))
        per_group_outcomes.append(outcomes)
    for path in itertools.product(*per_group_outcomes):
        shift = sum(item[0] for item in path)
        prob = Fraction(1)
        for item in path:
            prob *= item[1]
        dist[shift] = dist.get(shift, Fraction(0)) + prob
    return {s: p for s, p in dist.items() if p != 0}


def brute_force_posterior(payload: dict):
    """Enumerate every ion sampling sequence for every (candidate, offset).

    For each candidate and each shared offset delta, sums the probability of
    all ion draw sequences whose measured shifts (drawn shift + delta) tally
    exactly the listed observation counts, with every unlisted measured shift
    counted in "other". Returns (joint, candidate_marginal, offset_marginal)
    as nested Fraction lists, or None when the observations are impossible.
    """
    envelopes = [
        brute_force_peaks(candidate["atom_groups"])
        for candidate in payload["candidates"]
    ]
    priors = [parse_fraction(c["prior"]) for c in payload["candidates"]]
    radius = payload["calibration_radius"]
    ion_count = payload["ion_count"]
    obs_counts = {o["shift_mda"]: o["count"] for o in payload["observations"]}
    other_count = ion_count - sum(obs_counts.values())
    offsets = list(range(-radius, radius + 1))

    unnorm = []
    for envelope, prior in zip(envelopes, priors):
        shifts = sorted(envelope)
        row = []
        for delta in offsets:
            likelihood = Fraction(0)
            for draws in itertools.product(shifts, repeat=ion_count):
                tallies = dict.fromkeys(obs_counts, 0)
                other = 0
                prob = Fraction(1)
                for shift in draws:
                    prob *= envelope[shift]
                    measured = shift + delta
                    if measured in tallies:
                        tallies[measured] += 1
                    else:
                        other += 1
                if other == other_count and all(
                    tallies[s] == n for s, n in obs_counts.items()
                ):
                    likelihood += prob
            row.append(prior * likelihood / len(offsets))
        unnorm.append(row)

    total = sum((w for row in unnorm for w in row), Fraction(0))
    if total == 0:
        return None
    joint = [[w / total for w in row] for row in unnorm]
    candidate_marginal = [sum(row, Fraction(0)) for row in joint]
    offset_marginal = [
        sum((row[k] for row in joint), Fraction(0)) for k in range(len(offsets))
    ]
    return joint, candidate_marginal, offset_marginal


def symmetric_payload():
    # Candidate 0 envelope: 0 -> 1/2, 1 -> 1/2.
    # Candidate 1 envelope: 0 -> 1/2, 2 -> 1/2.
    return {
        "candidates": [
            {
                "atom_groups": [
                    {
                        "count": 1,
                        "isotopes": [
                            {"shift_mda": 0, "abundance": "1/2"},
                            {"shift_mda": 1, "abundance": "1/2"},
                        ],
                    }
                ],
                "prior": "1/2",
            },
            {
                "atom_groups": [
                    {
                        "count": 1,
                        "isotopes": [
                            {"shift_mda": 0, "abundance": "1/2"},
                            {"shift_mda": 2, "abundance": "1/2"},
                        ],
                    }
                ],
                "prior": "1/2",
            },
        ],
        "ion_count": 1,
        "observations": [{"shift_mda": 1, "count": 1}],
        "calibration_radius": 1,
    }


def assert_reduced(text: str) -> None:
    n, d = text.split("/")
    assert Fraction(int(n), int(d)).denominator == int(d)


def check_invariants(body: dict, n_candidates: int, radius: int) -> None:
    width = 2 * radius + 1
    joint = body["joint_posteriors"]
    assert len(joint) == n_candidates
    assert all(len(row) == width for row in joint)
    assert len(body["candidate_posteriors"]) == n_candidates
    assert len(body["offset_posteriors"]) == width
    for row in joint:
        for text in row:
            assert_reduced(text)
    for text in [*body["candidate_posteriors"], *body["offset_posteriors"]]:
        assert_reduced(text)
    joint_values = [[parse_fraction(t) for t in row] for row in joint]
    # The joint is a full distribution and both marginals are consistent.
    assert sum((v for row in joint_values for v in row), Fraction(0)) == 1
    for row, text in zip(joint_values, body["candidate_posteriors"]):
        assert sum(row, Fraction(0)) == parse_fraction(text)
    for k, text in enumerate(body["offset_posteriors"]):
        column = sum((row[k] for row in joint_values), Fraction(0))
        assert column == parse_fraction(text)


def test_symmetric_hand_computed():
    # Measured shift 1 is reachable as (candidate 0, delta 0), (candidate 0,
    # delta 1), (candidate 1, delta -1) and (candidate 1, delta 1), each with
    # likelihood 1/2 and prior weight 1/2: four equal posterior cells of 1/4,
    # so delta 1 (explained by both candidates) carries offset mass 1/2.
    resp = client.post("/v1/candidate-posterior", json=symmetric_payload())
    assert resp.status_code == 200
    assert resp.json() == {
        "possible": True,
        "joint_posteriors": [["0/1", "1/4", "1/4"], ["1/4", "0/1", "1/4"]],
        "candidate_posteriors": ["1/2", "1/2"],
        "offset_posteriors": ["1/4", "1/4", "1/2"],
    }


def test_zero_radius_identifies_candidate():
    # Without calibration freedom only candidate 0 can produce shift 1.
    payload = symmetric_payload()
    payload["calibration_radius"] = 0
    resp = client.post("/v1/candidate-posterior", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "possible": True,
        "joint_posteriors": [["1/1"], ["0/1"]],
        "candidate_posteriors": ["1/1", "0/1"],
        "offset_posteriors": ["1/1"],
    }


def test_other_bucket_absorbs_unlisted_ions():
    # Identical candidates (envelope 0 -> 1/2, 1 -> 1/2), two ions, one
    # observed at shift 1: the second ion falls into "other". Likelihood per
    # candidate: sequences (0, 1) and (1, 0), total 1/2; the posterior keeps
    # the 1/2-1/2 priors.
    payload = symmetric_payload()
    payload["candidates"][1] = dict(payload["candidates"][0])
    payload["ion_count"] = 2
    payload["calibration_radius"] = 0
    resp = client.post("/v1/candidate-posterior", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "possible": True,
        "joint_posteriors": [["1/2"], ["1/2"]],
        "candidate_posteriors": ["1/2", "1/2"],
        "offset_posteriors": ["1/1"],
    }


def test_multinomial_weighting_hand_computed():
    # Candidate 0 envelope: 0 -> 9/16, 1 -> 3/8, 2 -> 1/16 (two atoms);
    # candidate 1 envelope: 0 -> 1/2, 2 -> 1/2 (one atom). Two ions both
    # measured at shift 0: likelihoods (9/16)^2 = 81/256 and (1/2)^2 = 1/4,
    # prior-weighted 27/256 vs 1/6 = 128/768 vs 81/768 -> 81/209, 128/209.
    payload = symmetric_payload()
    payload["candidates"][0]["atom_groups"][0]["count"] = 2
    payload["candidates"][0]["atom_groups"][0]["isotopes"] = [
        {"shift_mda": 0, "abundance": "3/4"},
        {"shift_mda": 1, "abundance": "1/4"},
    ]
    payload["candidates"][0]["prior"] = "1/3"
    payload["candidates"][1]["prior"] = "2/3"
    payload["ion_count"] = 2
    payload["observations"] = [{"shift_mda": 0, "count": 2}]
    payload["calibration_radius"] = 0
    resp = client.post("/v1/candidate-posterior", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "possible": True,
        "joint_posteriors": [["81/209"], ["128/209"]],
        "candidate_posteriors": ["81/209", "128/209"],
        "offset_posteriors": ["1/1"],
    }


def test_offset_only_reachable_with_calibration():
    # Both candidates are deterministic at shift 0; the observation at
    # shift 2 is explained only by delta = 2, so the offset posterior is a
    # point mass while the candidate posterior keeps the priors.
    candidate = {
        "atom_groups": [
            {
                "count": 1,
                "isotopes": [
                    {"shift_mda": 0, "abundance": "1/1"},
                    {"shift_mda": 1, "abundance": "0/1"},
                ],
            }
        ]
    }
    payload = {
        "candidates": [
            dict(candidate, prior="1/3"),
            dict(candidate, prior="2/3"),
        ],
        "ion_count": 1,
        "observations": [{"shift_mda": 2, "count": 1}],
        "calibration_radius": 2,
    }
    resp = client.post("/v1/candidate-posterior", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "possible": True,
        "joint_posteriors": [
            ["0/1", "0/1", "0/1", "0/1", "1/3"],
            ["0/1", "0/1", "0/1", "0/1", "2/3"],
        ],
        "candidate_posteriors": ["1/3", "2/3"],
        "offset_posteriors": ["0/1", "0/1", "0/1", "0/1", "1/1"],
    }


def test_zero_prior_candidate_gets_zero_posterior():
    payload = symmetric_payload()
    payload["candidates"][0]["prior"] = "1/1"
    payload["candidates"][1]["prior"] = "0/1"
    payload["observations"] = [{"shift_mda": 0, "count": 1}]
    payload["calibration_radius"] = 0
    resp = client.post("/v1/candidate-posterior", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "possible": True,
        "joint_posteriors": [["1/1"], ["0/1"]],
        "candidate_posteriors": ["1/1", "0/1"],
        "offset_posteriors": ["1/1"],
    }


def test_zero_count_observation_has_no_effect():
    payload = symmetric_payload()
    payload["ion_count"] = 2
    with_extra = symmetric_payload()
    with_extra["ion_count"] = 2
    with_extra["observations"].append({"shift_mda": 7, "count": 0})
    base = client.post("/v1/candidate-posterior", json=payload).json()
    extra = client.post("/v1/candidate-posterior", json=with_extra).json()
    assert base == extra


def test_impossible_observation_returns_not_possible():
    # No candidate can produce measured shift 5 even with delta in [-1, 1].
    payload = symmetric_payload()
    payload["observations"] = [{"shift_mda": 5, "count": 1}]
    resp = client.post("/v1/candidate-posterior", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {"possible": False}


def test_counts_consuming_all_ions_leave_empty_other():
    # Both ions listed: the "other" count is 0, so unlisted mass is ignored.
    # Candidate 0 likelihood: sequences (0, 1) and (1, 0) -> 1/2; candidate 1
    # cannot produce shift 1 at delta 0 -> 0.
    payload = symmetric_payload()
    payload["ion_count"] = 2
    payload["observations"] = [
        {"shift_mda": 0, "count": 1},
        {"shift_mda": 1, "count": 1},
    ]
    payload["calibration_radius"] = 0
    resp = client.post("/v1/candidate-posterior", json=payload)
    assert resp.status_code == 200
    assert resp.json() == {
        "possible": True,
        "joint_posteriors": [["1/1"], ["0/1"]],
        "candidate_posteriors": ["1/1", "0/1"],
        "offset_posteriors": ["1/1"],
    }


def test_randomized_matches_brute_force():
    rng = random.Random(20260926)
    for _ in range(20):
        n_candidates = rng.randint(2, 3)
        cuts = sorted(rng.randint(0, 12) for _ in range(n_candidates - 1))
        parts = [b - a for a, b in zip([0] + cuts, cuts + [12])]
        candidates = []
        for part in parts:
            groups = []
            for _ in range(rng.randint(1, 2)):
                shifts = rng.sample(range(0, 5), 2)
                split = rng.randint(1, 23)
                groups.append(
                    {
                        "count": rng.randint(1, 2),
                        "isotopes": [
                            {"shift_mda": shifts[0], "abundance": f"{split}/24"},
                            {
                                "shift_mda": shifts[1],
                                "abundance": f"{24 - split}/24",
                            },
                        ],
                    }
                )
            candidates.append(
                {"atom_groups": groups, "prior": f"{part}/12"}
            )
        ion_count = rng.randint(1, 3)
        obs_shifts = rng.sample(range(0, 9), rng.randint(1, 3))
        counts = [0] * len(obs_shifts)
        for _ in range(rng.randint(0, ion_count)):
            counts[rng.randrange(len(counts))] += 1
        payload = {
            "candidates": candidates,
            "ion_count": ion_count,
            "observations": [
                {"shift_mda": shift, "count": count}
                for shift, count in zip(obs_shifts, counts)
            ],
            "calibration_radius": rng.randint(0, 2),
        }
        resp = client.post("/v1/candidate-posterior", json=payload)
        assert resp.status_code == 200
        body = resp.json()
        expected = brute_force_posterior(payload)
        if expected is None:
            assert body == {"possible": False}
            continue
        joint, candidate_marginal, offset_marginal = expected
        assert body["possible"] is True
        assert body["joint_posteriors"] == [
            [fmt(value) for value in row] for row in joint
        ]
        assert body["candidate_posteriors"] == [
            fmt(value) for value in candidate_marginal
        ]
        assert body["offset_posteriors"] == [
            fmt(value) for value in offset_marginal
        ]
        check_invariants(body, n_candidates, payload["calibration_radius"])


def test_supports_max_scale_candidates():
    # Largest admissible candidates (total count 80, max total shift 2000)
    # with the largest ion_count, calibration radius and observation count.
    group = {
        "count": 40,
        "isotopes": [
            {"shift_mda": 0, "abundance": "1/3"},
            {"shift_mda": 1, "abundance": "1/3"},
            {"shift_mda": 25, "abundance": "1/3"},
        ],
    }
    payload = {
        "candidates": [
            {"atom_groups": [group, group], "prior": "1/2"},
            {"atom_groups": [group, group], "prior": "1/2"},
        ],
        "ion_count": 80,
        "observations": [
            {"shift_mda": 3, "count": 1},
            {"shift_mda": 50, "count": 2},
            {"shift_mda": 100, "count": 1},
        ],
        "calibration_radius": 3,
    }
    resp = client.post("/v1/candidate-posterior", json=payload)
    assert resp.status_code == 200
    body = resp.json()
    assert body["possible"] is True
    check_invariants(body, 2, 3)
    # Identical candidates stay indistinguishable; the long exact fractions
    # are printed in full.
    assert (
        body["candidate_posteriors"][0] == body["candidate_posteriors"][1]
    )
    longest = max(
        len(text) for row in body["joint_posteriors"] for text in row
    )
    assert longest > 100


@pytest.mark.parametrize(
    "mutation",
    [
        # candidates list size / missing
        lambda p: p.update(candidates=p["candidates"][:1]),
        lambda p: p.update(candidates=(p["candidates"] * 3)[:5]),
        lambda p: p.pop("candidates"),
        # prior violations
        lambda p: p["candidates"][0].update(prior="1/3"),  # sum != 1
        lambda p: p["candidates"][0].update(prior="3/2"),
        lambda p: p["candidates"][0].update(prior="1/0"),
        lambda p: p["candidates"][0].update(prior="abc"),
        lambda p: p["candidates"][0].update(prior=0.5),
        lambda p: p["candidates"][0].pop("prior"),
        lambda p: p["candidates"][0].update(extra=1),
        lambda p: p["candidates"][0].pop("atom_groups"),
        # atom_groups violations reused verbatim from /v1/envelope
        lambda p: p["candidates"][0]["atom_groups"][0]["isotopes"][1].update(
            abundance="1/3"
        ),
        lambda p: p["candidates"][0]["atom_groups"][0]["isotopes"][1].update(
            shift_mda=0
        ),
        lambda p: p["candidates"][0]["atom_groups"][0].update(count=0),
        lambda p: p["candidates"][0].update(atom_groups=[]),
        # ion_count out of range / wrong type / missing
        lambda p: p.update(ion_count=0),
        lambda p: p.update(ion_count=81),
        lambda p: p.update(ion_count=1.5),
        lambda p: p.update(ion_count="2"),
        lambda p: p.pop("ion_count"),
        # observations list size / missing / duplicate shifts
        lambda p: p.update(observations=[]),
        lambda p: p.update(observations=p["observations"] * 9),
        lambda p: p.pop("observations"),
        lambda p: p["observations"].append({"shift_mda": 1, "count": 0}),
        # observation field violations
        lambda p: p["observations"][0].update(shift_mda=-1),
        lambda p: p["observations"][0].update(shift_mda=2001),
        lambda p: p["observations"][0].update(shift_mda=1.5),
        lambda p: p["observations"][0].update(shift_mda="1"),
        lambda p: p["observations"][0].pop("shift_mda"),
        lambda p: p["observations"][0].update(count=-1),
        lambda p: p["observations"][0].update(count=1.5),
        lambda p: p["observations"][0].update(count="1"),
        lambda p: p["observations"][0].pop("count"),
        lambda p: p["observations"][0].update(extra=1),
        # counts sum above ion_count
        lambda p: p["observations"][0].update(count=2),
        # calibration_radius out of range / wrong type / missing
        lambda p: p.update(calibration_radius=-1),
        lambda p: p.update(calibration_radius=4),
        lambda p: p.update(calibration_radius=1.5),
        lambda p: p.update(calibration_radius="1"),
        lambda p: p.pop("calibration_radius"),
        # unknown top-level field
        lambda p: p.update(extra=1),
    ],
)
def test_invalid_payloads_return_422(mutation):
    payload = symmetric_payload()
    mutation(payload)
    resp = client.post("/v1/candidate-posterior", json=payload)
    assert resp.status_code == 422


def test_malformed_json_returns_422():
    resp = client.post(
        "/v1/candidate-posterior",
        content="{not json",
        headers={"Content-Type": "application/json"},
    )
    assert resp.status_code == 422
