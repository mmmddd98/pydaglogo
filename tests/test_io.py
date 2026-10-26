from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy.stats.contingency import odds_ratio

matplotlib.use("Agg")

from pydaglogo import (
    build_background_model_from_samples,
    build_windows_from_flanks,
    build_windows_from_positions,
    format_sequence,
    make_daglogo,
    plot_logo,
    Proteome,
    run_daglogo,
    test_dau,
)


def test_build_windows_from_flanks_filters_invalid_noncanonical_rows():
    frame = pd.DataFrame(
        {
            "five_res_before": ["ABCDE", "ABXDE", None],
            "five_res_after": ["FGHIJ", "KLMNP", "QRSTU"],
        }
    )

    result = build_windows_from_flanks(frame)

    assert list(result["window"]) == ["ACDEFGHI"]
    assert len(result) == 1


def test_build_windows_from_flanks_can_keep_invalid_rows_without_canonical_filter():
    frame = pd.DataFrame(
        {
            "five_res_before": [" abcdE ", "VWXYZ"],
            "five_res_after": ["fghiJ ", "QRSTU"],
        }
    )

    result = build_windows_from_flanks(
        frame,
        canonical_only=False,
        drop_invalid=False,
        output_col="custom_window",
    )

    assert list(result["custom_window"]) == ["BCDEFGHI", "WXYZQRST"]
    assert len(result) == 2


def test_build_windows_from_positions_filters_invalid_noncanonical_rows():
    frame = pd.DataFrame(
        {
            "P4": ["D", "D", "D"],
            "P3": ["E", "E", "E"],
            "P2": ["V", "X", None],
            "P1": ["D", "D", "D"],
            "P1'": ["G", "G", "G"],
            "P2'": ["A", "A", "A"],
            "P3'": ["A", "A", "A"],
            "P4'": ["A", "A", "A"],
        }
    )

    result = build_windows_from_positions(frame)

    assert list(result["window"]) == ["DEVDGAAA"]
    assert len(result) == 1


def test_format_sequence_supports_unmapped_mode_without_proteome():
    result = format_sequence(
        ["ACDEFGHI", "LMNPQRST"],
        proteome=None,
        upstream_offset=3,
        downstream_offset=4,
        map_to_proteome=False,
    )

    assert result.type == "unmapped"
    assert list(result.data["IDs"]) == ["UNMAPPED", "UNMAPPED"]
    assert list(result.data["peptide"]) == ["ACDEFGHI", "LMNPQRST"]


def test_format_sequence_reports_mapped_anchor_position():
    proteome = Proteome(
        proteome=pd.DataFrame(
            {
                "SEQUENCE": ["XXACDEFGHIYY"],
                "ID": ["protein-1"],
                "LEN": [12],
            }
        ),
        type="fasta",
        species="test",
    )

    result = format_sequence(
        ["ACDEFGHI"],
        proteome=proteome,
        upstream_offset=3,
        downstream_offset=4,
        map_to_proteome=True,
    )

    assert result.data.loc[0, "anchorAA"] == "E"
    assert result.data.loc[0, "anchorPos"] == 6


def test_format_sequence_rejects_noncanonical_window_characters():
    with pytest.raises(ValueError, match="non-canonical"):
        format_sequence(
            ["ACX"],
            proteome=None,
            upstream_offset=1,
            downstream_offset=1,
            map_to_proteome=False,
        )


def test_run_daglogo_infers_left_of_center_anchor_for_even_width_windows():
    windows = ["DEVDGAAA", "DEVDSAAA", "DAVDGAAA"]
    background = build_background_model_from_samples(
        [
            pd.DataFrame(
                [
                    list("DEVDGAAA"),
                    list("DEVDSAAA"),
                    list("DAVDGAAA"),
                ],
                columns=[f"V{i}" for i in range(1, 9)],
            ),
            pd.DataFrame(
                [
                    list("DEVDGAAA"),
                    list("DAVDGAAA"),
                    list("DAVDSAAA"),
                ],
                columns=[f"V{i}" for i in range(1, 9)],
            ),
        ],
        test_type="ztest",
    )

    result = run_daglogo(
        windows,
        background_model=background,
        input_mode="windows",
        map_to_proteome=False,
    )

    assert result.upstream_offset == 3
    assert result.downstream_offset == 4
    assert result.test_type == "ztest"


def test_make_daglogo_supports_position_dataframe_input():
    frame = pd.DataFrame(
        {
            "P4": ["D", "D", "D"],
            "P3": ["E", "E", "A"],
            "P2": ["V", "V", "V"],
            "P1": ["D", "D", "D"],
            "P1'": ["G", "S", "G"],
            "P2'": ["A", "A", "A"],
            "P3'": ["A", "A", "A"],
            "P4'": ["A", "A", "A"],
        }
    )
    background = build_background_model_from_samples(
        [
            pd.DataFrame(
                [
                    list("DEVDGAAA"),
                    list("DEVDSAAA"),
                    list("DAVDGAAA"),
                ],
                columns=[f"V{i}" for i in range(1, 9)],
            ),
            pd.DataFrame(
                [
                    list("DEVDGAAA"),
                    list("DAVDGAAA"),
                    list("DAVDSAAA"),
                ],
                columns=[f"V{i}" for i in range(1, 9)],
            ),
        ],
        test_type="ztest",
    )

    result, ax, matrix = make_daglogo(
        frame,
        background_model=background,
        input_mode="positions",
        labels=["P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"],
        map_to_proteome=False,
        pvalue_cutoff=None,
    )

    assert result.upstream_offset == 3
    assert ax is not None
    assert ax.get_ylabel() == "Δ residue frequency"
    assert matrix.shape == (8, 20)


def test_make_daglogo_supports_flank_dataframe_input():
    frame = pd.DataFrame(
        {
            "five_res_before": ["ADEVD", "ADEVD", "ADEVD"],
            "five_res_after": ["GAAAK", "SAAAK", "GSSSK"],
        }
    )
    background = build_background_model_from_samples(
        [
            pd.DataFrame(
                [
                    list("AQADGAAA"),
                    list("DKAASAAA"),
                    list("NEADGASA"),
                    list("DRAAGAAA"),
                ],
                columns=[f"V{i}" for i in range(1, 9)],
            ),
            pd.DataFrame(
                [
                    list("AEADGAAA"),
                    list("NQVASAAA"),
                    list("DEADGAAA"),
                    list("GKVAGAAA"),
                ],
                columns=[f"V{i}" for i in range(1, 9)],
            ),
        ],
        test_type="ztest",
    )

    result, ax, matrix = make_daglogo(
        frame,
        background_model=background,
        input_mode="flanks",
        labels=["P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'"],
        map_to_proteome=False,
        pvalue_cutoff=None,
    )

    assert result.upstream_offset == 3
    assert ax is not None
    assert ax.get_ylabel() == "Δ residue frequency"
    assert matrix.shape == (8, 20)


def test_plot_logo_allows_custom_y_axis_label():
    matrix = pd.DataFrame(
        [
            {"A": 0.2, "C": 0.0},
            {"A": 0.0, "C": 0.1},
        ]
    )

    ax, _ = plot_logo(matrix, y_label="Custom label")

    assert ax.get_ylabel() == "Custom label"
    plt.close(ax.figure)


def test_plot_logo_can_suppress_y_axis_label():
    matrix = pd.DataFrame(
        [
            {"A": 0.2, "C": 0.0},
            {"A": 0.0, "C": 0.1},
        ]
    )

    ax, _ = plot_logo(matrix, y_label=None)

    assert ax.get_ylabel() == ""
    plt.close(ax.figure)


def test_build_background_model_from_samples_requires_single_matrix_for_fisher():
    sample = pd.DataFrame([["A", "C"], ["D", "E"]], columns=["V1", "V2"])

    try:
        build_background_model_from_samples(
            [sample, sample],
            test_type="fisher",
        )
        raised = False
    except ValueError as exc:
        raised = True
        assert "exactly one background matrix" in str(exc)

    assert raised


def test_test_dau_runs_for_fisher_with_single_background_matrix(monkeypatch):
    monkeypatch.setenv("PYDAGLOGO_FISHER_BACKEND", "scipy")
    dag_peptides = format_sequence(
        ["ACD", "AED", "ACD"],
        proteome=None,
        upstream_offset=1,
        downstream_offset=1,
        map_to_proteome=False,
    )
    background = build_background_model_from_samples(
        [pd.DataFrame([["A", "C", "D"], ["A", "E", "D"], ["A", "C", "D"]], columns=["V1", "V2", "V3"])],
        test_type="fisher",
    )

    result = test_dau(dag_peptides, background, grouping_scheme="no")

    assert result.test_type == "fisher"
    assert result.difference.shape == (20, 3)
    assert result.statistics.shape == (20, 3)
    assert result.pvalue.shape == (20, 3)


def test_test_dau_rejects_out_of_range_bg_noise():
    dag_peptides = format_sequence(
        ["ACD", "AED", "ACD"],
        proteome=None,
        upstream_offset=1,
        downstream_offset=1,
        map_to_proteome=False,
    )
    sample = pd.DataFrame(
        [["A", "C", "D"], ["A", "E", "D"], ["A", "C", "D"]],
        columns=["V1", "V2", "V3"],
    )
    background = build_background_model_from_samples(
        [sample, sample],
        test_type="ztest",
    )

    with pytest.raises(ValueError, match="bg_noise"):
        test_dau(dag_peptides, background, grouping_scheme="no", bg_noise=1.1)


def test_test_dau_fisher_uses_conditional_odds_ratio(monkeypatch):
    monkeypatch.setenv("PYDAGLOGO_FISHER_BACKEND", "scipy")
    dag_peptides = format_sequence(
        ["AAD", "ACD", "ACD"],
        proteome=None,
        upstream_offset=1,
        downstream_offset=1,
        map_to_proteome=False,
    )
    background = build_background_model_from_samples(
        [pd.DataFrame([["A", "A", "D"], ["A", "A", "D"], ["A", "C", "D"]], columns=["V1", "V2", "V3"])],
        test_type="fisher",
    )

    result = test_dau(dag_peptides, background, grouping_scheme="no")

    expected = odds_ratio(np.array([[1, 2], [2, 1]], dtype=np.int64), kind="conditional").statistic
    observed = result.statistics.loc["A", "AA0"]

    assert observed == expected
