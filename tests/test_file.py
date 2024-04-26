import pyranges as pr
from pyranges_plot.data_preparation import make_subset


def test_subset():
    # Data has one interval per id
    df_1 = pr.PyRanges(
        {
            "Chromosme": [1, 1, 2, 2, 2],
            "Start": [10, 20, 10, 20, 30],
            "End": [11, 21, 11, 21, 31],
            "transcript_id": ["T1", "T2", "T3", "T4", "T5"],
        }
    )

    result_subset_1, _ = make_subset(df_1, "transcript_id", 3)
    expected_subset_1 = pr.PyRanges(
        {
            "Chromosme": [1, 1, 2],
            "Start": [10, 20, 10],
            "End": [11, 21, 11],
            "transcript_id": ["T1", "T2", "T3"],
        }
    )

    assert len(result_subset_1) == len(expected_subset_1)

    # Data has one id per chromosome
    df_2 = pr.PyRanges(
        {
            "Chromosme": ["a", "b", "c", "d"],
            "Start": [1, 2, 3, 4],
            "End": [10, 20, 30, 40],
            "transcript_id": ["T1", "T2", "T3", "T4"],
        }
    )

    result_subset_2, _ = make_subset(df_2, "transcript_id", 2)
    expected_subset_2 = pr.PyRanges(
        {
            "Chromosme": [
                "a",
                "a",
            ],
            "Start": [1, 2],
            "End": [10, 20],
            "transcript_id": ["T1", "T2"],
        }
    )

    assert len(result_subset_2) == len(expected_subset_2)

    # Data has as many ids as max_shown in the same chromosome
    df_3 = pr.PyRanges(
        {
            "Chromosme": ["a", "a", "a"],
            "Start": [1, 2, 3],
            "End": [10, 20, 30],
            "transcript_id": ["T1", "T2", "T3"],
        }
    )

    result_subset_3, _ = make_subset(df_3, "transcript_id", 3)
    expected_subset_3 = pr.PyRanges(
        {
            "Chromosme": [
                "a",
                "a",
                "a",
            ],
            "Start": [1, 2, 3],
            "End": [10, 20, 30],
            "transcript_id": ["T1", "T2", "T3"],
        }
    )

    assert len(result_subset_3) == len(expected_subset_3)

    # Data has genes with several lines in different chromosomes
    df_4 = pr.PyRanges(
        {
            "Chromosme": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "Start": [i for i in range(9)],
            "End": [i + 10 for i in range(9)],
            "transcript_id": ["T1", "T1", "T1", "T2", "T2", "T3", "T3", "T4", "T4"],
        }
    )

    result_subset_4, _ = make_subset(df_4, "transcript_id", 3)
    expected_subset_4 = pr.PyRanges(
        {
            "Chromosme": [
                1,
                1,
                1,
                1,
                1,
                2,
                2,
            ],
            "Start": [i for i in range(7)],
            "End": [i + 10 for i in range(7)],
            "transcript_id": ["T1", "T1", "T1", "T2", "T2", "T3", "T3"],
        }
    )

    assert len(result_subset_4) == len(expected_subset_4)

    # Data contains fewer ids than max_shown and has genes with several lines
    df_5 = pr.PyRanges(
        {
            "Chromosme": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "Start": [i for i in range(9)],
            "End": [i + 10 for i in range(9)],
            "transcript_id": ["T1", "T1", "T1", "T2", "T2", "T3", "T3", "T4", "T4"],
        }
    )

    result_subset_5, _ = make_subset(df_5, "transcript_id", 5)
    expected_subset_5 = pr.PyRanges(
        {
            "Chromosme": [1, 1, 1, 1, 1, 2, 2, 2, 2],
            "Start": [i for i in range(9)],
            "End": [i + 10 for i in range(9)],
            "transcript_id": ["T1", "T1", "T1", "T2", "T2", "T3", "T3", "T4", "T4"],
        }
    )

    assert len(result_subset_5) == len(expected_subset_5)
