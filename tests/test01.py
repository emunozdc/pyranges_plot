import pyranges as pr

from pyranges_plot.data_preparation import make_subset


def test_correct_subset():
    df_1 = pr.PyRanges(
        {
            "Chromosme": [1, 1, 2, 2, 2],
            "Start": [10, 20, 10, 20, 30],
            "End": [11, 21, 11, 21, 31],
            "transcript_id": ["T1", "T2", "T3", "T4", "T5"],
        }
    )

    result_subset_1, _ = make_subset(df_1, "transcript_id", 3)
    expected_subset_1, _ = pr.PyRanges(
        {
            "Chromosme": [1, 1, 2],
            "Start": [10, 20, 10],
            "End": [11, 21, 11],
            "transcript_id": ["T1", "T2", "T3"],
        }
    )

    assert len(result_subset_1) == len(expected_subset_1)
