import pyranges as pr

from pyranges_plot.data_preparation import make_subset


def subset_test():
    df = pr.from_dict(
        {
            "Chromosme": [1, 1, 2, 2, 2],
            "Start": [10, 20, 10, 20, 30],
            "End": [11, 21, 11, 21, 31],
            "transcript_id": ["T1", "T2", "T3", "T4", "T5"],
        }
    ).df

    result_subset, _ = make_subset(df, "transcript_id", 3)
    expected_subset, _ = pr.from_dict(
        {
            "Chromosme": [1, 1, 2],
            "Start": [10, 20, 10],
            "End": [11, 21, 11],
            "transcript_id": ["T1", "T2", "T3"],
        }
    ).df

    assert len(result_subset) == len(expected_subset)
