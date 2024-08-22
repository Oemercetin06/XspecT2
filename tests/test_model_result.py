from xspect.models.result import ModelResult


def test_get_scores():
    """Test the get_scores method of ModelResult."""
    # Create a ModelResult instance with sample data
    hits = {
        "subsequence1": {"label1": 10, "label2": 5},
        "subsequence2": {"label1": 8, "label2": 3},
    }
    num_kmers = {"subsequence1": 100, "subsequence2": 50}
    model_result = ModelResult("test_slug", hits, num_kmers)

    # Expected scores
    expected_scores = {
        "subsequence1": {"label1": 0.1, "label2": 0.05},
        "subsequence2": {"label1": 0.16, "label2": 0.06},
        "total": {"label1": 0.12, "label2": 0.05},
    }

    # Get the scores from the model_result
    scores = model_result.get_scores()

    # Assert the scores are as expected
    assert scores == expected_scores
