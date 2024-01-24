import random
from numpy import array
from numpy import sum


def bootstrap(data, sample_amount, size):
    samples = []
    for i in range(size):
        sample = []
        for j in range(sample_amount):
            sample.append(random.choice(data))
        sample = array(sample)
        temp = sum(sample, 0)
        samples.append(list(temp))
    return samples


def bootstrap_scores(samples, number_of_kmeres, number_of_filters):
    scores = []
    # calculates float for each value in [hits per filter]
    for i in range(len(samples)):
        score = []
        for j in range(number_of_filters):
            if samples[i][j] == 0:
                score.append(0.0)
            else:
                score.append(round(float(samples[i][j]) / float(number_of_kmeres), 2))
        scores.append(score)
    return scores
