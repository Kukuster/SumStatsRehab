from functools import lru_cache
from math import sqrt as math_sqrt
from typing import Tuple

import numpy as np
from numpy import floor as np_floor, ceil as np_ceil
from scipy.stats import norm as normal_distribution, binom as binomial_distribution # type: ignore


@lru_cache(10000)
def normal_z_score(p: float = 0.95) -> float:
    """
    Answers the question: "At what point in terms of number of standard deviations from the mean
    the area spanning from -inf to that point on the normal distribution cover `p` of its area

    https://en.wikipedia.org/wiki/Standard_score#/media/File:Z_score_for_Students_A.png
    https://mat117.wisconsin.edu/wp-content/uploads/2014/12/Sec03.-z-score-5.png
    https://cdn1.byjus.com/wp-content/uploads/2017/09/word-image2.png
    """
    return normal_distribution.ppf(p)

@lru_cache(10000)
def normal_z_score_two_tailed(p: float = 0.95) -> float:
    """
    Answers the question: "How many standard deviations from the mean we have to span
    equally in both sides in a normal distribution to cover `p` of the area"

    https://www.freecodecamp.org/news/content/images/2020/08/normal_dist_68_rule.jpg
    https://upload.wikimedia.org/wikipedia/commons/thumb/8/8c/Standard_deviation_diagram.svg/1920px-Standard_deviation_diagram.svg.png
    """
    return normal_distribution.ppf(1-(p/2)) # (1+p)/2  is a more concervative formula

@lru_cache(10000)
def normal_p_area(z: float) -> float:
    """
    Answers the question: "How much area under the normal distribution curve does
    the range from -inf to `z` cover (`z` - number of standard deviations from the mean)"

    https://en.wikipedia.org/wiki/Standard_score#/media/File:Z_score_for_Students_A.png
    https://mat117.wisconsin.edu/wp-content/uploads/2014/12/Sec03.-z-score-5.png
    https://cdn1.byjus.com/wp-content/uploads/2017/09/word-image2.png
    """
    return normal_distribution.cdf(z)

@lru_cache(10000)
def normal_p_area_two_tailed(z: float) -> float:
    """
    Answers the question: "How much area does the range cover that spans
    `z` standard deviations from the mean in equally both sides in a normal distribution"

    https://www.freecodecamp.org/news/content/images/2020/08/normal_dist_68_rule.jpg
    https://upload.wikimedia.org/wikipedia/commons/thumb/8/8c/Standard_deviation_diagram.svg/1920px-Standard_deviation_diagram.svg.png
    """
    return (-normal_distribution.cdf(z)+1) * 2  # (2*z)-1  is a more concervative formula


