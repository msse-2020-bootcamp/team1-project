"""
Tests for python standard library implementation of MC code.
"""

from monte_carlo import calculate_distance
from monte_carlo import accept_or_reject
import pytest

def test_calculate_distance():

    point1 = [0,0,0]
    point2 = [0,0,1]

    expected_value = 1

    observed_value = calculate_distance(point1, point2)

    assert expected_value == observed_value


def test_accept_or_reject():
    delta_energy = -1
    beta = 1

    assert accept_or_reject(delta_energy, beta) is True


@pytest.mark.parametrize("point1, point2, expected_distance, box_longth", [
    ([0,0,0], [1,0,0], 1, None),
    ([0,0,0], [8,0,0], 8, None),
    ([0,0,0], [8,0,0], 2, 10),
])
def test_calculate_distance_many(point1, point2, expected_distance, box_longth):
    observed_distance = calculate_distance(point1, point2, box_longth)
    
    assert expected_distance == observed_distance
    