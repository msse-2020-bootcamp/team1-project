"""
Tests for Python Standard Library implementation of MC code.
"""
import os
import math
import random
import pytest

import monte_carlo as mc

def test_calculate_LJ():
    rij = 1 
    expected_value = 0

    observed_value = mc.calculate_LJ(rij)
    assert expected_value == observed_value

def test_calculate_distance():
    point1 = [0, 0, 0]
    point2 = [0, 0, 1]

    expected_value = 1

    observed_value = mc.calculate_distance(point1, point2)

    assert expected_value == observed_value

def test_calculate_total_energy():
    coordinates = [[0, 0, 0] , [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]]
    box_length = 10
    cutoff = 3

    expected_value = -2.031005859375
    observed_value = mc.calculate_total_energy(coordinates, box_length, cutoff)
    
    assert expected_value == observed_value

def test_calculate_tail_correction():
    num_particles = 800
    box_length = 10
    cutoff = 3

    expected_tail_correction = -198.4888837441566
    observed_tail_correction = mc.calculate_tail_correction(num_particles, box_length, cutoff)    

    assert observed_tail_correction == expected_tail_correction

def test_read_xyz():
    path = os.path.join('../..', 'lj_sample_configurations', 'lj_sample_config_periodic1.txt')
    expected_coord_len = 800 
    expected_box_len = 10

    calculated_coords, calculated_box_length = mc.read_xyz(path)

    assert expected_coord_len == len(calculated_coords)

def test_accept_or_reject():
    delta_energy = -1
    beta = 1
    assert mc.accept_or_reject(delta_energy, beta) is True

def test_calculate_distance2():
    point_3 = [0, 0, 0]
    point_4 = [0, 1, 1]

    dist2 = mc.calculate_distance(point_3, point_4)
    assert math.sqrt(2) == dist2

def test_calculate_distance3():
    point_1 = [0,  0, 0]
    point_2 = [0, 0, 8]
    box_length = 10

    expected_distance = 2
    observed_distance = mc.calculate_distance(point_1, point_2, box_length=box_length)
    assert expected_distance == observed_distance

def test_accept_or_reject_positive_false():
    try:
        random.seed(0)
        delta_e = 1 
        beta = 1
        assert mc.accept_or_reject(delta_e, beta) is False
    finally:
        random.seed()

def test_calculate_pair_energy():
    coordinates = [[0, 0, 0] , [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]]
    i_particle = 1
    box_length = 10
    cutoff = 3
    expected_value = -2

    calculated_value = mc.calculate_pair_energy(coordinates, i_particle, box_length, cutoff)
    assert math.isclose(calculated_value, expected_value, rel_tol=0.01)

def test_generate_config():
    num_particles = 800
    box_length = 10
    expected_dimensions = 3

    calculated_coords, calculated_box_length = mc.generate_config(num_particles, box_length)

    assert len(calculated_coords) == num_particles

    for e in calculated_coords:
        assert len(e) == expected_dimensions

def test_caclulate_total_energy_NIST():
    path = os.path.join('../..', 'lj_sample_configurations', 'lj_sample_config_periodic1.txt')
    cutoff = 3
    expected_box_length = 10

    expected_energy = -4.3515E+03

    nist_coords, nist_box_length = mc.read_xyz(path)

    calculated_energy = mc.calculate_total_energy(nist_coords, nist_box_length, cutoff)

    assert math.isclose(calculated_energy, expected_energy, rel_tol=0.01)


"""
@pytest.mark.parametrize('variable_name1, variable_name2, ... , variable_nameN', [
    (variable_value1, variable_value2, ...,  variable_valueN)
]
)
def test_function(variable_name1, variable_name2, ... , variable_nameN):
    *** TEST CODE HERE ***
"""

@pytest.mark.parametrize("point1, point2, expected_distance, box_length", 
    [([0, 0, 0], [1, 0, 0], 1, None),
    ([0, 0, 0], [8, 0, 0], 8, None),
    ([0, 0, 0], [8, 0, 0], 2, 10),
    ])
def test_calculate_distance_many(point1, point2, expected_distance, box_length):
    observed_distance = mc.calculate_distance(point1, point2, box_length)
    assert observed_distance == expected_distance