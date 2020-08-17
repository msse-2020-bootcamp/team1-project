"""
Tests for Python Standard Library implementation of MC code.
"""
import os
import math
import random
import pytest

import numpy as np

import monte_carlo_np as mc

def test_calculate_LJ():
    rij = 1 
    expected_value = 0

    observed_value = mc.calculate_LJ(rij)
    assert expected_value == observed_value

def test_calculate_distance():
    coord_1 = np.array([0, 0, 0])
    coord_2 = np.array([0, 0, 1])

    expected_distance = np.array([1])

    calculated_distance = mc.calculate_distance(coord_1, coord_2)

    assert np.array_equal(expected_distance, calculated_distance)

def test_calculate_total_energy():
    coordinates = np.array([[0, 0, 0] , [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]])
    box_length = 10
    cutoff = 3

    expected_value = -2.031005859375
    calculated_value = mc.calculate_total_energy(coordinates, box_length, cutoff)
    
    assert math.isclose(calculated_value, expected_value, rel_tol=0.01)


def test_accept_or_reject():
    delta_energy = -1
    beta = 1
    assert accept_or_reject(delta_energy, beta) is True

def test_calculate_distance2():
    coord_3 = np.array([[0, 0, 0], [1, 0, 0]])
    coord_4 = np.array([[0, 1, 0], [0, 0,0 ]])

    expected_distance = np.array([1,1])

    calculated_distance = mc.calculate_distance(coord_3, coord_4)

    assert np.array_equal(expected_distance, calculated_distance)
  

def test_calculate_tail_correction():
    num_particles = 800
    box_length = 10
    cutoff = 3

    expected_tail_correction = -198.4888837441566
    observed_tail_correction = mc.calculate_tail_correction(num_particles, box_length, cutoff)    

    assert observed_tail_correction == expected_tail_correction


def test_read_xyz():
    path = os.path.join('..', 'lj_sample_configurations', 'lj_sample_config_periodic1.txt')
    expected_coord_len = 800 
    expected_box_len = 10

    calculated_coords, calculated_box_length = mc.read_xyz(path)

    assert expected_coord_len == len(calculated_coords)


def test_accept_or_reject():
    delta_energy = -1
    beta = 1
    assert mc.accept_or_reject(delta_energy, beta) is True

def test_calculate_distance2():
    coord_3 = np.array([[0, 0, 0], [0, 0, 9]])
    coord_4 = np.array([[0, 0, 9], [0, 0, 0]])

    expected_distance = np.array([9,9])

    calculated_distance = mc.calculate_distance(coord_3, coord_4)

    assert np.array_equal(expected_distance, calculated_distance)

def test_calculate_distance3():
    coord_1 = np.array([[0, 0, 0], [0,0,2]])
    coord_2 = np.array([0, 0, 5])
    box_length = 10
    
    expected_distance = np.array([5, 3])
    
    calculated_distance = mc.calculate_distance(coord_1, coord_2, box_length=10)
    
    assert np.array_equal(expected_distance, calculated_distance)

def test_accept_or_reject_positive_false():
    try:
        random.seed(0)
        delta_e = 1 
        beta = 1
        assert mc.accept_or_reject(delta_e, beta) is False
    finally:
        random.seed()

def test_calculate_pair_energy():
    coordinates = np.array([[0, 0, 0] , [0, 0, 2**(1/6)], [0, 0, 2*(2**(1/6))]])
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
    path = os.path.join('..', 'lj_sample_configurations', 'lj_sample_config_periodic1.txt')
    cutoff = 3
    expected_box_length = 10

    expected_energy = -4.3515E+03


    nist_coords, nist_box_length = mc.read_xyz(path)

    nist_coords = np.array(nist_coords)

    calculated_energy = mc.calculate_total_energy(nist_coords, nist_box_length, cutoff)

    assert math.isclose(calculated_energy, expected_energy, rel_tol=0.01)



