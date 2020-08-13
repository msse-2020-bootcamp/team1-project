import pytest

import MC_w_raise_exception_camille as mc

def test_run_sim_exception():
    with pytest.raises(ZeroDivisionError):
        coordinates = [[0, 0, 0], [0, 0, 0]]
        box_length = 10
        cutoff = 3
        reduced_temperature = 1.4
        num_steps = 10000
        mc.run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps)

def test_run_sim_no_exception():
    #test will fail if unexpected exception occurs
    coordinates = [[1, 0, 0], [0, 0, 0]]
    box_length = 10
    cutoff = 3
    reduced_temperature = 1.4
    num_steps = 10000
    mc.run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps)

        