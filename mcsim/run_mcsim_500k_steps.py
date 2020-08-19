import monte_carlo as mc
import os 

num_steps = 500000

path = os.path.join('..', 'lj_sample_configurations', 'lj_sample_config_periodic1.txt')
coords, box_length = mc.read_xyz(path)
mc.run_simulation(coords, box_length, 3, 1.4, num_steps)