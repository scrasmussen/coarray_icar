import matplotlib.pyplot as plt
import pandas as pd
import sys

if (len(sys.argv) < 2):
    sys.exit("Error: too few arguments for `plot.py [graph_data.txt]`")


plt.rc('font', family='serif')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2
# --- be sure to add ---
# AND fig = plt.figure(figsize=(4,3))
# AND plt.tight_layout() before plt.show()


# ---- read input data ----
f = open(sys.argv[1])
l = f.readline()
dims = [int(n) for n in l.split()]
nx = dims[0]; nz = dims[1]; ny = dims[2]
ximages = dims[3]; yimages = dims[4]
dz_value = 500
time=0

nx += 1
ny += 1
nz += 1

header = ['image','timestep','identifier', 'exists','moved', 'x', 'y', 'z',
          'u', 'v', 'w', 'z_meters', 'z_interface', 'pressure','temperature',
          'potential_temperature', 'velocity', 'water_vapor','cloud_water',
          'relative_humidity']


particles = pd.read_csv(f, sep='\s+',header=None, names=header)
num_t = particles['timestep'].max()
unique_particles = particles.identifier.unique()
num_particles = len(list(particles.identifier.unique()))
