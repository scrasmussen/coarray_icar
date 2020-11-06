from mpl_toolkits import mplot3d
from matplotlib.ticker import MaxNLocator
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import sys

plt.rc('font', family='serif')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 2



if (len(sys.argv) < 2):
    sys.exit("Error: too few arguments for `plot.py [graph_data.txt]`")

# ---- read input data ----
f = open(sys.argv[1])
header = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time']
header = ['nx','nz','ny','np','x_images','y_images','n_particles','timesteps',
          'time', 'is_dry']
df = pd.read_csv(f, sep='\s+',header=None, names=header)


# --- order data ---
df.sort_values('n_particles', inplace=True)


dry = df[df.is_dry == 'T']
label = "With dry particles"
plt.plot(dry.n_particles, dry.time, marker = '.',
         label=label)

wet = df[df.is_dry == 'F']
label = "With saturated particles"
plt.plot(wet.n_particles, wet.time, marker = '.',
         label=label)

# base line
baseline = df[df.is_dry == 'B']
plt.plot([0,df.n_particles.max()], [baseline.time,baseline.time], marker = '.',
             label="Baseline: particles turned off", color='red')


# plt.legend(title="Dimensions")
plt.legend()
plt.xlabel("number of particles per image (million)")
plt.ylabel("time (seconds)")

plt.title("Particle Scaling, 200 timesteps, 44 images")

plt.show()
sys.exit()

def format_xaxis(value, tick_num):
    return int(value / 1000000.0)

ax = plt.gca()
ax.ticklabel_format(axis="x",style='plain')
ax.xaxis.set_major_formatter(plt.FuncFormatter(format_xaxis))
plt.title("Particle Scaling, 500 timesteps, 16 images")

model = LinearRegression()
model.fit(df.n_particles.values.reshape(-1,1), df.time.values.reshape(-1,1))
print("Linear Regression: y = " +str(model.coef_[0][0])+"x + " +
      str(model.intercept_[0]))
r_squared = model.score(df.n_particles.values.reshape(-1,1),
                      df.time.values.reshape(-1,1))
print("Coefficient of determination R^2 = ", str(r_squared))
# plt.plot(df.n_particles, model.predict(df.n_particles.values.reshape(-1,1)), color='red')


a = model.coef_[0][0]
b = model.intercept_[0]
y = a * 40000000 + b
time = 12 * 60 * 60
print("y = ", str(y), "data = ", time, "error =",(time-y)/y)

# plt.xscale('log', base=2)
# ax.set_xticklabels([])
# ax.set_yticklabels([])

plt.show()
print("Fin!")
