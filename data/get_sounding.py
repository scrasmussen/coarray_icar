import pandas as pd
import numpy as np
import requests
import sys
from bs4 import BeautifulSoup
from datetime import date
from scipy.interpolate import griddata


# ------------------------------------------------------------------------------
# Get today's date and URL setup
# ------------------------------------------------------------------------------
today = date.today()
url   = 'http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&'
year  = 'YEAR=' + today.strftime("%Y")
month = '&MONTH=' + today.strftime("%m")
day_i = today.strftime("%d")
day   = '&FROM=' + day_i + '12&TO=' + day_i + '12'
end   = '&STNM=72469'
url  += year + month + day + end
url = 'http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR=2020&MONTH=10&FROM=0612&TO=0612&STNM=72597'


# ------------------------------------------------------------------------------
# Scrape atmospheric sounding, edit dataframe
# ------------------------------------------------------------------------------
req = requests.get(url).text
soup = BeautifulSoup(req, 'lxml')
result = soup.find('pre').find(text=True).split('\n')
names = result[2:3][0].split()
data_array = result[5:-1]
data = [row.split() for row in data_array]
df = pd.DataFrame(data,columns=names).fillna('-999')

print("Data retrieved for " + today.isoformat())


# ------------------------------------------------------------------------------
# Interpolate and write to file
# ------------------------------------------------------------------------------
def write_data(filename, data, start, end):
    filename = 'sounding/' + filename
    f = open(filename, 'w')
    f.write(str(data.size) +" " + str(start) + " " + str(end) + '\n')
    data.tofile(f, "\n")

def interpolate_and_write(filename,  values, points):
    start=500
    end=15500
    grid_x = np.mgrid[start:end+1:1]
    data = griddata(points, values, grid_x, method='cubic')
    write_data(filename, data, start, end)


n_type = 'float'
points = df.HGHT[1:].to_numpy(n_type)
theta_values = df.THTA[1:].to_numpy(n_type)
pres_values  = df.PRES.to_numpy(n_type)

interpolate_and_write('sounding-potential-temp.txt', theta_values, points)
interpolate_and_write('sounding-pressure.txt', pres_values, points)






# ------------------------------------------------------------------------------
# Save dataframe to file
# ------------------------------------------------------------------------------
# f = open('sounding/sounding-README.txt', 'w')
# f.write("Sounding from the date " + str(date.today()))
# for col in df.columns:
#     # filename = 'sounding/sounding-' + str(date.today()) + '-' + col + '.txt'
#     filename = 'sounding/sounding-' + col + '.txt'
#     f = open(filename, 'w')
#     f.write(str(df.shape[0]) + '\n')
#     f.write(df[col].to_csv(header=False, index=False))

print('File written')
