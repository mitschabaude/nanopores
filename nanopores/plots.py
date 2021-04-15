from matplotlib import rcParams, rc
rcParams.update({
    "font.size" : 7,
    "axes.titlesize" : 7,
    "font.family" : "sans-serif",
    "font.sans-serif" : ["CMU Sans Serif"],
    "lines.linewidth" : 1,
    "lines.markersize" : 4,
})
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def addMinorTicks(x=2):
  ax = plt.axes()
  ax.xaxis.set_minor_locator(AutoMinorLocator(x))
  ax.yaxis.set_minor_locator(AutoMinorLocator(x))

def removeTopRightFrame():
  ax = plt.axes()
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)

rotateRed = lambda blue: "#" + blue[5:7] + blue[1:5]
rotateGreen = lambda blue: "#" + blue[3:7] + blue[1:3]

class Colors:
  muted = '#788abd'
  medium = '#3c59bd'
  intense = '#0224bd'

  lightmuted = '#9aaae6'
  darkintense = '#061a73'

  complement = '#c2b342'
  gold = '#c7a708'
  orange = '#cc9c43'
  
  pink = '#be40c7'

  simulation = medium
  experiment = rotateRed(simulation)

  protein = muted
  receptor = gold


colors = Colors()
