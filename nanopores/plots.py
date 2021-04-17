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

def hex2rgb(hex):
  return tuple(int(hex[i:i+2], 16) for i in (1, 3, 5))
def rgb2hex(rgb):
  return '#%02x%02x%02x' % rgb

def applyAlpha(color, alpha):
  print(color)
  rgb = hex2rgb(color)
  print(rgb)
  rgb = tuple(map(lambda r: r*alpha + 255.*(1.-alpha), rgb))
  print(rgb)
  return rgb2hex(rgb)

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
  lightpink = applyAlpha(pink, 0.4)

  simulation = muted
  experiment = rotateRed(medium)

  protein = muted
  receptor = gold


colors = Colors()
