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

def removeTopRightFrame(ax=None):
  if ax is None:
    ax = plt.axes()
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)

def hex2rgb(hex):
  return tuple(int(hex[i:i+2], 16) for i in (1, 3, 5))

def rgb2hex(rgb):
  return '#%02x%02x%02x' % rgb

def applyAlpha(color, alpha):
  rgb = hex2rgb(color)
  rgb = tuple(map(lambda r: r*alpha + 255.*(1.-alpha), rgb))
  return rgb2hex(rgb)
  
def lighten(c, a):
  return applyAlpha(c, a)

def darken(color, black):
  rgb = hex2rgb(color)
  rgb = tuple(map(lambda r: r*black, rgb))
  return rgb2hex(rgb)

def pureRGB(color):
  rgb = hex2rgb(color)
  (imin, cmin), (imid, cmid), (imax, cmax) = tuple(sorted(enumerate(rgb), key=lambda x: x[1]))
  c = [0,0,0]
  c[imax] = 255
  c[imid] = int(255.999 * (cmid - cmin) / (cmax - cmin))
  return tuple(c)

def pure(color):
  return rgb2hex(pureRGB(color))

# tools for HSL representation
# (H is implicitly set by choosing a base color, setSL conserves H)

def get_sl(color):
  rgb = hex2rgb(color)
  cmax = max(rgb) / 255.999
  cmin = min(rgb) / 255.999
  l = 0.5 * (cmax + cmin)
  if cmax == cmin:
    s = 0.
  elif l > 0.5:
    s = 0.5 * (cmax - cmin) / (1. - l)
  else:
    s = 0.5 * (cmax - cmin) / l
  return s, l

def set_sl(color, s=None, l=None):
  # fill out missing params
  (s_, l_) = get_sl(color)
  #print("color", hex2rgb(color))
  #print("sl (old)", s_, l_)
  s = s_ if s is None else s
  l = l_ if l is None else l
  # actually set sl
  rgbpure = pureRGB(color)
  #print("pure", rgbpure)
  #print("sl (new)", s, l)
  ll = 1. - l if l > 0.5 else l
  rgb = tuple(map(lambda r: int(255.999 * (l - ll * s * (1. - 2.*r / 255.999))), rgbpure))
  #print("final", rgb)
  return rgb2hex(rgb)

rotateRed = lambda blue: "#" + blue[5:7] + blue[1:5]
rotateGreen = lambda blue: "#" + blue[3:7] + blue[1:3]

class Colors:
  pure = '#0040ff'
  muted = set_sl(pure, .33, .6)
  medium = set_sl(pure, .5, .5)
  intense = set_sl(pure, 1., .37)

  # for three shades in 1 fig
  dark = set_sl(pure, 1., .15)
  mediumintense = set_sl(pure, 1., .6)
  light = set_sl(pure, 1., .85)

  # for 2 shades
  mediumdark = set_sl(pure, .9, .3)
  mediumlight = set_sl(pure, .7, .7)

  # for overlap with .experiment
  lightmediumlight = set_sl(pure, .8, .8)

  purepurple = "#6600ff"
  purple = set_sl("#6600ff", .8, .7)
  
  #light = applyAlpha(intense, 0.5)
  #lightlight = applyAlpha(intense, 0.3)

  #dark = '#061a73'

  #lightmuted = '#9aaae6'
  #darkintense = '#061a73'

  complement = '#c2b342'
  gold = '#c7a708'
  orange = '#cc9c43'
  
  pink = '#be40c7'
  lightpink = applyAlpha(pink, 0.4)

  simulation = medium
  experiment = rotateRed(medium) # #bf3f5f, used to be #bd3c59
  purered = "ff003f" # pure version #ff

  protein = muted
  receptor = gold

  adam1 = '#3069ab'
  adam2 = '#5993d0'
  adam3 = '#93cbf2'



colors = Colors()
