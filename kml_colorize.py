from simplekml import Kml, Style, Color
import numpy as np
import pandas as pd
import utm

import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt

#Create color wheel - https://bsou.io/posts/color-gradients-with-python
def RGB_to_hex(RGB):
  ''' [255,255,255] -> "#FFFFFF" '''
  # Components need to be integers for hex to make sense
  RGB = [int(x) for x in RGB]
  return "#"+"".join(["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in RGB])
  
def hex_to_RGB(hex):
  ''' "#FFFFFF" -> [255,255,255] '''
  # Pass 16 to the integer function for change of base
  return [int(hex[i:i+2], 16) for i in range(1,6,2)]
  
def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}

def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return color_dict(RGB_list)  
   
def rgb(minimum, maximum, value):
    ''' Alternate way to set rgb values on a range, here with green set always to zero '''
    minimum, maximum = float(minimum), float(maximum)
    ratio = 2 * (value-minimum) / (maximum - minimum)
    b = int(max(0, 255*(1 - ratio)))
    r = int(max(0, 255*(ratio - 1)))
    #g = 255 - b - r
    g = 0
    return r, g, b
    
    
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    import matplotlib.colors as mcolors
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

#Load Data
base = '/media/tsmith/DATAPART1/Dropbox/Scripts/Bodo_KML/'
base_fid = base + 'R16_projected_r1.0m_centerline_smooth_10.csv'
df = pd.read_csv(base_fid)
  
min_col = df[' ChannelWidth_m'].min() #Get the minimum channel width
max_col = df[' ChannelWidth_m'].max() #Max channel width

cols = np.linspace(min_col, max_col, 20) #Divide the channels into 20 bins
RGB = linear_gradient('#4682B4','#FFB347', 20) #Get 20 RGB/hex codes on a linear gradient between two hex codes
#NOTE that RGB is a dictionary with 4 lists: hex, r, g, and b. However, the html hex this generates doesn't play well with kml hex, which is defined differently

 #This takes a normal matplotlib colormap and turns it into tuples on the linspace [r, g, b, alpha] x 20
color = cm.viridis(np.linspace(0,0.9,20))*255


#This is a more customizable version that lets you choose from named matplotlib colors on an interval
c = colors.ColorConverter().to_rgb #Make a mapper from names to rgb tuples
rvb = make_colormap([c('midnightblue'), c('cornflowerblue'), c('lavender'), 0.495, c('black'), 0.505, c('bisque'), c('salmon'), c('indianred')]) #Make a colormap from the color on the left to the one on the far right passing through each of the named colors. You can also put in value ranges, as I have here, to make everything right around 0.5 be 'black'. Otherwise it runs from datamin to datamax
#color = rvb(np.linspace(0,0.9,20))*255 #You then do the same thing to get a set of tuples

#Create KML
kml = Kml()
fol = kml.newfolder(name="Channel Width")

for i in range(len(cols)-1):
    minc = cols[i] #Get the bin min and max
    maxc = cols[i+1]
    sharedstyle = Style()
    #sharedstyle.iconstyle.color = Color.rgb(RGB['r'][i],RGB['g'][i],RGB['b'][i]) #Define color from linear gradient from HEX codes
    sharedstyle.iconstyle.color = Color.rgb(color[i][0],color[i][1],color[i][2]) #Define color from matplotlib colorbar
    #r, g, b = rgb(min_col, max_col, minc) #Get RGB from alternate function
    #sharedstyle.iconstyle.color = Color.rgb(r, g, b) #Set RGB from alternate function
    #sharedstyle.iconstyle.color = RGB['hex'][i]
    sharedstyle.iconstyle.icon.href = 'https://maps.google.com/mapfiles/kml/pal2/icon26.png' #Choose a custom icon, see here: https://kml4earth.appspot.com/icons.html
    sharedstyle.iconstyle.scale = 0.25 #Shrink the icons so they look a little nicer
    pts = df[np.logical_and(df[' ChannelWidth_m'] < maxc, df[' ChannelWidth_m'] >= minc)] #Grab all points within the bin range
    for pt in pts.get_values():
        lat, lon = utm.to_latlon(pt[0], pt[1], 44, 'U') #Reproject
        #pnt = fol.newpoint(name="{0},{1}".format(lon, lat), coords=[(lon,lat)]) 
        pnt = fol.newpoint(name='', coords=[(lon,lat)]) #plot each point
        pnt.style = sharedstyle #Set the point to have the above shared style
        
kml.save(base + 'Channel_Width.kml')

#Create a linear feature
kml = Kml()
fol = kml.newfolder(name="Line String")

lins = []
for pt in df.get_values():
    lat, lon = utm.to_latlon(pt[0], pt[1], 44, 'U')
    lins.append([lon,lat])   

lin = kml.newlinestring(name="Channel", description="Channel",coords=lins)

kml.save(base + 'Channel_Line.kml')
#NOTE - THIS DOESNT DO ANY SMOOTHING, it just creates a line from all of the points in the csv. Could do some sort of interpolation on the x/y pairs maybe, or clean it after it is a linestring. 



#Red – ff0000ff.
#Yellow – ff00ffff.
#Blue – ffff0000.
#Green – ff00ff00.
#Purple – ff800080.
#Orange – ff0080ff.
#Brown – ff336699.
#Pink – ffff00ff.
#
#Color and opacity (alpha) values are expressed in hexadecimal notation. The range of values for any one color is 0 to 255 (00 to ff). For alpha, 00 is fully transparent and ff is fully opaque. The order of expression is aabbggrr, where aa=alpha (00 to ff); bb=blue (00 to ff); gg=green (00 to ff); rr=red (00 to ff). For example, if you want to apply a blue color with 50 percent opacity to an overlay, you would specify the following: <color>7fff0000</color>, where alpha=0x7f, blue=0xff, green=0x00, and red=0x00.