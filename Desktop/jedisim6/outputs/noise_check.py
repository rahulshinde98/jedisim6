## Sam Ferraro 7/3/23 13:34:00
## Analyze how which background galaxies are left after analysis

# Import Libraries
import pandas as pd
import numpy as np
import sys
#import plotille

#Define how to make a histogram
#def hist(mags):
  #Graph a histogram
  #fig = plotille.Figure()
  #fig.width = 60
  #fig.height = 50
  #fig.set_x_limits(min_=np.min(mags), max_=np.max(mags))
  #fig.histogram(mags, bins=20)
  #print(fig.show(legend=False))

# Define a search algorthim that finds the nearest filled pixel in a 2d array
def ringsearch(grid, x, y, depth=0):
  #3 for loops make a square around the point *depth* pixels away from the center
  for dx in range(depth +1):
    for i in [-1,1]:
      for j in [-1,1]:
        #For the horizontal bars
        #Make sure that the the space being checked exists
        if x +i*dx >= 0 and x +i*dx<len(grid) and y +j*depth >= 0 and y +j*depth <len(grid):
          #Return the value up if it is not -1
          val = grid[x +i*dx][y +j*depth]
          if val != -1:
            return val
        #Repeat for the vertical bars
        if y +i*dx >= 0 and y +i*dx<len(grid) and x +j*depth >= 0 and x +j*depth <len(grid):
          val = grid[x +j*depth][y +i*dx]
          if val != -1:
            return val
  return ringsearch(grid, x, y, depth+1)

# Define Processing
def process(bg, width, trim):
  # Get files
  print(bg + ":")
  preloc = bg + '/' + bg + '_catalog.txt'
  # Pre-processing catalog
  pre = pd.read_csv(preloc, sep = '\t', 
                    names=['file','x','y','deg','z','pixarc','t1','t2',
                           'mag','t3','stamp','dist'])
  postloc = bg + '/final_' + bg + '.csv'
  # Post-processing catalog
  post = pd.read_csv(postloc)
  
  # Resize pre catalog
  conv = 0.06/0.263
  pre['x'] = pre['x'].apply(lambda x: (x -trim)*conv)
  pre['y'] = pre['y'].apply(lambda x: (x -trim)*conv)
  
  # Set up pixelated background. I did this to search because it is more efficient.
  resolution = 0.1 #pixels
  #width = 41913 #pixels
  width -= 2*trim
  width *= conv
  width = round(width) +2
  print("Width: " + str(width))
  size = int(width/resolution)
  grid = [[-1 for k in range(size)] for i in range(size)]
  
  loss = 0
  # Fill the grid in with magnitudes
  for index, row in pre.iterrows():
    if int(round(row['x']/resolution)) < 0 or int(round(row['y']/resolution)) < 0:
      continue
    else:
      # Check how many overlap/error handling
      try:
        if grid[int(round(row['x']/resolution))][int(round(row['y']/resolution))] != -1 \
           and grid[int(round(row['x']/resolution))][int(round(row['y']/resolution))] != row['mag']:
          loss += 1
      except:
        print(int(round(row['x']/resolution)),int(round(row['y']/resolution)))
        print(len(grid))
        print(len(grid[1]))
      # Fill grid
      grid[int(round(row['x']/resolution))][int(round(row['y']/resolution))] = row['mag']
  print("Loss: " + str(loss))
  
  #Find the magnitude each row is most likely to be
  mags = []
  for index, row in post.iterrows():
    x = int(round(row['x']/resolution))
    y = int(round(row['y']/resolution))
    mags.append(ringsearch(grid, x, y))
  
  #Analyze the magnitudes
  divs = 12
  cats = [[0 for i in range(divs)] for j in range(2)]
  low = int(np.min(mags))
  high = int(np.max(mags)) +1
  diff = high - low
  
  #Find out how many of each bracket exist in pre
  for index, row in pre.iterrows():
    i = row['mag']
    i -= low
    i = int(i /diff *divs)
    try:
      cats[1][i] += 1
    except:
      print(i)
  
  #Find how many of each bracket exist in post
  for i in mags:
    i -= low
    i = int(i /diff *divs)
    cats[0][i] += 1
  
  #Print out the info
  print("-------------------------")
  print("x<" + str(low+diff/float(divs)) + "      | " + str(cats[0][0]) + "/" + str(cats[1][0]))
  for i in range(1,divs -1):
    print(str(low + diff/float(divs)*i) + "<x<" + str(low + diff/float(divs)*(i+1)) + "   | " + str(cats[0][i]) + "/" + str(cats[1][i]))
  print(str(low + diff/float(divs)*(divs-1)) + "<x      | " + str(cats[0][-1]) + "/" + str(cats[1][-1]))
  print("-------------------------")
  num = (25 - low) // (diff/float(divs))
  s = 0
  for i in range(int(num)):
    s += cats[0][i]
  print("Total <25: " + str(s))
  s = 0
  s1 = 0
  for i in range(int(divs - num)+1):
    s += cats[0][-i]
    s1 += cats[1][-i]
  print(">25: " + str(s) + "/" + str(s1))
  print("")
  
  return mags


#Actual Script

with open(sys.argv[2], 'r') as f:
  for line in f:
    lsep = line.split("=")
    if lsep[0] == 'nx':
      width = int(lsep[1].split(' ')[0])
    if lsep[0] == 'x_trim':
      trim = int(lsep[1].split('\t')[0])
      break

mags = process(sys.argv[1], width, trim)
