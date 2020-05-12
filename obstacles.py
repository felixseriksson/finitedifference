import numpy as np

# shape function return the set of coordinates defining their edge/border,
# not necessarily in order
# (mainly for readibility)

'''
def square(dx, dy, x1=2, y1=1.3, x2=3.4, y2=2.7):
    # default: hörn i (2, 1.3), (3.4, 1.3), (2, 2.7) och (3.4, 2.7)
    # (anpassat för medium i 4 x 16 ruta)
    coords = []
    for xval in np.arange(int(x1), int(x2), dx):
        coords.append((xval, y1))
        coords.append((xval, y2))
    for yval in np.arange(int(y1), int(y2), dy):
        coords.append((x1, yval))
        coords.append((x2, yval))
    return coords

print(square(0.1, 0.1, 0,0,1,1))
'''

'''
obstacledict = {
            None: lambda *args: None, 
            #"smallsquare": smallsquare,
            "mediumsquare": mediumsquare #,
            # "largesquare": largesquare,
            # "smallcylinder": smallcylinder,
            # "mediumcylinder": mediumcylinder,
            # "largecylinder": mediumcylinder,
            # "airfoil": airfoil,
            # "circlewithtail": circlewithtail,
            # "halfwall": halfwall
            
            }
'''