from visual import *
import pandas as pd,numpy as np

###############################################################
#CONSTANTS
###############################################################
RAGC=(17+45./60+40.04/3600)*15
DECGC=-(29+0./60+28.1/3600)
ETAGC=58.5986320306

###############################################################
#ROUTINES
###############################################################
def sky2dir(RA,DEC,stride=90):
    x=np.cos(radians(DEC))*np.cos(radians(stride+RA))
    z=np.cos(radians(DEC))*np.sin(radians(stride+RA))
    y=np.sin(radians(DEC))
    return vector(x,y,z)

###############################################################
#SCENE AND CAMMERA
###############################################################
scene = display(title='Interstellar Simulation',
                x=0,y=0,width=1080, height=720,
                center=(0,0,0),background=(0,0,0),
                forward=(+1,0,0),scale=(0.5,0.5,0.5),
                ambient=color.gray(0.5),
                lights= []
                )

###############################################################
#AXIS
###############################################################
xaxis=cylinder(pos=(0,0,0),axis=(1,0,0),radius=0.005)
xaxis=cylinder(pos=(0,0,0),axis=(-1,0,0),radius=0.005)
label(text='x',align='center',pos=(1,0,0),box=False,opacity=0)
yaxis=cylinder(pos=(0,0,0),axis=(0,1,0),radius=0.005)
yaxis=cylinder(pos=(0,0,0),axis=(0,-1,0),radius=0.005)
label(text='y',align='center',pos=(0,1,0),box=False,opacity=0)
zaxis=cylinder(pos=(0,0,0),axis=(0,0,1),radius=0.005)
zaxis=cylinder(pos=(0,0,0),axis=(0,0,-1),radius=0.005)
label(text='z',align='center',pos=(0,0,1),box=False,opacity=0)

###############################################################
#MILKYWAY BACKGROUND
###############################################################
f=frame()
mradius=10.0
mopacity=1.0

data=materials.loadTGA("photo-milkyway")
tgrid=materials.texture(data=data,mapping="spherical")
sky1=sphere(frame=f,pos=(0,0,0),radius=mradius,opacity=mopacity,material=tgrid)
sky1.rotate(angle=radians(180+RAGC),origin=(0,0,0),axis=(0,1,0))
sky1.rotate(angle=radians(-DECGC),origin=(0,0,0),axis=(0,0,1))
sky1.rotate(angle=radians(ETAGC),origin=(0,0,0),axis=(1,0,0))

##circle1=shapes.circle(pos=(0,0,0),radius=1.1*mradius)
##circle2=shapes.circle(pos=(0,0,0),radius=1.2*mradius)
##straight=[(0.0,-mradius/5,0.0),(0.0,+mradius/5,0.0)]
##extrusion(pos=straight,shape=circle2-circle1,color=color.black)

data=materials.loadTGA("tycho-milkyway")
tgrid=materials.texture(data=data,mapping="spherical")
sky2=sphere(frame=f,pos=(0,0,0),radius=1.3*mradius,opacity=1.0,material=tgrid)
sky2.rotate(angle=radians(90+RAGC),origin=(0,0,0),axis=(0,1,0))
sky2.rotate(angle=radians(-DECGC),origin=(0,0,0),axis=(0,0,1))
sky2.rotate(angle=radians(ETAGC),origin=(0,0,0),axis=(1,0,0))

###############################################################
#OBJECTS
###############################################################
sphere(pos=(0,0,0),radius=0.1,material=materials.emissive,opacity=0.5)
lamp = local_light(pos=(0,0,0), color=color.red)


