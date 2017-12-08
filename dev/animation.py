from visual import *
import pandas as pd,numpy as np
from sys import exit

###############################################################
#CONSTANTS
###############################################################
BASE_DIR="../"

RAGC=(17+45./60+40.04/3600)*15
DECGC=-(29+0./60+28.1/3600)
ETAGC=58.5986320306

###############################################################
#ROUTINES AND CLASSES
###############################################################
def sky2dir(RA,DEC,stride=90):
    x=np.cos(radians(DEC))*np.cos(radians(stride+RA))
    z=np.cos(radians(DEC))*np.sin(radians(stride+RA))
    y=np.sin(radians(DEC))
    return vector(x,y,z)

class star(object):
    #rlayers=[1,2,4,6]
    rlayers=[1,1,1,1]    
    olayers=[1,2,5,10]

    def __init__(self,pos=(0,0,0),radius=0.1,opacity=0.5,color=color.yellow):
        self.spheres=[]
        for fr,fo in zip(self.rlayers,self.olayers):
            self.spheres+=[sphere(pos=pos,radius=fr*radius,opacity=opacity/fo,
                                  color=color,material=materials.emissive)]
    def setp(self,prop,value):
        for s in self.spheres:exec("s.%s=%s"%(prop,value))
    def set_pos(self,pos):
        for s in self.spheres:s.pos=pos
    def set_radius(self,radius):
        for i,s in enumerate(self.spheres):
            s.radius=self.rlayers[i]*radius
    def set_opacity(self,opacity):
        for i,s in enumerate(self.spheres):
            s.opacity=opacity/self.olayers[i]

###############################################################
#SCENE PROPERTIES
###############################################################
scale=300.0 #PC
fov=2*scale #PC
rstar=scale/200 #PC
qaxis=True #False

#"""
###############################################################
#SCENE AND CAMMERA
###############################################################
sc=1/fov;
scene = display(title='Interstellar Simulation',
                x=0,y=0,width=1080, height=720,
                center=(0,0,0),background=(0,0,0),
                forward=(+1,-0.2,-0.2),scale=(sc,sc,sc),
                ambient=color.gray(0.5),
                lights= []
                )

###############################################################
#AXIS
###############################################################
axradius=fov/500
if qaxis:
    xaxis=cylinder(pos=(0,0,0),axis=(scale,0,0),radius=axradius)
    xaxis=cylinder(pos=(0,0,0),axis=(-scale,0,0),radius=axradius)
    label(text='x',align='center',pos=(scale,0,0),box=False,opacity=0)
    yaxis=cylinder(pos=(0,0,0),axis=(0,scale,0),radius=axradius)
    yaxis=cylinder(pos=(0,0,0),axis=(0,-scale,0),radius=axradius)
    label(text='y',align='center',pos=(0,scale,0),box=False,opacity=0)
    zaxis=cylinder(pos=(0,0,0),axis=(0,0,scale),radius=axradius)
    zaxis=cylinder(pos=(0,0,0),axis=(0,0,-scale),radius=axradius)
    label(text='z',align='center',pos=(0,0,scale),box=False,opacity=0)

###############################################################
#MILKYWAY BACKGROUND
###############################################################
f=frame()
mradius=100*scale
mopacity=1.0

data=materials.loadTGA("photo-milkyway")
tgrid=materials.texture(data=data,mapping="spherical")
sky1=sphere(frame=f,pos=(0,0,0),radius=mradius,opacity=mopacity,material=tgrid)
sky1.rotate(angle=radians(180+RAGC),origin=(0,0,0),axis=(0,1,0))
sky1.rotate(angle=radians(-DECGC),origin=(0,0,0),axis=(0,0,1))
sky1.rotate(angle=radians(ETAGC),origin=(0,0,0),axis=(1,0,0))
#"""


###############################################################
#READ OBJECTS
###############################################################
info=pd.read_csv(BASE_DIR+"simstars-Oumuamua.csv")
data=pd.read_csv(BASE_DIR+"simulation-Oumuamua.csv")
ntimes=len(data)
nobjs=int((data.shape[1]-2)/6)
print(data.shape)
print("Ntimes = ",ntimes)
print("Nobjs = ",nobjs)
ts=data["t"].values
robjs=np.zeros((nobjs,ntimes,3))
for i in range(nobjs):
    if i==2:continue
    pref="x%d"%i
    robjs[i]=np.vstack((data["%s_0"%pref].values,
                        data["%s_1"%pref].values,
                        data["%s_2"%pref].values)).transpose()
rsuns=robjs[0]
rnoms=robjs[1]

###############################################################
#CREATE OBJECTS
###############################################################
#GET THE RANGES
sun=star(pos=(0,0,0),radius=rstar)
tlabel=label(por=(0,0,0),text="t=0",yoffset=-5*rstar)


wanderer=sphere(pos=robjs[1,0,:],radius=2*rstar,color=color.red,
                make_trail=True,retain=500)

maxnum=200
stars=[]
labels=[]
slabels=[]
dmins=[]
k=0
for i in range(3,nobjs):
    #print("Initial position:",robjs[i,0,:])
    stars+=[star(pos=robjs[i,0,:],radius=rstar,color=color.white)]
    labels+=[label(pos=robjs[i,0,:],text="",xoffset=-5*rstar,line=0,box=False,opacity=0.0)]
    
    #Get minimum distances
    dmin=min([np.linalg.norm(robjs[i,it,:]-robjs[1,it,:]) for it in range(len(ts))])
    #print("Minimum distance star %d = "%i,dmin)
    dmins+=[dmin]
        
    #Stellar info
    infostar=info.iloc[k]
    s1=str(infostar.hip)
    s2=str(infostar.tycho2_id)
    s3=str(infostar.name_simbad)
    if s1=="nan":s1=""
    if s2=="nan":s2=""
    if s3=="nan":s3=""
    slabel=s1+" "+s2+" "+s3
    slabels+=[slabel]
    k+=1
    
    if i>maxnum:break
    #break

#"""
###############################################################
#ANIMATE
###############################################################
twait=0.05 #/maxnum
trate=50
sct=1/fov
for it,t in enumerate(ts):
    rate(trate)
    tlabel.text="t=%.2f"%(t/1e6)

    #print("t(%d) = "%it,t)
    pwanderer=robjs[1,it,:]
    wanderer.pos=pwanderer
    scene.center=pwanderer
    #scene.range=(100,100,100)

    k=0
    for i in range(3,nobjs):
        #print("pos=",robjs[i,it,:])

        #Position of the star
        pstar=robjs[i,it,:]

        #Distance star-wanderer
        dstar=np.linalg.norm(pstar-pwanderer)

        #Update star position
        stars[k].set_pos(pstar)

        #Set label
        #if dstar<2*dmin:
        #if dstar<1.1*dmins[k] and dmins[k]<5:
        #if dstar==dmins[k] and dmins[k]<5:
        if ("200325" in slabels[k]) and dstar==dmins[k]:
            labels[k].opacity=1.0
            labels[k].box=True
            labels[k].text=slabels[k]+", d=%.1f"%dmins[k]
            labels[k].pos=pstar
            sleep(10.0)
        else:
            labels[k].box=False
            labels[k].opacity=0.0
            labels[k].text=""
 
        k+=1
        if i>maxnum:break
        #break
    #break
#"""
