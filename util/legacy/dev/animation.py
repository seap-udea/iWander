from visual import *
import pandas as pd,numpy as np
from sys import exit

###############################################################
#CONSTANTS
###############################################################
BASE_DIR="../"

#ASTRONOMICAL CONSTANTS
RAGC=(17+45./60+40.04/3600)*15
DECGC=-(29+0./60+28.1/3600)
ETAGC=58.5986320306

#OBJECT
#nwanderer="Oumuamua";fwanderer="Oumuamua";twanderer="1.8 kyr"
nwanderer="Voyayer 1";fwanderer="Voyager1";twanderer="1.8 kyr"

#SCENE PROPERTIES
scale=200.0 #PC
srange=2*scale #PC
rstar=scale/150.0 #PC
axradius=1
retain=1e100

wcolor=color.yellow

#VIEW PARAMETERS
phi0=radians(165.0)
dphi=-np.pi

theta0=radians(85.0)
thetamax=radians(10.0)
dtheta=np.pi

#BACKGROUND PROPERTIES
mradius=100*scale
mopacity=1.0

#WANDERER TRAJECTORY
options=dict(q5retain=500)

#SHOW NO SHOW
qaxis=1

#NUMBER OF OBJECTS
maxnum=200

#ANIMATION
twait=1.0
trate=100
dmax=5.0

#FONT
font="monospace"

###############################################################
#ROUTINES AND CLASSES
###############################################################
def sky2dir(RA,DEC,stride=90):
    x=np.cos(radians(DEC))*np.cos(radians(stride+RA))
    z=np.cos(radians(DEC))*np.sin(radians(stride+RA))
    y=np.sin(radians(DEC))
    return vector(x,y,z)


def setView(i,n,phi0=0.0,theta0=np.pi/2,dphi=2*np.pi,dtheta=0.0,thetamax=0.0):

    #Equivalen unit time
    t=(1.*i)/n

    #Angles
    theta=theta0+thetamax*np.sin(dtheta*t)
    phi=phi0+dphi*t

    #Coordinates on unit sphere
    y=-np.cos(theta)
    x=-np.sin(theta)*np.cos(phi)
    z=-np.sin(theta)*np.sin(phi)
    return dict(phi=phi,theta=theta,forward=(x,y,z))

def circle(r=1,color=color.white,lw=0.1):
    qs=np.linspace(0,2*np.pi,100)
    path=[(r*np.cos(q),0,r*np.sin(q)) for q in qs]
    c=curve(pos=path,radius=lw)
    return c
    
class star(object):
    #rlayers=[1,2,4,6]
    rlayers=[1,1.2,2.2,2.5,3.0]
    #rlayers=[1,1,1,1]    
    olayers=[1,2,5,10,20]

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
    def set_color(self,color):
        for i,s in enumerate(self.spheres):
            s.color=color
    def set_opacity(self,opacity):
        for i,s in enumerate(self.spheres):
            s.opacity=opacity/self.olayers[i]

###############################################################
#SCENE AND CAMMERA
###############################################################
view=setView(0,1,phi0=phi0,theta0=theta0)
scene = display(title='Interstellar Simulation',
                x=0,y=0,width=1080, height=720,
                center=(0,0,0),
                background=(0,0,0),
                forward=view["forward"],
                range=srange,
                ambient=color.gray(0.5),
                lights=[]
                )
zlabel=label(text="Zuluaga et al. (2017) [arXiv:1711.09397]",
      align='center',
      pos=(0,scene.center[1]-1.2*scale,0),
      box=False,opacity=1)

invisible=[]

tsize=scale/10

title=label(text="%s galactic adventure"%nwanderer,
            align='center',
            pos=(0,scene.center[1]+scale,0),
            box=False,opacity=0,height=50)

###############################################################
#AXIS AND SCALES
###############################################################
if qaxis:
    xaxis1=cylinder(pos=(0,0,0),axis=(scale,0,0),radius=axradius)
    xaxis2=cylinder(pos=(0,0,0),axis=(-scale,0,0),radius=axradius)
    xl=label(text='x',align='center',pos=(scale,0,0),box=False,opacity=0)
    yaxis1=cylinder(pos=(0,0,0),axis=(0,scale,0),radius=axradius)
    yaxis2=cylinder(pos=(0,0,0),axis=(0,-scale,0),radius=axradius)
    zl=label(text='z',align='center',pos=(0,scale,0),box=False,opacity=0)
    zaxis1=cylinder(pos=(0,0,0),axis=(0,0,scale),radius=axradius)
    zaxis2=cylinder(pos=(0,0,0),axis=(0,0,-scale),radius=axradius)
    yl=label(text='y',align='center',pos=(0,0,+scale),box=False,opacity=0)
    invisible+=[xaxis1,xaxis2,yaxis1,yaxis2,zaxis1,zaxis2,xl,yl,zl]
    
    for r in np.arange(10,scale,50):
        c=circle(r=r,lw=scale/500)
        invisible+=[c]
    l=label(pos=(0,0,160),text="%.0f pc"%r,box=False,
            line=0,opacity=0,height=10,yoffset=5*rstar)
    invisible+=[l]
    l=label(pos=(0,0,-160),text="-%.0f pc"%r,box=False,
            line=0,opacity=0,height=10,yoffset=5*rstar)
    invisible+=[l]
    

###############################################################
#MILKYWAY BACKGROUND
###############################################################
data=materials.loadTGA("photo-milkyway")
tgrid=materials.texture(data=data,mapping="spherical")
sky1=sphere(pos=(0,0,0),radius=mradius,opacity=mopacity,material=tgrid)
sky1.rotate(angle=radians(180+RAGC),origin=(0,0,0),axis=(0,1,0))
sky1.rotate(angle=radians(-DECGC),origin=(0,0,0),axis=(0,0,1))
sky1.rotate(angle=radians(ETAGC),origin=(0,0,0),axis=(1,0,0))

###############################################################
#READ OBJECTS
###############################################################
info=pd.read_csv(BASE_DIR+"simstars-%s.csv"%fwanderer)
data=pd.read_csv(BASE_DIR+"simulation-%s.csv"%fwanderer)
ntimes=len(data)
nobjs=int((data.shape[1]-2)/6)
#print("Ntimes = ",ntimes)
#print("Nobjs = ",nobjs)
ts=data["t"].values
robjs=np.zeros((nobjs,ntimes,3))
for i in range(nobjs):
    if i==2:continue
    pref="x%d"%i
    robjs[i]=np.vstack((+data["%s_0"%pref].values,
                        +data["%s_2"%pref].values,
                        +data["%s_1"%pref].values)).transpose()
rsuns=robjs[0]
rnoms=robjs[1]

###############################################################
#CREATE OBJECTS
###############################################################
#GET THE RANGES
sun=star(pos=(0,0,0),radius=2*rstar)
sunlabel=label(pos=(0,0,0),text="The Solar System",yoffset=-10*rstar,
             color=color.yellow,opacity=1,box=False)
invisible+=[sunlabel]

oort=sphere(pos=(0,0,0),radius=0.5,opacity=0.20,color=color.blue)

tlabel=label(text="t=%s"%twanderer,
            align='center',
            pos=(0,scene.center[1]-1*scale,0),
            box=False,opacity=1,height=20)

wanderer=sphere(pos=robjs[1,0,:],radius=0.5*rstar,
                color=wcolor,make_trail=True,**options)
olabel=label(pos=wanderer.pos,text=nwanderer,yoffset=+10*rstar,
             color=color.red,opacity=1,box=False)
invisible+=[olabel]

stars=[]
labels=[]
slabels=[]
dmins=[]
k=0
for i in range(3,nobjs):
    stars+=[star(pos=robjs[i,0,:],radius=rstar,color=color.white)]
    labels+=[label(pos=robjs[i,0,:],text="",height=14,xoffset=-5*rstar,line=0,box=False,opacity=0.0)]
    
    #Get minimum distances
    dmin=min([np.linalg.norm(robjs[i,it,:]-robjs[1,it,:]) for it in range(len(ts))])
    dmins+=[dmin]
        
    #Stellar info
    infostar=info.iloc[k]
    s1=str(infostar.hip)
    s2=str(infostar.tycho2_id)
    s3=str(infostar.name_simbad)
    if s1=="nan":s1=""
    if s2=="nan":s2=""
    if s3=="nan":s3=""
    if s1!="":slabel="HIP "+s2
    if s2!="":slabel="TYC "+s2
    if s3!="":slabel=s3.replace("_"," ")
    #slabel=s1+" "+s2+" "+s3
    slabels+=[slabel]
    k+=1
    
    if i>maxnum:break

###############################################################
#ANIMATE
###############################################################
#Wait for starting
ev=scene.waitfor('keydown')

for inv in invisible:
    inv.visible=False
    pass
tlabel.opacity=0
zlabel.opacity=0

oldphi=view["phi"]
for it,t in enumerate(ts):
    rate(trate)
    tlabel.text="t=%.2f Myr"%(t/1e6)

    #Wanderer
    pwanderer=robjs[1,it,:]
    wanderer.pos=pwanderer

    #Scene transformation
    scene.center=pwanderer
    view=setView(it,ntimes,
                 phi0=phi0,dphi=dphi,
                 theta0=theta0,dtheta=dtheta,thetamax=thetamax)
    scene.forward=view["forward"]

    #Label location
    title.pos=scene.center+vector(0,+1*scale,0)
    zlabel.pos=scene.center+vector(0,-1.2*scale,0)
    tlabel.pos=scene.center+vector(0,-1*scale,0)

    k=0
    for i in range(3,nobjs):
        
        #Position of the star
        pstar=robjs[i,it,:]
        
        #Update star position
        stars[k].set_pos(pstar)

        #Distance star-wanderer
        dstar=np.linalg.norm(pstar-pwanderer)

        #Set label
        #"""
        if dstar==dmins[k] and dmins[k]<dmax:
            labels[k].opacity=1.0
            labels[k].box=True
            labels[k].text=slabels[k]+", d=%.1f pc"%dmins[k]
            labels[k].pos=pstar
            stars[k].set_color(wcolor)
            vicinity=sphere(pos=pstar,radius=dmins[k],color=wcolor,
                            opacity=0.3)
            ev=scene.waitfor('keydown')
            #vicinity.visible=False
            #sleep(1.0)
            labels[k].visible=False
        #"""
        
        k+=1
        if i>maxnum:break
