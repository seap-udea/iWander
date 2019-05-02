#########################################################
#    _ _       __                __         		#
#   (_) |     / /___ _____  ____/ /__  _____		#
#  / /| | /| / / __ `/ __ \/ __  / _ \/ ___/		#
# / / | |/ |/ / /_/ / / / / /_/ /  __/ /    		#
#/_/  |__/|__/\__,_/_/ /_/\__,_/\___/_/     		#
# Dynamics of Interestellar Wanderers			#
# Jorge I. Zuluaga et al. [)] 2017			#
# http://github.com/seap-udea/iWander.git		#
#########################################################
# Python routines for iWander
#########################################################

############################################################
#EXTERNAL MODULES
############################################################
import sys,glob,re,os
import numpy as np, scipy as sp, pandas as pd,time
from matplotlib import use
use('Agg')
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats as st,signal
from collections import OrderedDict


############################################################
#MACROS
############################################################
gaussian=st.norm
maxwell=st.maxwell
exit=sys.exit
argv=sys.argv

############################################################
#CONSTANTS AND MACROS
############################################################

#BASEDIR 
BD="."

from matplotlib.ticker import FormatStrFormatter

#CONSTANTS
AU=1.496e8 #km

DEG=np.pi/180
RAD=1/DEG

YEAR=365.25*86400
PARSEC=1/(1/3600.0*np.pi/180)*AU

DEG=np.pi/180
RAD=1/DEG

FIGDIR="figures/"

############################################################
#ROUTINES
############################################################
#ROUTINES
def figure(qax=True,**figopts):
    #default=dict(qax=True)
    fig=plt.figure(**figopts)
    if qax:ax=fig.add_subplot(111)
    else:ax=None
    return fig,ax

def subPlots(panels,ncols=1,l=0.1,b=0.1,w=0.8,figsize=(8,6),dh=None,dw=None,fac=2.0):
    """
    Subplots
    """
    import numpy
    npanels=len(panels)
    spanels=sum(panels)

    # GET SIZE OF PANELS
    b=b/npanels
    if dh is None:dh=[b/2]*npanels
    elif type(dh) is not list:dh=[dh]*npanels
    else:
        dh+=[0]

    if dw is None:dw=w/5
        
    # EFFECTIVE PLOTTING REGION
    hall=(1-fac*b-sum(dh))
    hs=(hall*numpy.array(panels))/spanels
    fach=(1.0*max(panels))/spanels

    # SEE IF THERE IS MORE THAN ONE COLUMN
    wc=w/2 if ncols>1 else w
    
    # CREATE AXES
    fig=plt.figure(figsize=(figsize[0],figsize[1]/fach))
    axs=[]
    
    bo=b
    axc=[]
    for i in range(npanels):
        ax=fig.add_axes([l,b,wc,hs[i]])
        axc+=[ax]
        b+=hs[i]+dh[i]
    axs+=[axc]
        
    b=bo
    if ncols>1:
        for j in range(ncols-1):
            b=bo
            l+=wc+dw
            axc=[]
            for i in range(npanels):
                ax=fig.add_axes([l,b,wc,hs[i]])
                axc+=[ax]
                b+=hs[i]+dh[i]
            axs+=[axc]
            
    return fig,axs

def uniform(a,b):
    return a+np.random.rand()*(b-a)

#Maxwellian
def maxwellDistrib(x,mu=1):
    f=np.sqrt(2/np.pi)
    a=mu/(2*f)
    p=f*x**2*np.exp(-x**2/(2*a**2))/a**3
    return p

def wFunction(d,h):
    """                                                                                                                                                                              
    Schoenber B-spline function                                                                                                                                                      
    See: https://arxiv.org/pdf/1012.1885.pdf                                                                                                                                         
                                                                                                                                                                                     
    Plot:                                                                                                                                                                            
        h=0.1                                                                                                                                                                        
        sigma=wNormalization(h)                                                                                                                                                      
        fig=plt.figure()                                                                                                                                                             
        ax=fig.gca()                                                                                                                                                                 
        ds=np.linspace(0,5*h,100)                                                                                                                                                    
        ws=np.array([sigma*wFunction(d,h) for d in ds])                                                                                                                              
        ax.plot(ds,ws)                                                                                                                                                               
        fig.savefig("scratch/weighting-shoenberg.png")                                                                                                                               
                                                                                                                                                                                     
    Test it:                                                                                                                                                                         
        from scipy.integrate import quad                                                                                                                                             
        wnorm=lambda d:wFunction(d,h)*sigma                                                                                                                                          
        print quad(wnorm,0,2*h)                                                                                                                                                      
    """
    q=d/h
    if q<1:w=0.25*(2-q)**3-(1-q)**3
    elif q<2:w=0.25*(2-q)**3
    else:w=0
    return w

def wNormalization(h):
    from scipy.integrate import quad
    sigma=1/quad(wFunction,0,2*h,args=(h,))[0]
    return sigma

def wFunction2(d,h):
    q=d/h
    w=1/h**3*np.exp(-d**2/h**2)
    return w

def wNormalization2(h):
    from scipy.integrate import quad
    sigma=1/quad(wFunction2,0,2*h,args=(h,))[0]
    return sigma

import re
RE_CHILD_ARRAY = re.compile(r'(.*)\[(.*)\]')
RE_INTERNAL_ATTR = re.compile('__.*__')
def child_attrs_of(klass):
    """
    Given a Node class, get a set of child attrs.
    Memoized to avoid highly repetitive string manipulation
    """
    non_child_attrs = set(klass.attr_names)
    all_attrs = set([i for i in klass.__slots__ if not RE_INTERNAL_ATTR.match(i)])
    return all_attrs - non_child_attrs

def parse2dict(node):
    """ Recursively convert an ast into dict representation. """


    klass = node.__class__

    result = {}

    # Metadata
    result['_nodetype'] = klass.__name__

    # Local node attributes
    for attr in klass.attr_names:
        result[attr] = getattr(node, attr)

    # Coord object
    if node.coord:
        result['coord'] = str(node.coord)
    else:
        result['coord'] = None

    # Child attributes
    for child_name, child in node.children():
        # Child strings are either simple (e.g. 'value') or arrays (e.g. 'block_items[1]')
        match = RE_CHILD_ARRAY.match(child_name)
        if match:
            array_name, array_index = match.groups()
            array_index = int(array_index)
            # arrays come in order, so we verify and append.
            result[array_name] = result.get(array_name, [])
            if array_index != len(result[array_name]):
                raise CJsonError('Internal ast error. Array {} out of order. '
                    'Expected index {}, got {}'.format(
                    array_name, len(result[array_name]), array_index))
            result[array_name].append(parse2dict(child))
        else:
            result[child_name] = parse2dict(child)

    # Any child attributes that were missing need "None" values in the json.
    for child_attr in child_attrs_of(klass):
        if child_attr not in result:
            result[child_attr] = None

    return result

verbose=0
def readValue(val):
    if verbose:print("Analysing:",val)
    result=""

    try:val["_nodetype"]
    except:return result

    if val["_nodetype"]=="Constant":
        result=eval("%s"%val["value"])
        
    if val["_nodetype"]=="ID":
        result=val["name"]

    if val["_nodetype"]=="UnaryOp":
        pass

    if val["_nodetype"]=="BinaryOp":
        op=val["op"]
        term1=readValue(val["left"]);
        term2=readValue(val["right"]);
        result=term1+op+term2

    if val["_nodetype"]=="UnaryOp":
        op=val["op"]
        term1=val["expr"]["value"];
        result=op+term1

    if verbose:print("Result:",result)
    return result

def readConf(filename):
    from pycparser import parse_file,c_generator
    os.system("grep -v '^/' %s > .ccode"%filename)
    confori=parse2dict(parse_file(".ccode"))
    conf=dict()
    for val in confori["ext"]:
        name=val["name"]
        if verbose:print("Reading node:",name)
        if verbose:print("String:",val["init"])
        conf[val["name"]]=readValue(val["init"])
        if verbose:print("Value:",conf[val["name"]])
        if verbose:print("*"*50)
    return conf

#Intelligent Shell script execution
def _run(cmd):
    import sys,subprocess
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    while True:
        line = p.stdout.readline().rstrip()
        if not line:
            break
        yield line
    (output,error)=p.communicate()
    yield p.returncode,error

def System(cmd,verbose=True):
    out=[]
    for path in _run(cmd):
        try:
            if verbose:print(path.decode("utf-8"))
            out+=[path.decode("utf-8")]
        except:
            out+=[(path[0],path[1].decode("utf-8"))]
            pass
    return out
