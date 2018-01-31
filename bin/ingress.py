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
# Generate a LaTeX table with the ingress parameters
#########################################################
#!/usr/bin/env python
############################################################
#PACKAGES
############################################################
from iwander import *

############################################################
#READ GENERAL CONFIGURATION FILE
############################################################
conf=readConf("iwander.conf");

############################################################
#READ INGRESS FILE
############################################################
exec(open("scratch/ingress-%s.dat"%conf["WANDERER"]).read());

############################################################
#READ WANDERERS FILE
############################################################
wanderer=pd.read_csv("scratch/wanderer-%s.csv"%conf["WANDERER"])

RAm=wanderer.RA.mean()
dRA=wanderer.RA.std()

DECm=wanderer.DEC.mean()
dDEC=wanderer.DEC.std()

lm=wanderer.l.mean()
dl=wanderer.l.std()

bm=wanderer.b.mean()
db=wanderer.b.std()

Um=wanderer.vxgal.mean()
dU=wanderer.vxgal.std()

Vm=wanderer.vygal.mean()
dV=wanderer.vygal.std()

Wm=wanderer.vzgal.mean()
dW=wanderer.vzgal.std()

############################################################
#CREATE LATEX TABLE
############################################################
mu_nom/=1e11
mu_asy/=1e11
parts=date_asy.split()
dms=parts[3].split(":")
day=(float(dms[0])+float(dms[1])/60.0+float(dms[2])/3600.)/24.
date_asy="%s %.2f %s"%(parts[1],float(parts[2])+day,parts[0])
fac=1e-6

table="""
%%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
%%TABLE: ASYMPTOTIC ELEMENTS
\\begin{table}
  \\centering
  \\scriptsize
  \\begin{tabular}{ll}
  \\hline Property & Value \\\\\\hline\\hline
  Reference Epoch & %.1f TDB = %s\\\\
  Nominal elements (%s)   
                      & $q = %.9f$ AU \\\\
                      & $e = %.9f$ \\\\
                      & $i = %.6f$ deg \\\\
                      & $\\Omega = %.6f$ deg \\\\
                      & $\\omega = %.6f$ deg \\\\
                      & $M = %.6f$ deg \\\\
                      & $\\mu = %.9f\\sci{11}$ km$^3$/s$^2$\\\\
  Epoch of Asymptotic elements & %.2f TDB = %s\\\\
  Asymptotic elements & $q = %.9f$ AU \\\\
                      & $e = %.9f$ \\\\
                      & $i = %.6f$ deg \\\\
                      & $\\Omega = %.6f$ deg \\\\
                      & $\\omega = %.6f$ deg \\\\
                      & $M = %.6f$ deg \\\\
                      & $\\mu = %.9f\\sci{11}$ km$^3$/s$^2$\\\\
  Asymptotic covariance 
                      & Eccentricity \\\\
  $\\times 10^{%d}$   & $C_{ee}$ = %.3f\\\\
                      & $C_{eq}$ = %.3f\\\\
                      & $C_{et}$ = %.3f\\\\
                      & $C_{e\\Omega}$ = %.3f\\\\
                      & $C_{e\\omega}$ = %.3f\\\\
                      & $C_{ei}$ = %.3f\\\\

                      & Perihelion distance \\\\
                      & $C_{qq}$ = %.3f\\\\
                      & $C_{et}$ = %.3f\\\\
                      & $C_{q\\Omega}$ = %.3f\\\\
                      & $C_{q\\omega}$ = %.3f\\\\
                      & $C_{qi}$ = %.3f\\\\

                      & Periapsis time \\\\
                      & $C_{tt}$ = %.3f\\\\
                      & $C_{t\\Omega}$ = %.3f\\\\
                      & $C_{t\\omega}$ = %.3f\\\\
                      & $C_{ti}$ = %.3f\\\\

                      & Long. ascending node \\\\
                      & $C_{\\Omega\\Omega}$ = %.3f\\\\
                      & $C_{\\Omega\\omega}$ = %.3f\\\\
                      & $C_{\\Omega i}$ = %.3f\\\\

                      & Argument of periapsis\\\\
                      & $C_{\\omega\\omega}$ = %.3f\\\\
                      & $C_{\\omega i}$ = %.3f\\\\

                      & Inclination\\\\
                      & $C_{ii}$ = %.3f\\\\

  Truncation radius & $%.1f$ AU \\\\
  Time of ingress & $%.1f$ years \\\\
  Radiant at ingress & RA = $%.2f\\pm %.2f$ deg\\\\
                     & DEC = $%.2f\\pm %.2f$ deg \\\\
                     & l =  $%.2f\\pm %.2f$ deg \\\\
                     & b = $%.2f\\pm %.2f$ deg \\\\
  Velocity at ingress & U = $%.3f\\pm %.3f$ km/s\\\\
                      & V = $%.3f\\pm %.3f$ km/s\\\\
                      & W = $%.3f\\pm %.3f$ km/s\\\\\\hline
  \\end{tabular}
\\caption{Properties of the \\%s's orbit in the Solar System.\\label{tab:AsymptoticElements}}
\\end{table}
%%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
"""%(Epoch,Epoch_Date,
     Reference,
     q_nom,e_nom,i_nom,W_nom,w_nom,M_nom,mu_nom,
     t_asy,date_asy,
     q_asy,e_asy,i_asy,W_asy,w_asy,M_asy,mu_asy,
     int(np.log10(fac)),
     cov_ee[0]/fac,
     cov_eq[0]/fac,
     cov_et[0]/fac,
     cov_eW[0]/fac,
     cov_ew[0]/fac,
     cov_ei[0]/fac,
     cov_qq[0]/fac,
     cov_qt[0]/fac,
     cov_qW[0]/fac,
     cov_qw[0]/fac,
     cov_qi[0]/fac,
     cov_tt[0]/fac,
     cov_tW[0]/fac,
     cov_tw[0]/fac,
     cov_ti[0]/fac,
     cov_WW[0]/fac,
     cov_Ww[0]/fac,
     cov_Wi[0]/fac,
     cov_ww[0]/fac,
     cov_wi[0]/fac,
     cov_ii[0]/fac,
     truncation,t_ing,
     RAm,dRA,DECm,dDEC,
     lm,dl,bm,dl,
     Um,dU,Vm,dV,Wm,dW,
     conf["WANDERER"])     

print(table)
