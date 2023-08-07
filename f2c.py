#-----------------------------------------------------------------------------
#f2c.py: fluorescence signals to [Ca+2]
#Description: given Mag-Fluo-4's dF signals, calculate and plot the [Ca2+].
#The units of the calculations are in uM and ms; the temperature is 23oC. 
#-----------------------------------------------------------------------------
#import requiered libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.ticker as ticker
from fontTools.ttLib import TTFont
#-----------------------------------------------------------------------------    
#read in the Mag-Fluo-4 measurement in mouse fibers activated by a single AP    
f2b=pd.read_excel('inputs/f2b.xlsx')
f2x=pd.read_excel('inputs/f2x.xlsx')
f2a=pd.read_excel('inputs/f2a.xlsx')
f1=pd.read_excel('inputs/f1.xlsx')
#convert f dataframes to numpy array
f2b=f2b.to_numpy()
f2x=f2x.to_numpy()
f2a=f2a.to_numpy()
f1=f1.to_numpy()
#select time points 200 to 1490 for the analysis
t2b_start=200;t2b_end=1490
t2x_start=200;t2x_end=1490
t2a_start=200;t2a_end=1490
t1_start=200;t1_end=1490
#the time (t) used to interp base goes from 50 (t_start) to 110 ms (t_end)
#shifted -50 ms to start simulation at 0 ms
t2b=f2b[t2b_start:t2b_end,0]-50
t2x=f2x[t2x_start:t2x_end,0]-50
t2a=f2x[t2a_start:t2a_end,0]-50
t1=f2x[t1_start:t1_end,0]-50
#the time used for simulation (ts) goes from 50 to 110 ms
ts=np.arange(0,60,0.1)
#get the fluorescence values for the selected time base
f2b=f2b[t2b_start:t2b_end,1]
f2x=f2x[t2x_start:t2x_end,1]
f2a=f2a[t2a_start:t2a_end,1]
f1=f1[t1_start:t1_end,1]
#-----------------------------------------------------------------------------
#set minimum fluorescence (A.U.)
Fmin=0.14
#set peak of fluorescence during one AP (A.U.)
pF2b=6.3635
pF2x=6.3635
pF2a=4.9929
pF1=4.0139
#set total Dye concentration (uM)
DT=229.1
#set maximum fluorescence (A.U.)
Fmax=150.9
#set disociation constant (uM2)
Kd=1.652e+05
#-----------------------------------------------------------------------------
#reescale the fluorescence signal  
F2b=(pF2b*((f2b)/np.max(f2b)))+Fmin
F2x=(pF2x*((f2x)/np.max(f2x)))+Fmin
F2a=(pF2a*((f2a)/np.max(f2a)))+Fmin
F1=(pF1*((f1)/np.max(f1)))+Fmin
#get the concentration from the F trace
dCa2b=(Kd/DT)*(((F2b-Fmin)*(Fmax-Fmin))/(2*((Fmax-F2b)**2)))
dCa2x=(Kd/DT)*(((F2x-Fmin)*(Fmax-Fmin))/(2*((Fmax-F2x)**2)))
dCa2a=(Kd/DT)*(((F2a-Fmin)*(Fmax-Fmin))/(2*((Fmax-F2a)**2)))
dCa1=(Kd/DT)*(((F1-Fmin)*(Fmax-Fmin))/(2*((Fmax-F1)**2)))
#-----------------------------------------------------------------------------
#Ca(t) equals resting [Ca2+](=Ca0) plus the interpolated dCa, where dCa is
#the active [Ca] measurement  
def Ca(xs,x,y): 
    return np.interp(xs,x,y)
#-----------------------------------------------------------------------------
#calculate change of [Ca2+] as a funtion of ts
ca2bs=Ca(ts,t2b,dCa2b)
ca2xs=Ca(ts,t2x,dCa2x)
ca2as=Ca(ts,t2a,dCa2a)
ca1s=Ca(ts,t1,dCa1)
#-----------------------------------------------------------------------------
#plotted traces will have a dotted baseline
#get x-values of a baseline for all plots
base_x=np.array([ts[0],ts[-1]])
#get y-values of a baseline for all plots
base_y=np.zeros(2) 
#-----------------------------------------------------------------------------
#set the figure paramenter
arialfont1=TTFont('fonts/arial.ttf') #file arial.tff include charaters used in the
#plots such as \u207a (plus as superincex)
plt.rcParams["font.family"] = "arialfont1"#set arialfont as default font 
#font size (pt)
plt.rcParams["font.size"] = 12
#lines and axes spines width (pt)
wd=1.2 
#length of mayor ticks (pt)
lma=7
#length of minor ticks (pt)
lmi=4
#remove top axes spines 
plt.rcParams['axes.spines.top'] = False
#remove rigth axes spines 
plt.rcParams['axes.spines.right'] = False
#lines with (pt)
plt.rcParams['lines.linewidth'] = wd
#set the value globally 
plt.rcParams['axes.linewidth'] = wd
#length of major x-axis ticks (pt)
plt.rcParams['xtick.major.size'] = lma
#length of minor x-axis ticks (pt)
plt.rcParams['xtick.minor.size'] = lmi
#length of major y-axis ticks (pt)
plt.rcParams['ytick.major.size'] = lma
#length of minor y-axis ticks (pt)
plt.rcParams['ytick.minor.size'] = lmi
#width of major x-axis ticks (pt)
plt.rcParams['xtick.major.width'] = wd
#width of minor x-axis ticks (pt)
plt.rcParams['xtick.minor.width'] = wd
#width of major y-axis ticks(pt)
plt.rcParams['ytick.major.width'] = wd
#width of minor y-axis ticks(pt)
plt.rcParams['ytick.minor.width'] = wd
#-----------------------------------------------------------------------------
#now do the plotting of the Ca2+ measurements during single AP
fig, ax = plt.subplots(figsize=(4,2))
ax.plot(ts,ca2bs,'k')
ax.plot(ts,ca2xs,'r')
ax.plot(ts,ca2as,'b')
ax.plot(ts,ca1s,'g')
ax.set_ylim([-2,19])
ax.set_xlim([0,60])
ax.set_xlabel('Time (ms)')
ax.set_ylabel('\u0394[Ca\u00b2\u207a] (\u03BCM)')
ax.set_xticks([10,30,50])
ax.yaxis.set_minor_locator(ticker.FixedLocator([4.5,13.5]))
ax.xaxis.set_minor_locator(ticker.FixedLocator([20,40]))
ax.set_yticks([0,9,18])
fig.savefig("figures/cae.pdf",bbox_inches='tight',pad_inches=0,transparent=\
True)
#-----------------------------------------------------------------------------
#join the vectors that generate the figure cae.pdf
cae=[ts,ca1s,ca2as,ca2xs,ca2bs]
cae=pd.DataFrame(cae).T
cae.columns =['t (ms)','I \u0394[Ca\u00b2\u207a] (\u03BCM)',\
'IIA \u0394[Ca\u00b2\u207a] (\u03BCM)','IIX \u0394[Ca\u00b2\u207a] (\u03BCM)',\
'IIB \u0394[Ca\u00b2\u207a] (\u03BCM)']
cae.to_excel('outputs/cae.xlsx',index=False)
#-----------------------------------------------------------------------------    
#read in the Mag-Fluo-4 measurement in mouse fibers activated by tetanus    
f1t=pd.read_excel('inputs/f1t.xlsx')
f2bt=pd.read_excel('inputs/f2bt.xlsx')
#convert f dataframes to numpy array
f1t=f1t.to_numpy()
f2bt=f2bt.to_numpy()
#select time points 200 to 5752 for the analysis
t1t_start=200;t1t_end=5752
t2bt_start=200;t2bt_end=12752
#the time (t) used to interp base goes from 50 (t_start) to 600 ms (t_end)
#shifted -50 ms to start simulation at 0 ms
t1t=f1t[t1t_start:t1t_end,0]-92
t2bt=f2bt[t2bt_start:t2bt_end,0]-40
#get the fluorescence values for the selected time base
f1t=f1t[t1t_start:t1t_end,1]
f2bt=f2bt[t2bt_start:t2bt_end,1]
#-----------------------------------------------------------------------------
#reescale the fluorescence signal  
F1t=(pF1*((f1t)/8.17462))+Fmin
F2bt=(pF2b*((f2bt)/np.max(f2bt)))+Fmin
#get the concentration from the F trace
dCa1t=(Kd/DT)*(((F1t-Fmin)*(Fmax-Fmin))/(2*((Fmax-F1t)**2)))
dCa2bt=(Kd/DT)*(((F2bt-Fmin)*(Fmax-Fmin))/(2*((Fmax-F2bt)**2)))
#the time used for simulation (tts) goes from 50 to 600 ms
tts=np.arange(0,600,0.1)
#-----------------------------------------------------------------------------
#calculate change of [Ca2+] as a funtion of ts
ca1ts=Ca(tts,t1t,dCa1t)
ca2bts=Ca(tts,t2bt,dCa2bt)
#-----------------------------------------------------------------------------
#now do the plotting of the Ca2+ measurements during tetanus
fig,ax=plt.subplots(figsize=(4,2))
ax.plot(tts,ca2bts,'k')
ax.plot(tts,ca1ts,'g')
ax.set_ylim([-2,35])
ax.set_xlim([0,600])
ax.set_ylabel('\u0394[Ca\u00b2\u207a] (\u03BCM)')
ax.set_xlabel('Time (ms) ⁺')
ax.set_yticks([0,16,32])
ax.set_xticks([200,400])
ax.yaxis.set_minor_locator(ticker.FixedLocator([8,24]))
ax.xaxis.set_minor_locator(ticker.FixedLocator([100,300]))
fig.savefig("figures/caet.pdf",bbox_inches='tight',pad_inches=0,transparent=\
True)
#-----------------------------------------------------------------------------
#join the vectors that generate the figure caet.pdf
caet=[tts,ca1ts,ca2bts]
caet=pd.DataFrame(caet).T
caet.columns =['t (ms)','I \u0394[Ca\u00b2\u207a] (\u03BCM)',\
'IIB \u0394[Ca\u00b2\u207a] (\u03BCM)']
caet.to_excel('outputs/caet.xlsx',index=False)
#-----------------------------------------------------------------------------
