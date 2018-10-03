
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as CM
import matplotlib.colors as COL
from matplotlib import ticker


# 5 continuous variables - x,y, color, area, symbol

areaRange = [100,360]  # sets the dynamic range for the plot symbol areas, increase to make larger symbols
rotRange = [0.,0.5]    # sets the range of phase for the notch rotations, ex. [0.5,1] would give notches below points
notchScale = 2.0       # sets the length of the notch relative to the radius of the plot circle, 2 is for twice as long, 3 for three times, etc.


# general plotting choices formatting choices
labelSize1 = 7
fontSize1 = 9
legendSize1 = 5

widthST = 3.45 #inches
heightST = 2.5875  # fig.set_figwidth assumes inches
Rect = (-0.04,-0.05,1,1) # makes the plot take up more of the plotting space



def Deconfusogram(xin,yin,colin,areain,rotin,
              cmap='Spectral_r',
              ax=None,
              xlims=None,
              ylims=None,
              xlabel="Period (days)",
              ylabel="Mass ($M_{\odot}$)",
              collabel="Rotation: axisymmetry       Color: poloidal",
              arealabel="Size: $<B^{2}>$ (kG$^{2}$)",
              file=None,
              AbsCol=True,   # set color and symbol rotation on 0-1 range, assumes inputs are in (0,1)
              AbsRot=True):

    if ax is None:
        fig, axis = plt.subplots()
    else:
        axis = ax
        fig = axis.get_figure()


    if AbsCol:
        colout = CM.ScalarMappable(norm=COL.Normalize(vmin=np.minimum(0,min(colin)),
                                                  vmax=np.maximum(max(colin),1)),
                               cmap=plt.get_cmap(cmap)).to_rgba
    else:
        colout = CM.ScalarMappable(norm=COL.Normalize(vmin=min(colin),
                                              vmax=max(colin)),
                           cmap=plt.get_cmap(cmap)).to_rgba
    
    areaout = COL.Normalize(vmin=min(areain),vmax=max(areain))(areain).data*(areaRange[1]-areaRange[0]) + areaRange[0]

    if AbsRot:
        rotout = COL.Normalize(vmin=np.minimum(0,min(rotin)),vmax=np.maximum(max(rotin),1))(rotin).data * (0.5)
    else:
        rotout = COL.Normalize(vmin=min(rotin),vmax=max(rotin))(rotin).data * (0.5) 

    num = len(xin)

    fig.set_figwidth(widthST*1.1)  #make a little wider to accomadate legend axis
    fig.set_figheight(heightST)
    axis.tick_params(which='both',labelsize=labelSize1)
    
    axis.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    axis.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    for i in range(num):
        theta = np.linspace(0,2*np.pi,1000)
        radius = 50   # arbitrary scale for plot symbol
        verts = np.zeros((len(theta),2))
        verts = np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])

        phase = rotout[i]
        shift = np.int(phase*len(theta))
        # verts here creates the fletched circle symbol
        verts[shift,:] = [radius*notchScale * np.cos(2*np.pi*phase),radius*notchScale*np.sin(2*np.pi*phase)]

        scaleA = np.maximum(np.cos(2*np.pi*phase)**2,np.sin(2*np.pi*phase)**2) # scale needed to make sure area of circles correponds to data 
        
        axis.scatter(xin[i],yin[i],c=colout(colin[i]),edgecolors="black",
                     s=areaout[i]*scaleA,alpha=0.6,verts=verts,linewidths=0.5)

    axis.set_ylabel(ylabel,size=fontSize1)
    axis.set_xlabel(xlabel,size=fontSize1)

    axis.set_position([0.12,0.15,0.55,0.8])  #defines the locaiton of the main data plot

    colAx = fig.add_axes([0.78,0.15,0.1,0.8]) #defines the locaiton of the colorbar
    legAx = colAx.twinx()

    #legAx.set_ylim([min(areain),max(areain)])

    if AbsCol:
        colAx.set_ylim([-0.1,1.1])
    else:
        colAx.set_ylim([min(colin),max(colin)])

    legAx.set_ylim([-0.1,1.1])
    colAx.set_xlim([0.5,1.5])
    colAx.set_xticklabels([])
    colAx.set_xticks([])

    colAx.tick_params(labelsize=legendSize1,length=2)
    legAx.tick_params(labelsize=legendSize1,length=2)

    legAx.tick_params(which='minor',labelsize=legendSize1)
    colAx.tick_params(which='minor',labelsize=legendSize1)

    legAx.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    colAx.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    legAx.set_ylabel(arealabel,size=legendSize1,rotation=-90.,labelpad=10)
    colAx.set_ylabel(collabel,size=legendSize1,labelpad=3)

    legendx = np.ones(6)
    legendy = np.linspace(0,1,6)
    areasleg = legendy*(areaRange[1]-areaRange[0]) + areaRange[0]
    rotleg = legendy*(rotRange[1]-rotRange[0]) + rotRange[0]

    legy = legAx.get_yticks()
    #pdb.set_trace()

    #legAx.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    areaticklab = legy*(max(areain) - min(areain)) + areain[0]
    legAx.set_yticklabels(["{0:5.3G}".format(ll) for ll in areaticklab])

    numleg = len(legendx)

    # populate the colorbar legend
    for i in range(numleg):
        theta = np.linspace(0,2*np.pi,1000)
        radius = 50   # seems arbitrary?
        verts = np.zeros((len(theta),2))
        verts = np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])

        phase = rotleg[i]
        shift = np.int(phase*len(theta))
        verts[shift,:] = [radius*notchScale* np.cos(2*np.pi*phase),radius*notchScale*np.sin(2*np.pi*phase)]
        scaleA = np.maximum(np.cos(2*np.pi*phase)**2,np.sin(2*np.pi*phase)**2)
        colAx.scatter(legendx[i],legendy[i],c=colout(legendy[i]),
                         edgecolors="black",s=areasleg[i]*scaleA,alpha=0.6,verts=verts,linewidths=0.5)

    if xlims is not None:
        axis.set_xlim(xlims)

    if ylims is not None:
        axis.set_ylim(ylims)

    #fig.set_tight_layout({'rect':[0,0,1,1]})

    if file is not None:
        fig.savefig(file)

    return fig,axis





