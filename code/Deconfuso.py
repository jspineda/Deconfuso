
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as CM
import matplotlib.colors as COL
from matplotlib import ticker

#import pdb


# If using this code or derivatives thereof, please cite:


# 5 continuous variables - x,y, color, area, symbol

areaRange = [100,500]  # sets the dynamic range for the plot symbol areas, increase to make larger symbols
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
              collabel="Color: Poloidal Fraction",
              arealabel="Size: $<B^{2}>$ (kG$^{2}$)",
              rotlabel = "Fletch Rotation: \n Perc. Axisymmetric",
              rottextscale= 1.3,
              file=None,
              AbsCol=True,   # set color and symbol rotation on 0-1 range, assumes inputs are in (0,1)
              AbsRot=True,
              Manual=False,
              addColorStrip=True):
    """Plotting routine to display ZDI data or more generally 5 variables on single 2-D plot.
        
        Returns a figure and an axis object, with general keywords to adjusting plot contents, all variable arrays should have same length.
        
        Parameters
        ----------
        xin: array_like
             Array pertaining to the X-axis positions of the output scatter points. Defaults set to interpret this as Rotation Period.
        
        yin: array_like
             Array pertaining to the Y-axis positions of the output scatter points. Defaults set to interpret this as Mass.
            
        colin: array_like
               Array pertaining to variable corresponding to scatter point colors. Ideally data in range [0,1]. Defaults set to interpret this as Poloidal degree of magnetic field map.
        
        areain: array_like
                Array pertaining to variable mapped to area of circles for plotting data points. Plot sizes set by areaRange. Defaults to intepret this as Squared magnetic field.
                
        rotin: array_like
               Array pertaining to variable displayed according to flag orientation on each circle scatter point. Range of flags set by rotRange, with 0 to 0.5 corresponding to top 180 degrees of circle. Ideally data in range [0,1]. Defaults set to interpret this as degree of field axisymmetry.
        
        ax: matplotlib axis object
            Pass existing axis object to plot data in.
        
        cmap: string
              Pass string to indicate name of matplotlib color map to use defaults to 'Spectral_r'; works best with continuos maps.
        
        xlims: None or 2 element list
               Passed to axis.set_xlim
        
        ylims: None or 2 element list
               Passed to axis.set_ylim
        
        xlabel: string
                Passed to axis.set_xlabel
                
        ylabel: string
                Passed to axis.set_ylabel
 
        collabel: string
                Passed to colAx.set_ylabel ; Used label the left side of the colorbar axis
                
        arealabel: string
                Passed to legAx.set_ylabel ; Used label the right side of the colorbar axis
        
        rotlabel: string
                Passed to rotAx.set_title ; Used label the top of the rotation legend
        
        rottextscale: float
                Used to scale spacing of rotation label entries to scatter point in legend
        
        file: string
              String as name of output file for saving figure.
              
        AbsCol: boolean
                Default is True, places colin array to plot on absolute range [0,1] even if max or min do not extend to 0 or 1.
        
        AbsRot: boolean
                Default is True, places rotin array to plot on absolute range [0,1] even if max or min do not extend to 0 or 1.
        
        Manual: boolean
                Default is False, for disabling default plotting placements
                
        addColorStrip: boolean
                If True, adds color strip to colorbar
        
        Returns
        -------
        fig: matplotlib figure object
             figure object containing axis of plot
             
        ax: matplotlib axis object
            axis object pertaining to plot data
        
        colA: matplotlib axis object
            axis object pertaining colorbar/ size bar

        rotA: matplotlib axis object
            axis object pertaining fletch legend
        
            
        Notes
        -----
        Be careful with plotting data for colin and rotin that are not already normalized to range of [0,1], to make sure legend bar accurately reflects the data.
        
        Examples
        --------
        
        
        >>> from astropy.table import Table
        >>> import Deconfuso
        
        >>> dat = Table.read("../zdi_data/MasterTable_ZDI.csv")
        >>> fig, ax, colA, rotA = Deconfuso.Deconfusogram(np.log10(dat['Prot(d)']), dat['Mass(Msun)'], dat['Pol.'], 0.1*dat['<B2>(1e5 G2)'],dat['Axisym'],xlabel="$\log_{10}$ Period (d)",file="../plots/Deconfuso_01_Alt.pdf")
        
        
        """

    if ax is None:
        fig, axis = plt.subplots()
    else:
        axis = ax
        fig = axis.get_figure()


    if AbsCol:
        col0 = CM.ScalarMappable(norm=COL.Normalize(vmin=np.minimum(0,min(colin)),
                                                  vmax=np.maximum(max(colin),1)),
                               cmap=plt.get_cmap(cmap))
        colout = col0.to_rgba
    else:
        col0 = CM.ScalarMappable(norm=COL.Normalize(vmin=min(colin),
                                              vmax=max(colin)),
                           cmap=plt.get_cmap(cmap))
        colout = col0.to_rgba
    
    areaout = COL.Normalize(vmin=min(areain),vmax=max(areain))(areain).data*(areaRange[1]-areaRange[0]) + areaRange[0]

    if AbsRot:
        rotout = COL.Normalize(vmin=np.minimum(0,min(rotin)),vmax=np.maximum(max(rotin),1))(rotin).data * (0.5)
    else:
        rotout = COL.Normalize(vmin=min(rotin),vmax=max(rotin))(rotin).data * (0.5) 

    num = len(xin)

    if not Manual:
        fig.set_figwidth(widthST*1.1)  #make a little wider to accomadate legend axis
        fig.set_figheight(heightST)
        axis.tick_params(which='both',labelsize=labelSize1)
    
        axis.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        axis.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    zordindex = np.argsort(areaout)[::-1]

    for l in range(num):
        i = zordindex[l]
        theta = np.linspace(0,2*np.pi,1000)
        radius = 50   # arbitrary scale for plot symbol
        verts = np.zeros((len(theta),2))
        verts = np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])

        phase = rotout[i]
        shift = np.int(phase*(len(theta)-1))
        # verts here creates the fletched circle symbol
        verts[shift,:] = [radius*notchScale * np.cos(2*np.pi*phase),radius*notchScale*np.sin(2*np.pi*phase)]

        scaleA = np.maximum(np.cos(2*np.pi*phase)**2,np.sin(2*np.pi*phase)**2) # scale needed to make sure area of circles correponds to data 
        
        axis.scatter(xin[i],yin[i],c=[colout(colin[i])],edgecolors="black",
                     s=areaout[i]*scaleA,alpha=0.6,marker=verts,linewidths=0.5)

    axis.set_ylabel(ylabel,size=fontSize1)
    axis.set_xlabel(xlabel,size=fontSize1)

    if not Manual:
        axis.set_position([0.12,0.15,0.55,0.8])  #defines the locaiton of the main data plot


    #colAx = fig.add_axes([0.78,0.15,0.1,0.8],label="Colorbar_{}".format(len(fig.axes))) #defines the location of the colorbar, with arbitrary identifier
    colAx = fig.add_axes([0.78,0.15,0.1,0.65],label="Colorbar_{}".format(len(fig.axes))) #defines the location of the colorbar, with arbitrary identifier
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
    legAx.set_yticklabels(["{0:5.2G}".format(ll) for ll in areaticklab])

    numleg = len(legendx)
    notch0 = 1 #  size colorbar has no fletches, nm
    radius0 = radius # set to 50 above



    for i in range(numleg):
        theta = np.linspace(0,2*np.pi,1000)
        #radius = 50   # seems arbitrary?
        verts = np.zeros((len(theta),2))
        verts = np.column_stack([radius0 * np.cos(theta), radius0 * np.sin(theta)])
        
        phase = np.median(rotleg[i])
        shift = np.int(phase*len(theta))
        verts[shift,:] = [radius0*notch0* np.cos(2*np.pi*phase),radius0*notch0*np.sin(2*np.pi*phase)]
        scaleCB = (notch0 / notchScale)**2 #np.maximum(np.cos(2*np.pi*phase)**2,np.sin(2*np.pi*phase)**2)
        colAx.scatter(legendx[i],legendy[i],c=[colout(legendy[i])],
              edgecolors="black",s=areasleg[i]*scaleCB,marker=verts,linewidths=0.5,zorder=100+i)


    if addColorStrip:
        cxs = np.linspace(0,1,100)
        cXX = np.transpose(np.stack((cxs,cxs)))
        legAx.contourf(cXX,cmap=cmap,levels=len(cxs),zorder=0,alpha=0.6,extent=[0.91,1.09,0,1])
        legAx.set_zorder(0)
        colAx.set_zorder(100)
        colAx.patch.set_visible(False)

    # for rotation legend
    rotAx = fig.add_axes([0.78,0.15+0.65,0.1,0.1],label="Rotation_{}".format(len(fig.axes))) #defines the location of the rotation legend
    rotAx.set_title(rotlabel,fontsize=legendSize1)

    rotAx.set_xticklabels([])
    rotAx.set_xticks([])
    rotAx.set_yticklabels([])
    rotAx.set_yticks([])
    rotAx.set_xlim([-5,5])
    rotAx.set_ylim([-5,5])

    rotaxleg = np.array([0,0.5,1])*(rotRange[1]-rotRange[0]) + rotRange[0]

    # populate the rotation legend
    for i in range(len(rotaxleg)):
        theta = np.linspace(0,2*np.pi,1000)
        radius = 50   # seems arbitrary?
        verts = np.zeros((len(theta),2))
        verts = np.column_stack([radius * np.cos(theta), radius * np.sin(theta)])

        phase = rotaxleg[i]
        shift = np.int(phase*(len(theta)-1))
        verts[shift,:] = [radius*notchScale* np.cos(2*np.pi*phase),radius*notchScale*np.sin(2*np.pi*phase)]
        scaleA = np.maximum(np.cos(2*np.pi*phase)**2,np.sin(2*np.pi*phase)**2)
        rotAx.scatter(0,0,c=[colout(0.5)],
                         edgecolors="black",s=areaRange[0]*scaleA,alpha=0.6,marker=verts,linewidths=0.5)

    rotAx.text(2.75*rottextscale,-0.75*rottextscale,"0",fontsize=legendSize1)
    rotAx.text(-5*rottextscale,-0.75*rottextscale,"100",fontsize=legendSize1)
    rotAx.text(-1*rottextscale,3*rottextscale,"50",fontsize=legendSize1)

    rotAx.set_frame_on(False)

    if xlims is not None:
        axis.set_xlim(xlims)

    if ylims is not None:
        axis.set_ylim(ylims)

    #fig.set_tight_layout({'rect':[0,0,1,1]})

    if file is not None:
        fig.savefig(file)

    return fig,axis,colAx, rotAx





