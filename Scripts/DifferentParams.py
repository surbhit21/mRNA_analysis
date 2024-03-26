from fractions import Fraction
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import numpy as np
import os
from SNSPlottingWidget import SNSPlottingWidget
import Total_AMPAR
colors = {"surf":'#005f73',"cyto" : '#9b2226',"spine" : '#CA6702'}
v_pfit = 1.4e-3
d_pfit = 0.22
def getautostring(val):
    print(val,val.as_integer_ratio(),Fraction(val).limit_denominator(100000).as_integer_ratio())
    num,den = Fraction(val).limit_denominator(100000).as_integer_ratio()
    if num > den:
        return "x {}".format(num)
    elif num == den:
        return ""
    else:
        return "/ {}".format(den)

def_colo = "k"
colors_arr = ["#f72585","#b5179e","#7209b7","#560bad","#480ca8","#3a0ca3","#3f37c9","#4361ee","#4895ef","#4cc9f0"]
def ExploreParameter(V_p_list = [],
                     D_p_list = [],
                     t_half_list = [],
                     V_pinit=v_pfit,
                     D_pinit=d_pfit,
                     t_half_init = 3.12,
                     Norm_by_own_origin=True):
    op_folder = "./Figures/Protein/"
    plt_widget = SNSPlottingWidget()
    plt_widget.CreateFolderRecursive(op_folder)
    fsize = 20
    x_init, yi_init = Total_AMPAR.RunSSProtein(D_pinit, V_pinit)
    ymin = 0
    ymax = -10
    xmin,xmax = 0, 500
    suff = "Own"
    if not V_p_list == []:
        fig,ax = plt.subplots(figsize = (8,6),nrows = 1,ncols = 1)
        ax.plot(x_init,yi_init/yi_init[0],label=r"$V_p$",color=def_colo)

        for vdx,vp in enumerate(V_p_list):
            x, yi = Total_AMPAR.RunSSProtein(D_pinit, V_pinit*vp)
            if Norm_by_own_origin:
                ax.plot(x, yi / yi[0], label=r"$V_p {}$".format(getautostring(vp)),color=colors_arr[vdx])
            else:
                ax.plot(x, yi/yi_init[0] , label=r"$V_p {}$".format(getautostring(vp)),color=colors_arr[vdx])
                suff = "Original"
            if max(yi/yi_init) > ymax:
                ymax = max(yi/yi_init)
            # getautostring(vp)
        plt.legend()

        plt.xlim([xmin,xmax])
        # y_ticks = np.arange(ymin, 1.2*ymax, 0.5)
        # ax.set_yticks(y_ticks)
        # plt.ylim([ymin,1.2 ])
        # plt.yscale("log")
        ax.set_ylabel("Normalized Protein concentration", fontsize=fsize)
        ax.set_xlabel(r"Dendritic distance ($\mu m$)", fontsize=fsize)
        plt_widget.SaveFigures(os.path.join(op_folder,"Effect_VP_Norm_by_{}".format(suff)))
        plt.show()
    if not D_p_list == []:
        fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)

        ax.plot(x_init, yi_init/yi_init[0], label=r"$D_p$",color=def_colo)
        ymax = -10
        for vdx,dp in enumerate(D_p_list):
            x, yi = Total_AMPAR.RunSSProtein(D_pinit*dp, V_pinit )
            if max(yi/yi_init) > ymax:
                ymax = max(yi/yi_init)
            if Norm_by_own_origin:
                ax.plot(x, yi / yi[0], label=r"$D_p {}$".format(getautostring(dp)),color=colors_arr[vdx])
            else:
                ax.plot(x, yi/yi_init[0], label=r"$D_p {}$".format(getautostring(dp)),color=colors_arr[vdx])
                suff = "Original"
        plt.legend()
        plt.xlim([xmin, xmax])
        # y_ticks = np.arange(ymin, 1.2 * ymax, 0.5)
        # ax.set_yticks(y_ticks)
        # plt.ylim([ymin, 1.2 ])
        # plt.yscale("log")
        ax.set_ylabel("Normalized Protein concentration", fontsize=fsize)
        ax.set_xlabel(r"Dendritic distance ($\mu m$)", fontsize=fsize)
        plt_widget.SaveFigures(os.path.join(op_folder, "Effect_DP_Norm_by_{}".format(suff)))
        plt.show()
    if not t_half_list == []:
        fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=1)

        ax.plot(x_init, yi_init/yi_init[0] , label=r"${}$".format(t_half_init),color=def_colo)
        ymax = -10
        for vdx,t_half in enumerate(t_half_list):
            x, yi = Total_AMPAR.RunSSProtein(D_pinit , V_pinit,t_half=t_half)
            if max(yi/yi_init) > ymax:
                ymax = max(yi/yi_init)
            if Norm_by_own_origin:
                ax.plot(x, yi/yi[0] , label=r"${:.2f}$".format(t_half),color=colors_arr[vdx])
            else:
                ax.plot(x, yi / yi_init[0], label=r"${:.2f}$".format( t_half),color=colors_arr[vdx])
                suff = "Original"
        leg = t_half_list.insert(0,t_half_init)
        plt.legend(handles=leg, title="$T_{1/2}(Days)$",
                   fontsize='small', fancybox=True)
        # plt.yscale("log")
        plt.xlim([xmin, xmax])
        # y_ticks = np.arange(ymin, 1.2 * ymax, 0.5)
        # ax.set_yticks(y_ticks)
        # plt.ylim([ymin, 1.2 ])
        ax.set_ylabel("Normalized Protein concentration", fontsize=fsize)
        ax.set_xlabel(r"Dendritic distance ($\mu m$)", fontsize=fsize)
        plt_widget.SaveFigures(os.path.join(op_folder, "Effect_t_half_Norm_by_{}".format(suff)))
        plt.close()
        # plt.show()


own_norm = False
ExploreParameter(D_p_list = [1/10,10,100],V_pinit=0,D_pinit=d_pfit,Norm_by_own_origin=own_norm)
ExploreParameter(V_p_list = [1/2,2,3],D_pinit=d_pfit,Norm_by_own_origin=own_norm)
ExploreParameter(t_half_list = [1/2,2,5],D_pinit=d_pfit,V_pinit=v_pfit,Norm_by_own_origin=own_norm)
own_norm = True
ExploreParameter(D_p_list = [1/10,10,100],V_pinit=0,D_pinit=d_pfit,Norm_by_own_origin=own_norm)
ExploreParameter(V_p_list = [1/2,2,3],D_pinit=d_pfit,Norm_by_own_origin=own_norm)
ExploreParameter(t_half_list = [0.5,1.95,4.35],D_pinit=d_pfit,V_pinit=v_pfit ,Norm_by_own_origin=own_norm)


