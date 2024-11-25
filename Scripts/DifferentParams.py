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
d_pdef = 0.1
def getautostring(val):
    print(val,val.as_integer_ratio(),Fraction(val).limit_denominator(100000).as_integer_ratio())
    num,den = Fraction(val).limit_denominator(100000).as_integer_ratio()
    if num > den:
        if den == 1:
            return "x {}".format(num)
        else:
            return r"x $\frac{{{}}}{{{}}}$".format(num,den)
    elif num == den:
        return r"x $\frac{}{}$".format(num,den)
    else:
        return "/ {}".format(den)

def_colo = "#2ca02c"
def_lw = 6
ot_lw = 3
colors_arr = ["#ADB5BD","#6C757D","#212529","#212529","#4895ef","#4cc9f0",def_colo]
def getPerChange(y_dist):
    return 100 * (y_dist[-1]/y_dist[0] -1 )
def ExploreParameter(V_p_list = [],
                     D_p_list = [],
                     t_half_list = [],
                     V_pinit=v_pfit,
                     D_pinit=d_pfit,
                     t_half_init = 3.12,
                     Norm_by_own_origin=True,
                     display_range=[0,100]):
    op_folder = "./Figures/Protein/"
    plt_widget = SNSPlottingWidget()
    plt_widget.CreateFolderRecursive(op_folder)
    fsize = plt_widget.fsize
    x_init, yi_init = Total_AMPAR.RunSSProtein(D_pinit, V_pinit)
    # breakpoint()
    ymin = 0
    ymax = -10
    xmin,xmax = 0, 500
    xlab,ylab = "Normalized GluA2 density",r"Dendritic distance ($\mu m$)"
    suff = "Own"
    left, bottom, width, height = [0.4, 0.45, 0.5, 0.5]
    fig,ax = plt.subplots(figsize = (8,6),nrows = 1,ncols = 1)
    ax.spines[['right', 'top']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(2)
    ax1 = ax.inset_axes(bounds=[left, bottom, width, height], zorder=4)
    # breakpoint()

    if not V_p_list == []:
        tics = np.arange(0, 1.2, 0.5)
        ax.plot(x_init[display_range],yi_init[display_range]/yi_init[0],label=r"$V_p$",color=colors_arr[-1],linewidth=def_lw)
        ax1.plot(x_init, yi_init / yi_init[0], color=def_colo)
        drop = getPerChange(yi_init[display_range])
        print("drop = ", drop)
        ax.text(x=display_range[-1], y=(yi_init/ yi_init[0])[display_range[-1]],
                s=r"{}%".format(int(drop)), color=def_colo,fontsize=fsize)
        print("change at ~ 300 mu m is: ", (1 - yi_init[300] / yi_init[1]))
        for vdx,vp in enumerate(V_p_list):
            x, yi = Total_AMPAR.RunSSProtein(D_pinit, V_pinit*vp)
            print("change at ~ 300 mu m is: ", (1 - yi[300] / yi[1]))
            if Norm_by_own_origin:
                yi /= yi[0]
            else:
                yi /= yi_init[0]
                suff = "Original"
            ax.plot(x[display_range], yi[display_range] , label=r"$V_p$ {}".format(getautostring(vp)),color=colors_arr[vdx],linewidth=ot_lw)
            ax1.plot(x, yi, color=colors_arr[vdx])
            drop = getPerChange(yi[display_range])
            print("drop = vp for {}".format(vp), drop)
            ax.text(x=display_range[-1], y=yi[display_range[-1]],
                    s=r"{}%".format(int(drop)), color=colors_arr[vdx],fontsize=fsize)
        if max(yi/yi_init) > ymax:
            ymax = max(yi/yi_init)
        # getautostring(vp)
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [1,2,0,3]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
                   fontsize=fsize,labelcolor='linecolor',markerscale=0)
        # plt.ylim([ymin, 1.2])
        # ax.set_xlim(display_range)
        plt.gca().set_ylim(bottom=0)
        ax.set_ylabel(xlab, fontsize=fsize)
        ax.set_xlabel(ylab, fontsize=fsize)
        plt_widget.SaveFigures(os.path.join(op_folder,"Effect_VP_Norm_by_{}".format(suff)))
        plt.show()
    if not D_p_list == []:
        tics = np.arange(0, 1.2, 0.5)
        ax.plot(x_init[display_range], yi_init[display_range]/yi_init[0], label=r"$D_p$",color=colors_arr[-1],linewidth=def_lw)
        ax1.plot(x_init, yi_init / yi_init[0], color=def_colo)
        drop = getPerChange(yi_init[display_range])
        print("drop = ", drop)
        ax.text(x=display_range[-1], y=(yi_init/ yi_init[0])[display_range[-1]],
                s=r"{}%".format(int(drop)), color=def_colo,fontsize=fsize)
        ymax = -10
        for vdx,dp in enumerate(D_p_list):
            x, yi = Total_AMPAR.RunSSProtein(D_pinit*dp, V_pinit )
            if max(yi/yi_init) > ymax:
                ymax = max(yi/yi_init)
            if Norm_by_own_origin:
                yi /= yi[0]
            else:
                yi /= yi_init[0]
                suff = "Original"
            ax.plot(x[display_range], yi[display_range] , label=r"$D_p {}$".format(getautostring(dp)), color=colors_arr[vdx],linewidth=ot_lw)
            ax1.plot(x, yi, color=colors_arr[vdx])
            drop = getPerChange(yi[display_range])
            print("drop = for dp {}".format(vdx), drop)
            ax.text(x=display_range[-1], y=yi[display_range[-1]],
                    s=r"{}%".format(int(drop)),color=colors_arr[vdx],fontsize=fsize)
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [1, 0, 2, 3]
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],labelcolor='linecolor',markerscale=0,fontsize=fsize)
        plt.gca().set_ylim(bottom=0)
        ax.set_ylabel(xlab, fontsize=fsize)
        ax.set_xlabel(ylab, fontsize=fsize)
        plt_widget.SaveFigures(os.path.join(op_folder, "Effect_DP_Norm_by_{}".format(suff)))
        plt.show()
    if not t_half_list == []:

        ax.plot(x_init[display_range], yi_init[display_range]/yi_init[0] , label=r"${}$".format(t_half_init),color=colors_arr[-1],linewidth=def_lw)
        ax1.plot(x_init, yi_init / yi_init[0], color=def_colo)
        drop = getPerChange(yi_init[display_range])
        print("drop = ", drop)
        ax.text(x=display_range[-1], y=(yi_init/ yi_init[0])[display_range[-1]],
                s=r"{}%".format(int(drop)), color=def_colo)
        ymax = -10
        print(t_half_list)
        for vdx,t_half in enumerate(t_half_list):
            x, yi = Total_AMPAR.RunSSProtein(D_pinit , V_pinit,t_half=t_half)
            if max(yi/yi_init) > ymax:
                ymax = max(yi/yi_init)
            if Norm_by_own_origin:
                yi /= yi[0]
            else:
                yi /= yi_init[0]
                suff = "Original"
            ax.plot(x[display_range], yi[display_range] , label=r"${:.2f}$".format(t_half), color=colors_arr[vdx],linewidth=ot_lw)
            ax1.plot(x, yi, color=colors_arr[vdx])
            drop = getPerChange(yi[display_range])
            print("drop = for t_half {}".format(t_half), drop)
            ax.text(x=display_range[-1], y=yi[display_range[-1]],
                    s=r"{}%".format(int(drop)), color=colors_arr[vdx],fontsize=fsize)

        # leg = t_half_list.insert(0,t_half_init)
        handles, labels = plt.gca().get_legend_handles_labels()
        # leg = ax.get_legend()
        # handles[0].set_color(colors_arr[-1])
        # for i in range(1,len(t_half_list)+1):
        #     handles[i].set_color(colors_arr[i])
        order = [1, 2,0, 3]
        plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order],title="$T_{1/2}(Days)$",
                   fontsize=fsize, fancybox=True,labelcolor='linecolor',markerscale=0)
        # plt.legend(handles=leg, title="$T_{1/2}(Days)$",
        #            fontsize=fsize, fancybox=True)
        plt.gca().set_ylim(bottom=0)
        ax.set_ylabel(xlab, fontsize=fsize)
        ax.set_xlabel(ylab, fontsize=fsize)
        plt_widget.SaveFigures(os.path.join(op_folder, "Effect_t_half_Norm_by_{}".format(suff)))
        # plt.close()
        plt.show()


own_norm = False
# do not change the order here or the label order would change in legend
dpl = [1/10,2,10]
vpl = [1/10,1/2,3/2]
thalfl = [0.5,1.95,4.35]
ExploreParameter(D_p_list = dpl,V_pinit=0,D_pinit=d_pdef,Norm_by_own_origin=own_norm)
ExploreParameter(t_half_list =thalfl,D_pinit=d_pdef,V_pinit=0,Norm_by_own_origin=own_norm)
ExploreParameter(V_p_list = vpl,D_pinit=d_pfit,Norm_by_own_origin=own_norm)
# do not change the order here or the label order would change in legend
thalfl = [0.5,1.95,4.35]
own_norm = True
ExploreParameter(D_p_list = dpl,V_pinit=0,D_pinit=d_pdef,Norm_by_own_origin=own_norm)
ExploreParameter(t_half_list =thalfl,D_pinit=d_pdef,V_pinit=0,Norm_by_own_origin=own_norm)
ExploreParameter(V_p_list = vpl,D_pinit=d_pfit,Norm_by_own_origin=own_norm)


