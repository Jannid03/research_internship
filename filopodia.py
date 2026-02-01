#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 10:57:50 2022

@author: eric
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Ellipse, Circle
import matplotlib.colors
from scipy import interpolate
from scipy.signal import fftconvolve
from scipy.spatial import distance as dist
from scipy.integrate import solve_ivp
import time
# import numba
import numexpr as ne
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from joblib import Parallel, delayed
import shapely.affinity
from shapely.geometry import Point
import descartes

#42
seed = 42
np.random.seed(seed=seed)
interest = [20]   #50
verfolgx = []
verfolgy = []
gil = True
update = (gil and False)
start = time.time()
collect = []


a = 0.7  	#0.7
b = 0.7	#0.8
tau = 100 	#12.5
I = 0.0		#0.5
tau_l = 12.5
def func(t,y):
	return [(y[0]-1/20*(y[0]**3)-y[1])/tau_l+I, (y[0]+a-b*y[1])/tau]

def lfunc(t):
    return 7*t*np.exp(-0.5*t)

def sample_gil(FRONT, type, angs = [], lens = [], tims = [], inhib = [], first = False):
    filo= angs.copy()
    # if(not first): print("Start: ", filo[interest])
    lengths = lens.copy()
    times = tims.copy()
    inhibs = inhib.copy()
    factor = 1
    tau = dt * 60
    sigma = np.sqrt(0.102*tau)
    c_l = 0.156
    a = 4.309

    xinds = []
    yinds = []

    for f in range(len(FRONT)):
        wink = np.round(360*startangs[type[f]-1]/(2*np.pi))
        wink = wink + 360 * (wink < 0)

        upper =  wink + 90
        lower = wink - 90

        # print("Upper: ", upper, "    Lower: ", lower)

        ### Initalizing filopodia
        if(first):
            loc_filo = []
            loc_len = []
            loc_tim = []
            loc_inh = []
            number = np.random.poisson(4.842)
            # if(f==interest): print("Num; ", num)
            for i in range(number):
                num = np.random.randint(lower,upper)
                if (num < 0): num += 360
                if (num > 360): num -= 360
                loc_filo.append(num)

                # loc_len.append(np.abs(np.random.normal(0,sigma)))
                loc_len.append(0.1) #0.1 for fitz
                # loc_len.append(4)
                loc_tim.append(0)
                loc_inh.append(-20)
                # if loc_len[-1] > 8: loc_len[-1] = 8
            # if(f in interest): print("First: ", loc_filo)
            filo.append(loc_filo)
            lengths.append(loc_len)
            times.append(loc_tim)
            inhibs.append(loc_inh)
            new = False
        else:
            # if(f == interest): print("Angs before: ", filo[interest])
            # if(f == interest): print("Lengths before: ", lengths[interest])
            # print("Inh: ", inhibs)
            #reaction parameters
            B50 = 0.00272
            c = 950.535
            F = len(filo[f])
            # if(f==interest): print("f: ", f)
            k = [B50/(F+B50)*c,0.124*F]

            prop = np.zeros((360,1))
            # print(prop)
            # print(prop[0])
            # print(filo)
            for i in range(360):
                if ((i > upper) and (i < lower)) or (((i-360) < lower) and (i > upper)) or (((i+360) > upper ) and (i < lower)):
                    prop[i] = 0
                    continue
                elif (len(filo[f]) == 0):
                    prop[i] = i/180
                else:
                    # print("List: ", [np.abs(i-akt) for akt in filo[f]])
                    prop[i] = np.min([np.abs(i-akt) for akt in filo[f]])/180

            prop = prop.flatten()
            # if(f == interest): print("K: ", k)
            # k[0] = B50/(F+B50)*c
            new = False
            if(np.random.uniform(0,1) < k[0]/np.sum(k)):
                # if(f == interest): print("Birth")
                s = np.sum(prop)
                prop = [p/s for p in prop]

                # print("Prop: ", prop)

                draw = np.random.choice(range(360),1,p=prop)
                # if(f==interest): print("Draw: ", draw)
                filo[f].append(draw[0])
                # lengths[f].append(np.abs(np.random.normal(0,sigma)))
                lengths[f].append(0.1)
                times[f].append(0)
                inhibs[f].append(-20)

                # if lengths[f][-1] > 8: lengths[f][-1] = 8

                new = True
                # new = False
                # if(f == interest): print("Draw len: ", lengths[f][-1])
            else:
                
                ind = np.random.randint(0,F)
                # if(f == interest): print("Death of ", ind)
                filo[f].pop(ind)
                lengths[f].pop(ind)
                collect.append(times[f][ind])
                times[f].pop(ind)
                inhibs[f].pop(ind)

        pops=[]
        # if(f==interest): print("Bevor: ", lengths[f])
        for i in range(len(filo[f])-new):
            ### Standard
            # lengths[f][i] += (lengths[f][i]-a) * c_l * tau + np.random.normal(0,sigma)

            ##Time wise
            times[f][i] += dt*60
            # lengths[f][i] = lfunc(times[f][i])

            ##Fitz
            sol = solve_ivp(func, [0,dt*60], [lengths[f][i],inhibs[f][i]])
            lengths[f][i] = sol.y[0][-1]
            inhibs[f][i] = sol.y[1][-1]

            # if(f==interest): print("New len: ", lengths[f][i])

            if (lengths[f][i] <= 0):
                if(first):
                    lengths[f][i] = 0
                else:
                    pops.append(i)
                    # lengths[f][i] = 0
            # if lengths[f][i] > 8:
            #     lengths[f][i] = 8
        
        counter = 0
        # if(f==interest): print("Pops: ", pops)
        for p in pops:
            lengths[f].pop(p-counter)
            filo[f].pop(p-counter)
            inhibs[f].pop(p-counter)
            counter += 1
    
        ###Berechnung wo
        tempx = []
        tempy = []

        for i in range(len(lengths[f])):
            tempx.append(FRONT[f][0]+np.cos(np.deg2rad(filo[f][i]))*lengths[f][i]*CONVFAC)
            tempy.append(FRONT[f][1]+np.sin(np.deg2rad(filo[f][i]))*lengths[f][i]*CONVFAC)

        xinds.append(tempx)
        yinds.append(tempy)

        # if(f == interest): print("End Angs: ", filo[interest])
        # if(f == interest): print("End Len: ", lengths[interest])


    return xinds, yinds, filo, lengths, inhibs, times

def update_parameters(FRONTLOC, filo, lengths, xfil, yfil):

    for f in range(len(FRONTLOC)):
        for i in range(len(lengths[f])):
            # if(f == interest): print("X: ", xfil[f][i], "     Real: ", FRONTLOC[f][0], "         Y: ", yfil[f][i], "    Real: ",FRONTLOC[f][1])
            # if(f == interest): print("X: ", xfil[f][i]-FRONTLOC[f][0], "         Y: ", yfil[f][i]-FRONTLOC[f][1])
            if((f == interest) and (i == 0)): print("X: ", xfil[f][i], "     Y: ", yfil[f][i])
            if(f == interest): print("X diff: ", (xfil[f][i]-FRONTLOC[f][0])**2, "    Y diff: ", (yfil[f][i]-FRONTLOC[f][1])**2)
            new_len = np.sqrt((xfil[f][i]-FRONTLOC[f][0])**2 + (yfil[f][i]-FRONTLOC[f][1])**2)
            if(f == interest): print("New length (UP): ", new_len)
            if(f == interest): print("Zahlen (UP): ", xfil[f][i], "     ", FRONTLOC[f][0], "       ", xfil[f][i]-FRONTLOC[f][0])
            if(f == interest): print("Bruch (UP): ", (xfil[f][i]-FRONTLOC[f][0])/new_len)
            new_ang = np.arccos((xfil[f][i]-FRONTLOC[f][0])/new_len) * np.sign(yfil[f][i]-FRONTLOC[f][1])

            if(np.arccos((xfil[f][i]-FRONTLOC[f][0])/new_len) < 0):
                new_ang = 360 - new_ang

            if(f == interest): print("New ang (UP): ", new_ang)
            lengths[f][i] = new_len/CONVFAC
            filo[f][i] = np.rad2deg(new_ang)       


def average_vec(angs, lengths):
    
    erg = []
    for i in range(len(angs)):
        if(len(angs[i]) == 0): 
            erg.append(0j+0)
            # print("Hier: ", i)
            continue

        temp = []
        ### imaginary number -> xj + y -> x for x coordinate and y for y coordinate growth
        for j in range(len(angs[i])):
            temp.append(np.exp(1j*np.deg2rad(angs[i][j]))*lengths[i][j]*CONVFAC)
        erg.append(np.mean(temp))
    
    # print("Erg: ", erg[interest])
    return erg

def getAcc(pos):
    """
    calculates the acceleration for the simulation of repulsive particles
    used for the dynamic interaction of the L-cell bundles
    """
    # positions r = [x,y,z] for all particles
    x = np.array([pos[:, 0]])
    y = np.array([pos[:, 1]])

    # matrix that stores all pairwise particle separations: r_j - r_i
    dx = (x.T - x)
    dy = (y.T - y)

    # matrix that stores 1/r^3 for all particle pairwise particle separations
    inv_r3 = (dx**2 + dy**2 + 1e-18)**(-1.5)

    ax = np.sum((dx * inv_r3), axis=1)
    ay = np.sum((dy * inv_r3), axis=1)

    # pack together the acceleration components
    a = 5e3 * np.vstack((ax, ay)).T
    return a


def heatmap(data, row_labels, col_labels, xlabel='', ylabel='', extralines=True, ax1=None,
            cbar_kw={}, cbarnorm=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax1:
        ax1 = plt.gca()

    # Plot the heatmap
    im1 = ax1.imshow(data, norm=cbarnorm, **kwargs)

    # Create colorbar
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = ax1.figure.colorbar(im1, cax=cax, **cbar_kw)
    # cbar.ax.set_yticklabels(['repulsion', 'attraction'])
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom", fontsize=10)

    # We want to show all ticks...
    ax1.set_xticks(np.arange(data.shape[1]))
    ax1.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax1.set_xticklabels(col_labels)
    ax1.set_yticklabels(row_labels)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    # Let the horizontal axes labeling appear on top.
    ax1.tick_params(top=False, bottom=True,
                    labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax1.get_xticklabels(), rotation=0, ha="center", va='center',
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for _, spine in ax1.spines.items():
        spine.set_visible(False)

    ax1.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax1.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax1.grid(which="minor", color="w", linestyle='-', linewidth=2)
    if extralines:
        ax1.axvline(10.5, linestyle='--', color='k')
        ax1.axvline(21.5, linestyle='--', color='k')
        ax1.axvline(32.5, linestyle='--', color='k')
        ax1.axvline(43.5, linestyle='--', color='k')
        ax1.axvline(54.5, linestyle='--', color='k')
        ax1.axvline(65.5, linestyle='--', color='k')

    ax1.tick_params(which="minor", bottom=False, left=False)

    return im1, cbar


def remove_ticks_and_box(ax1):
    """
    helper function for plotting
    """
    ax1.tick_params(axis='x', which='both', bottom=False,
                    top=False, labelbottom=False)
    ax1.tick_params(axis='y', which='both', right=False,
                    left=False, labelleft=False)
    for pos in ['right', 'top', 'bottom', 'left']:
        ax1.spines[pos].set_visible(False)


def potential_gauss(loc, mu_loc, cov, amp):
    """
    bivariate Gaussian to model repulsion from central L3 area
    ouput has shape (M, N) = (number of samples in loc,
                              number of mean positions in mu_loc)
    """
    return amp*np.exp(-0.5*(loc[0, :, None] - mu_loc[0, None, :])**2/cov[0, 0]
                      - 0.5*(loc[1, :, None] - mu_loc[1, None, :])**2/cov[1, 1])


def potential_parabola(loc, mu_loc, width, amp):
    """
    square parabola with support defined by width,
    different extent for x-axis and y-axis
    """
    return amp*np.maximum(-np.sum((loc[:, :, None]-mu_loc[:, None, :])**2
                                  / width[:, None, None]**2, axis=0) + 1, 0)


def rotate_coord(x, y, angle):
    """
    angle in rad
    """
    return x*np.cos(angle) - y*np.sin(angle), x*np.sin(angle) + y*np.cos(angle)


def create_starting_grid_noLcell(center_loc_x, center_loc_y):
    """
    create grid of bundles in the noLcell situation, elliptical bundles
    at 22 hAPF, without equator
    """
    alpha = -np.pi/3
    ell_a = 2.58 * 24 / 2
    ell_b = 2.96 * 24 / 2

    yrec = np.array([ell_a*np.cos(-(alpha + (n - 1)*(np.pi - 2*alpha)/5) + np.pi)
                    for n in np.arange(1, 7)])
    xrec = np.array([ell_b*np.sin(-(alpha + (n - 1)*(np.pi - 2*alpha)/5) + np.pi)
                    for n in np.arange(1, 7)])

    recep_y = yrec[:, None] + center_loc_y[None, :]
    recep_x = xrec[:, None] + center_loc_x[None, :]
    return recep_x, recep_y, center_loc_x, center_loc_y


def create_starting_grid2(center_loc_x, center_loc_y, hour, alpha=-0.23):
    """
    for real L cell grid
    create grid at 22 hAPF, without equator
    """
    if species == 'dros':
        if hour < 35:
            ell_a = 32  # 35#*np.mean(np.diff(Y)[np.diff(Y) > 0])
            ell_b = 55  # 60#*np.mean(np.diff(X)[np.diff(X) > 0])
        else:
            ell_a = 34  # 38#*np.mean(np.diff(Y)[np.diff(Y) > 0])
            ell_b = 60  # 65#*np.mean(np.diff(X)[np.diff(X) > 0])
        yrec = np.array([ell_a*np.cos(-(alpha + (n - 1) *
                        (np.pi - 2*alpha)/5) + np.pi) for n in np.arange(1, 7)])
        xrec = np.array([ell_b*np.sin(-(alpha + (n - 1) *
                        (np.pi - 2*alpha)/5) + np.pi) for n in np.arange(1, 7)])
        xrec_rot, yrec_rot = rotate_coord(
            xrec, yrec, -7*np.pi/180)  # - 20*np.pi/180)

    elif species == 'bib':
        r = 30
        # angles = np.linspace(-np.pi/2, 3/2*np.pi+np.pi/3, 6, endpoint=False)
        angles = np.linspace(0, 2*np.pi, 6, endpoint=False)
        xrec_rot = r*np.sin(angles)
        yrec_rot = r*np.cos(angles)

    recep_y = yrec_rot[:, None] + center_loc_y[None, :]
    recep_x = xrec_rot[:, None] + center_loc_x[None, :]

    return recep_x, recep_y, center_loc_x, center_loc_y


def calc_closest_point_on_ellipse(a, b, point):
    # assume that the center of the ellipse is at (0, 0)
    # a is the axis in x-direction, b in y-direction
    xr = np.sqrt(a**2 * b**2 / (b**2 + a**2 *
                 (point[:, :, 1]/point[:, :, 0])**2))
    yr = point[:, :, 1]/point[:, :, 0] * xr
    return np.sign(point[:, :, 0])*xr, np.sign(point[:, :, 1])*np.abs(yr)


def tanh_cust(x, x_half, sl):
    return 1/2 * (1 + np.tanh(2*sl*(x - x_half)))


def generate_indmat(xind, yind, fr_x, fr_y, r, frontang, REC_X_INTER, REC_Y_INTER, circ_width):
    rec_x, rec_y = np.ravel(REC_X_INTER), np.ravel(REC_Y_INTER)
    rec_x, rec_y = rec_x[(rec_x - fr_x)**2 + (rec_y - fr_y)**2 <= 1.2*r **
                         2], rec_y[(rec_x - fr_x)**2 + (rec_y - fr_y)**2 <= 1.2*r**2]
    ind = np.zeros(np.shape(xind), dtype='bool')

    # ind1 = ((xind - fr_x)**2 + (yind - fr_y)**2 <= r**2)
    ind1 = ne.evaluate('((xind-fr_x)**2 + (yind-fr_y)**2)<=r**2')
    ind2 = (np.cos(np.arctan2(yind[ind1]-fr_y, xind[ind1]-fr_x) - frontang)
            > np.cos(np.pi/180*circ_width))  # 0.643) #0.25882) #0.77)
    ind[yind[ind1][ind2], xind[ind1][ind2]] = True
    """radius = 0.5 * CONVFAC
    ind3 = np.zeros(np.shape(xind), dtype='bool')
    for rx, ry in zip(rec_x, rec_y):
        ind3 += ne.evaluate('((xind-rx)**2 + (yind-ry)**2)<=radius**2')
    return np.logical_and(ind, ~ind3)
    """
    return ind ### all indices that are withing a certain radius and in the right angle direction


# sample_roi(density, FRONTLOC, np.tile(startangs, (len(Xrec[0]), 1)).T + ang_noise, region, xind, yind, POS, Xrec, Yrec, np.tile(
#         radii, (len(Xrec[0]), 1)).T, np.tile(circ_width, (len(Xrec[0]), 1)).T, sampfac = 1)

def sample_roi(dat2_inter, FRONTLOC, frontang, region, xind, yind, POS, REC_X_INTER, REC_Y_INTER, r, c, sampfac):
    # xfil, yfil = np.zeros((n_fil, nr_of_rec)), np.zeros((n_fil, nr_of_rec))
    xfil, yfil = [], []
    allinds = np.zeros(1023*2047, dtype='bool')
    # print("Len of Frontloc: ", len(FRONTLOC))
    for irec in range(len(FRONTLOC)):
        if region == 'circle':
            ind = (
                np.sqrt((xind - FRONTLOC[irec, 0])**2 + (yind - FRONTLOC[irec, 1])**2) <= r)
        elif region == 'circ_segment':
            # frontang = np.arctan2(meanvec_old[irec].imag, meanvec_old[irec].real)
            ind = generate_indmat(xind, yind, FRONTLOC[irec, 0], FRONTLOC[irec, 1], np.ravel(
                r)[irec], np.ravel(frontang)[irec], REC_X_INTER, REC_Y_INTER, np.ravel(c)[irec])
        # if(irec == 42): print("Ind: ", ind)
        # if(irec == 42): print("Inter: ", dat2_inter[ind])
        histog, bins = np.histogram(dat2_inter[ind], bins=10000)
        # if(irec == 42): print("Histog num: ", np.sum(histog > 0))

        if np.sum(histog > 0) > 0:
            histog = np.r_[0, histog, 0]
            bins = np.r_[bins[0] - (bins[1]-bins[0]),
                         bins, bins[-1] + (bins[-1]-bins[-2])]

            rr = np.random.rand(int(n_fil/sampfac))
            if density_type == 'fil_attract':
                cs = (np.cumsum(histog)/np.sum(histog))**4
                ind_res = np.array([np.where(cs > rtemp)[0][0] for rtemp in rr])
                vals = ((bins[1:] + bins[:-1])/2)[ind_res + 1]
                # this +1 makes sure that vals is slightly above 1
                # for a flat density = 1 above the sub-1 no-go zones
                # so the chosen locations will in the flat area, not the edge of
                # a no-go zone at e.g. 0.9999
            else:
                cs = (1 - np.cumsum(histog)/np.sum(histog))**4
                ind_res = np.array([np.where(cs < rtemp)[0][0] for rtemp in rr])
                vals = ((bins[1:] + bins[:-1])/2)[ind_res - 1]

            # print(vals)
            D_vals = np.abs(vals[:, None] - dat2_inter[ind][None, :])
            # all occurences of max value
            a, b = np.where(D_vals == D_vals.min(axis=1)[:, None])
            maxind = np.array([np.random.choice(b[a == i])
                              for i in range(int(n_fil/sampfac))])
            xfil = np.hstack((xfil, POS[0][ind][maxind]))
            yfil = np.hstack((yfil, POS[1][ind][maxind]))
        elif np.sum(histog > 0) == 0:
            # flashlight outside density, happens in ablation condition, just put filopodia at front
            xfil = np.hstack((xfil, np.repeat(FRONTLOC[irec, 0], int(n_fil/sampfac))))
            yfil = np.hstack((yfil, np.repeat(FRONTLOC[irec, 1], int(n_fil/sampfac))))
        else:
            print('I think that should never happen')
            # many values in one bin
            # -> constant density, just pick random positions
            randint = np.random.randint(0, np.sum(ind), n_fil)
            xfil = np.hstack((xfil, POS[0][ind][randint]))
            yfil = np.hstack((yfil, POS[1][ind][randint]))

        if irec == 2*162 + 73:
            allinds = allinds + np.ravel(ind)
    return xfil, yfil, allinds


def distance_to_exp(firstpos, pos_eval, Xc, Yc, v1, v2, ind_abl, goal_loc):
    closest_bundle = np.zeros(np.shape(goal_loc))
    a, b = LCELL_SIZE[0]*CONVFAC/2
    for kk in range(len(goal_loc)):
        # D_bundle = dist.cdist([pos_eval[kk, :]], [[x, y] for x, y in zip(np.ravel(Xc), np.ravel(Yc))], metric='euclid')
        # closest_bundle[kk, 0] = np.ravel(Xc)[D_bundle.argmin(axis=1)[0]]
        # closest_bundle[kk, 1] = np.ravel(Yc)[D_bundle.argmin(axis=1)[0]]

        Xs_ell, Ys_ell = calc_closest_point_on_ellipse(
            a, b, pos_eval[None, kk, :] - np.array([[x, y] for x, y in zip(Xc, Yc)])[:, None, :])
        Xs_circ, Ys_circ = calc_closest_point_on_ellipse(1.2*b, 1.2*b, pos_eval[None, kk, :] - np.array(
            [[x, y] for x, y in zip(Xc - 0.333*v2[0], Yc - 0.333*v2[1])])[:, None, :])
        D_bundles_ell = dist.cdist([pos_eval[kk, :]], [[x, y] for x, y in zip(
            Xc + Xs_ell[:, 0], Yc + Ys_ell[:, 0])], metric='euclid')
        D_bundles_circ = dist.cdist([pos_eval[kk, :]], [[x, y] for x, y in zip(
            Xc - 0.333*v2[0] + Xs_circ[:, 0], Yc - 0.333*v2[1] + Ys_circ[:, 0])], metric='euclid')
        min_ell = D_bundles_ell.argmin(axis=1)
        min_circ = D_bundles_circ.argmin(axis=1)
        if D_bundles_circ.min(axis=1) < D_bundles_ell.min(axis=1):
            closest_bundle[kk, 0] = Xc[min_circ]
            closest_bundle[kk, 1] = Yc[min_circ]
        else:
            closest_bundle[kk, 0] = Xc[min_ell]
            closest_bundle[kk, 1] = Yc[min_ell]

    correct = np.isclose(closest_bundle, goal_loc).all(axis=1)
    # restrict output to the central bundles
    bundles = np.swapaxes(np.stack((np.repeat(Xc, 6), np.repeat(
        Yc, 6)), axis=1).reshape((len(Xc), 6, 2)), 0, 1)
    bundles = bundles.reshape((np.prod(np.shape(bundles)[:2]), 2))

    if noLcell:
        keep = np.where((bundles[:, 0] > 300) * (bundles[:, 0] < 1710)
                        * (bundles[:, 1] > 260) * (bundles[:, 1] < 785))[0]
    else:
        keep = np.where((bundles[:, 0] > 160) * (bundles[:, 0] < 1850)
                        * (bundles[:, 1] > 120) * (bundles[:, 1] < 925))[0]

    ind = np.zeros(len(pos_eval), dtype='bool')
    ind[keep] = True
    correct_filtered = correct[ind * ~ind_abl]
    correct_temp = np.array(correct.copy(), dtype='int')
    correct_temp[ind_abl] = -1
    correct_full = np.reshape(correct_temp[ind], (6, -1))
    return correct_filtered, correct_full


def run_gc(c_st, sl_st, A_magnet, perc_ablation, c_magnet=35, sl_magnet=0.4):
    # create starting grid
    v1, v2 = vec_opt[0, :2], vec_opt[0, 2:]
    # print("V1: ", v1)
    # print("V2: ", v2)
    mini, maxi = -20, 20
    n1, n2 = np.meshgrid(np.arange(mini, maxi), np.arange(mini, maxi))

    center_loc_x = (n1*v1[0] + n2*v2[0]).flatten()
    center_loc_y = (n1*v1[1] + n2*v2[1]).flatten()

    # print("X: ", center_loc_x)
    # print("Y: ", center_loc_y)
    # must reduce the area, kernel can't handle the whole thing
    if noLcell:
        ind = (center_loc_x > 0)*(center_loc_x < 2100) * \
            (center_loc_y > 0)*(center_loc_y < 1100)
    else:
        # ind = (center_loc_x > 0)*(center_loc_x < 2100) * \
        #     (center_loc_y > 0)*(center_loc_y < 1100)
        ind = (center_loc_x > 500)*(center_loc_x < 1000) * \
            (center_loc_y > 500)*(center_loc_y < 750)
    Xint, Yint = center_loc_x[ind], center_loc_y[ind]
    # bundle_pos = np.c_[Xint.copy(), Yint.copy()]
    # vel = np.zeros((len(Xint), 2))
    # acc = getAcc(bundle_pos)

    if species == 'bib':
        Xrec, Yrec, Xc, Yc = create_starting_grid2(
            Xint, Yint, 20, -45 * np.pi/180)
    elif species == 'dros':
        if noLcell:
            Xrec, Yrec, Xc, Yc = create_starting_grid_noLcell(Xint, Yint)
        else:
            Xrec, Yrec, Xc, Yc = create_starting_grid2(Xint, Yint, 20)


    if jitter_heels > 0:
        Xrec += jitter_heels * CONVFAC * np.random.randn(*np.shape(Xrec))
        Yrec += jitter_heels * CONVFAC * np.random.randn(*np.shape(Yrec))
    ang_noise = std_ang*np.random.randn(len(Xrec[0]))
    meanvec_old = np.ravel(np.tile(speeds, (len(Xrec[0]), 1)).T * 2/3*np.tile(radii, (len(Xrec[0]), 1)).T *
                           np.sin(np.pi/180 * np.tile(circ_width, (len(Xrec[0]), 1)).T) /
                           (np.pi/180 * np.tile(circ_width, (len(Xrec[0]), 1)).T) * np.exp(1j*(np.tile(startangs, (len(Xrec[0]), 1)).T + ang_noise)))
    FRONTLOC = np.c_[np.ravel(Xrec), np.ravel(Yrec)]
    types = np.repeat(np.array([1,2,3,4,5,6]), 9)
    remove = np.zeros(len(FRONTLOC), dtype='bool')
    if single_ablation:
        remove[15] = True
        FRONTLOC[remove, :] = 1e4, 1e4
    elif perc_ablation > 0:
        rem = np.random.choice(len(FRONTLOC), int(
            len(FRONTLOC) * perc_ablation/100), replace=False)
        remove[rem] = True
        FRONTLOC[remove, :] = 1e4, 1e4
    FRONTLOCS_TIME = FRONTLOC.copy()
    xind, yind = np.meshgrid(np.arange(2*Nx-1), np.arange(2*Ny-1))
    xmin2, xmax2, ymin2, ymax2 = 0, 2*Nx, 0, 2*Ny
    POS = np.meshgrid(np.linspace(xmin2, xmax2, 2*Nx-1),
                      np.linspace(ymin2, ymax2, 2*Ny-1))
    POS = np.array(POS)

    # calculate goal locations
    goal_loc = np.swapaxes(np.stack(
        (np.repeat(Xc, 6), np.repeat(Yc, 6)), axis=1).reshape((len(Xc), 6, 2)), 0, 1)
    goal_loc[0, :, :] += -v1[None, :]
    goal_loc[1, :, :] += -v1[None, :] + v2[None, :]
    goal_loc[2, :, :] += 2*v2[None, :] - v1[None, :]
    goal_loc[3, :, :] += v2[None, :]
    goal_loc[4, :, :] += v1[None, :]
    goal_loc[5, :, :] += -v2[None, :] + v1[None, :]
    goal_loc = goal_loc.reshape((np.prod(np.shape(goal_loc)[:2]), 2))

    ### Correcting positions of goal
    # for i in range(len(goal_loc)):
    #     if(goal_loc[i][0] < 500):
    #         goal_loc[i][0] = 1000 + goal_loc[i][0] - 500
    #     elif(goal_loc[i][0] > 1000):
    #         goal_loc[i][0] = 500 + goal_loc[i][0] - 1000
            
    #     if(goal_loc[i][1] < 500):
    #         goal_loc[i][1] = 750 + goal_loc[i][1] - 500
    #     elif(goal_loc[i][1] > 750):
    #         goal_loc[i][1] = 500 + goal_loc[i][1] - 750

    # print("Goal for interest: ", goal_loc[interest])

    # create initial filopodial distribution
    density_fil = np.ones((2*Ny-1, 2*Nx-1))
    # density_fil = np.random.rand(*(2*Ny-1, 2*Nx-1))
    POS2 = np.reshape(np.meshgrid(
        np.arange(2*Nx-1), np.arange(2*Ny-1)), (2, -1))
    if density_type == 'full':
        kernel_heel = potential_parabola(POS2, np.array([[Nx], [Ny]]),
                                         np.mean(rm_heels, axis=0)[0]*np.ones(2)*CONVFAC*2, 1).reshape((2*Ny-1, 2*Nx-1))
        hist_rec = np.histogram2d(np.ravel(Xrec), np.ravel(Yrec),
                                  bins=[np.linspace(xmin2, xmax2, 2*Nx),
                                        np.linspace(ymin2, ymax2, 2*Ny)])[0].T
        density_heels = fftconvolve(hist_rec, kernel_heel, mode='same')
        density_heels /= np.max(density_heels)

    match species:
        case 'dros':
            if noLcell:
                lcell_ell = potential_parabola(
                    POS2, np.c_[Xint, Yint].T, LCELL_SIZE[0, :]/2*CONVFAC, 0)
                density_lcell = np.sum(lcell_ell, axis=1).reshape(
                    (2*Ny-1, 2*Nx-1))
            else:
                lcell_ell = potential_parabola(
                    POS2, np.c_[np.append(Xint, Xint[0]-500+1000), np.append(Yint, Yint[0])].T, LCELL_SIZE[0, :]/2*CONVFAC, 5*np.max(density_fil))
                lcell_circ = potential_parabola(
                    POS2, np.c_[np.append(Xint, Xint[0]-500+1000) - 0.333*v2[0], np.append(Yint, Yint[0]) - 0.333*v2[1]].T, 1.2*np.ones(2)*LCELL_SIZE[0, 1]/2*CONVFAC, 5*np.max(density_fil))
                density_lcell = np.sum(lcell_ell, axis=1).reshape(
                    (2*Ny-1, 2*Nx-1)) + np.sum(lcell_circ, axis=1).reshape((2*Ny-1, 2*Nx-1))
        case 'bib':
            lcell_ell = potential_parabola(
                POS2, np.c_[Xint, Yint].T, np.array([30, 30]), 10)
            density_lcell = np.sum(lcell_ell, axis=1).reshape(
                (2*Ny-1, 2*Nx-1))

    if include_nogozones or density_type == 'full':
        density_nogo_heels = np.zeros(np.shape(xind))
        radius = 0.5 * CONVFAC
        for rx, ry in zip(np.ravel(Xrec), np.ravel(Yrec)):
            density_nogo_heels[((xind-rx)**2 + (yind-ry)**2) <= radius**2] = 100
    
    match density_type:
        case 'full':
            density = density_heels + density_lcell + density_nogo_heels + density_fil
        case 'fil_attract':
            if include_nogozones:
                density = density_fil - density_lcell - density_nogo_heels
            else:
                density = density_fil
        case 'fil_repuls':
            if include_nogozones:
                density = density_fil + density_lcell + density_nogo_heels
            else:
                density = density_fil
        case 'constant':
            if include_nogozones:
                density = density_fil + density_nogo_heels + density_lcell
            else:
                density = density_fil
        case _:
            raise ValueError('Density type not implemented')
    
    # print("Dens: ", density)
    # draw filopodia

    if(gil):
        xfil, yfil, angs, lengths, inhib, times = sample_gil(FRONTLOC, types, [], [], [], [], True)
    else:
        xfil, yfil, allinds = sample_roi(density, FRONTLOC, np.tile(startangs, (len(Xrec[0]), 1)).T + ang_noise, region, xind, yind, POS, Xrec, Yrec, np.tile(
            radii, (len(Xrec[0]), 1)).T, np.tile(circ_width, (len(Xrec[0]), 1)).T, sampfac = 1)
    
        xfil = np.reshape(xfil, (-1, n_fil))
        yfil = np.reshape(yfil, (-1, n_fil))
   

        angl = np.arctan2(yfil-FRONTLOC[:, 1][:, None],
                        xfil-FRONTLOC[:, 0][:, None])
        radi = np.sqrt((yfil-FRONTLOC[:, 1][:, None])
                    ** 2 + (xfil-FRONTLOC[:, 0][:, None])**2)

    # print("Xfil gil: ", xfil)
    # print("Yfil gil: ", yfil)
    # print("angs gil: ", angs)
    # print("lengths gil: ", lengths)
    # print("Xfil: ", xfil)
    # print("Frontloc: ", FRONTLOC)
    # print("Xfil length: ", len(xfil))
    # print("FRONTLOCS: ", len(FRONTLOC))               
    # print("radi: ", radi)
    # print("angs: ", angs)
    # print("angl: ", angl)
    
    # cubic splines for the size of the L-cells
    # ell_large = LCELL_SIZE[:, 0]
    # ell_small = LCELL_SIZE[:, 1]
    # ell_large_interp = interpolate.splrep(hours, ell_large)
    # ell_small_interp = interpolate.splrep(hours, ell_small)
    
    outside_fils = []
    # c_st = 22 + 8*np.random.rand(len(FRONTLOC))
    # sl_st = 0.45
    if species == 'dros':
        axis_lim = [800, 1300, 300, 650]
        # axis_lim = [600, 1500, 200, 750]
    elif species == 'bib':
        axis_lim = [700, 1100, 350, 600]

    for ii, tt in enumerate(np.arange(20, 40+dt, dt)):
        print('Current time ' + str(np.round(tt, 2)) + ' hAPF')
        print("Position: ", FRONTLOC[interest])
        # lcell_int_large, lcell_int_small = LCELL_SIZE[:, 0], LCELL_SIZE[:, 1]
        lcell_int_large, lcell_int_small = LCELL_SIZE[0, 0], LCELL_SIZE[0, 1]
        #interpolate.splev(tt, ell_large_interp), interpolate.splev(tt, ell_small_interp)
        ###############################
        # plot stuff and create movie
        if create_movie:
            fig, ax = plt.subplots(2, 1, figsize=(7, 6), gridspec_kw={
                                   'height_ratios': [1, 20]})
            ax[0].barh([1], [hours[-1]], color='w', edgecolor='k')
            ax[0].barh([1], [tt], color='k')
            ax[0].set_title('Developmental time (hAPF)')
            ax[0].tick_params(axis='y', which='both', right=False,
                              left=False, labelleft=False)
            ax[0].set_xticks(np.arange(20, hours[-1]+1, 5))
            for pos in ['right', 'top', 'bottom', 'left']:
                ax[0].spines[pos].set_visible(False)

            ax[0].set_xlim([19.9, hours[-1]+0.1])
            
            if include_nogozones or density_type == 'full':
                if density_type == 'fil_attract':
                    density_plot = density + density_nogo_heels + density_lcell
                else:
                    density_plot = density - density_nogo_heels #- density_lcell
            else:
                density_plot = density

            ax[1].imshow(density_plot, origin='lower',
                         vmin=np.min(
                             (density_plot)[axis_lim[2]:axis_lim[3], axis_lim[0]:axis_lim[1]]),
                         vmax=np.max((density_plot)[axis_lim[2]:axis_lim[3], axis_lim[0]:axis_lim[1]]))

            if(not gil):
                ax[1].imshow(np.reshape(allinds, (1023, 2047)),
                         origin='lower', alpha=0.2)

            
                color = np.repeat('r', 54)
                color[interest] = 'y'
                ax[1].scatter(*FRONTLOC.T, c=color)

                # ax[1].scatter(Xc-25, Yc, s=50, c='k')
                # print("Xfil: ", xfil)
                col = np.repeat('w', 54)
                col[interest] = 'y'

                col = np.repeat(col, 10)
                x_draw = []
                y_draw = []
                for v in xfil:
                    for i in v:
                        c = i.copy()
                        if(c < 500):
                            c = 1000 + c - 500
                        elif(c > 1000): 
                            c = 500 + c - 1000
                        x_draw.append(c)

                for v in yfil:
                    for i in v:
                        c = i.copy()
                        if(c < 500):
                            c = 750 + c - 500
                        elif(c > 750): 
                            c = 500 + c - 750
                        
                        y_draw.append(c)
                ax[1].scatter(x_draw, y_draw, s=5, c=col)
                # ax[1].scatter(xfil[interest], yfil[interest], s=5, c = 'y')
                # ax[1].scatter(xfil, yfil, s=5, c='w')
            else:

                ###Draw fronts
                x_drawf = []
                y_drawf = []
                for i in range(len(FRONTLOC)):
                # for i in [interest]:
                #     #####Hier
                    # ax[1].plot(FRONTLOCS_TIME[:, i, 0],
                    #            FRONTLOCS_TIME[:, i, 1], 'r-', lw=2)

                    locx = FRONTLOC[i][0].copy()
                    locy = FRONTLOC[i][1].copy()
                    if(locx < 500):
                        locx = 1000 + FRONTLOC[i][0] - 500
                    elif(locx > 1000):
                        locx = 500 + FRONTLOC[i][0] - 1000
                
                    if(locy < 500):
                        locy = 750 + FRONTLOC[i][1] - 500
                    elif(locy > 750):
                        locy = 500 + FRONTLOC[i][1] - 750
                    
                    x_drawf.append(locx)
                    y_drawf.append(locy)

                color = np.repeat('r', 54)
                color[interest] = 'y'
                ax[1].scatter(x_drawf,y_drawf,c=color)


                ## Draw filopodia
                x_draw = []
                y_draw = []
                cols = []
                # for i in range(len(xfil)):
                for i in interest:
                    for j in range(len(xfil[i])):

                        start_x = FRONTLOC[i][0]
                        start_y = FRONTLOC[i][1]
                        ax[1].plot([FRONTLOC[i][0], xfil[i][j]], [FRONTLOC[i][1], yfil[i][j]], c='w', lw=2, alpha=0.3)
                        cols.append('w')
                        
                        if(xfil[i][j] < 500):
                            x_draw.append(xfil[i][j]-500+1000)
                            start_x = start_x - 500 +1000
                        elif(xfil[i][j] > 1000):
                            x_draw.append(xfil[i][j]-1000+500)
                            start_x = start_x -1000+500
                        else:
                            x_draw.append(xfil[i][j])

                        if(yfil[i][j] < 500):
                            y_draw.append(yfil[i][j]-500+750)
                            start_y = start_y -500 + 750
                        elif(yfil[i][j] > 750):
                            y_draw.append(yfil[i][j]-750+500)
                            start_y = start_y -750 + 500
                        else:
                            y_draw.append(yfil[i][j])
                        
                        ax[1].plot([start_x, x_draw[-1]], [start_y, y_draw[-1]], c='w', lw=2, alpha=0.3)

                        # x_draw.append(xfil[i][j])
                        # y_draw.append(yfil[i][j])
                        if(i in interest):
                            verfolgx.append(xfil[i][j])
                            verfolgy.append(yfil[i][j])

                ax[1].scatter(x_draw, y_draw, s=5, c=cols)


                # ax[1].scatter(verfolgx,verfolgy,s=3, c='y')
                # ax[1].scatter([x for il in xfil for x in il], [y for il in yfil for y in il], s=5, c='w')
            #ax[1].scatter(xfil[2*162 + 73, :], yfil[2*162 + 73, :], s=5, c='r')
           

            # counter = 0
            # for i in interest:
            #     for j in range(len(xfil[i])):
            #         ax[1].plot([x_drawf[i], x_draw[j]], [y_drawf[i], y_draw[j]])
            #         counter
            # ax[1].scatter(x_draw,y_draw,c='r')

            if np.ndim(FRONTLOCS_TIME) == 3:
                for i in range(len(FRONTLOC)):
                    # ax[1].plot(FRONTLOCS_TIME[:, i, 0],FRONTLOCS_TIME[:, i, 1], 'r-', lw=2)
                # for i in [interest]:
                    x = [FRONTLOCS_TIME[0,i,0]]
                    y = [FRONTLOCS_TIME[0,i,1]]

                    if(x[0] < 500):
                        x[0] = x[0]-500+1000
                    elif(x[0] > 1000):
                        x[0] = x[0]-1000+500

                    if(y[0] < 500):
                        y[0] = y[0]-500+750
                    elif(y[0] > 750):
                        y[0] = y[0]-750+500

                    for j in range(1, len(FRONTLOCS_TIME)):
                        seperate = False
                        newx = []
                        newy = []

                        # if((i in interest)): print("Diff: ", np.abs(FRONTLOCS_TIME[j,i,0]))
                        # if(np.abs(FRONTLOCS_TIME[j,i,0]-FRONTLOCS_TIME[j-1,i,0]) > 10):
                        # if(((FRONTLOCS_TIME[j,i,0] > 1000) and (FRONTLOCS_TIME[j-1,i,0] < 1000)) or ((FRONTLOCS_TIME[j,i,0] < 500) and (FRONTLOCS_TIME[j-1,i,0] > 500))):
                        
                        if(FRONTLOCS_TIME[j,i,0] > 1000): 
                            if(FRONTLOCS_TIME[j-1,i,0] < 1000):
                                x.append(FRONTLOCS_TIME[j,i,0])
                                seperate = True
                                newx.append(FRONTLOCS_TIME[j-1,i,0]-1000+500)
                                newx.append(FRONTLOCS_TIME[j,i,0]-1000+500)
                            else:
                                x.append(FRONTLOCS_TIME[j,i,0]-1000+500)
                                newx.append(FRONTLOCS_TIME[j-1,i,0]-1000+500)
                                newx.append(FRONTLOCS_TIME[j,i,0]-1000+500)
                        elif(FRONTLOCS_TIME[j,i,0] < 500): 
                            if(FRONTLOCS_TIME[j-1,i,0] > 500):
                                x.append(FRONTLOCS_TIME[j,i,0])
                                seperate = True
                                newx.append(FRONTLOCS_TIME[j-1,i,0]+1000-500)
                                newx.append(FRONTLOCS_TIME[j,i,0]+1000-500)
                            else:
                                x.append(FRONTLOCS_TIME[j,i,0]+1000-500)
                                newx.append(FRONTLOCS_TIME[j-1,i,0]+1000-500)
                                newx.append(FRONTLOCS_TIME[j,i,0]+1000-500)
                        elif((FRONTLOCS_TIME[j,i,0] < 1000) and (FRONTLOCS_TIME[j-1,i,0] > 1000)):
                            seperate = True
                            x.append(FRONTLOCS_TIME[j,i,0]+500-1000)
                            newx.append(FRONTLOCS_TIME[j-1,i,0])
                            newx.append(FRONTLOCS_TIME[j,i,0])
                        elif((FRONTLOCS_TIME[j,i,0] > 500) and (FRONTLOCS_TIME[j-1,i,0] < 500)):
                            seperate = True
                            x.append(FRONTLOCS_TIME[j,i,0]+1000-500)
                            newx.append(FRONTLOCS_TIME[j-1,i,0])
                            newx.append(FRONTLOCS_TIME[j,i,0])
                        else:
                            x.append(FRONTLOCS_TIME[j,i,0])
                            newx.append(FRONTLOCS_TIME[j-1,i,0])
                            newx.append(FRONTLOCS_TIME[j,i,0])
                        
                        # print(np.abs(FRONTLOCS_TIME[j,i,1]-FRONTLOCS_TIME[j-1,i,1]))
                        # if(np.abs(FRONTLOCS_TIME[j,i,1]-FRONTLOCS_TIME[j-1,i,1]) > 10):
                        #     # print("Here: ", i)
                        #     if(np.abs(FRONTLOCS_TIME[j-1,i,1]-500) > np.abs(FRONTLOCS_TIME[j-1,i,1]-750)):
                        #         y.append(FRONTLOCS_TIME[j,i,1]+750-500)
                        #         seperate = True
                        #         newy = FRONTLOCS_TIME[j-1,i,1]-750+500
                        #     else:
                        #         y.append(FRONTLOCS_TIME[j,i,1]-750+500)
                        #         newy = FRONTLOCS_TIME[j-1,i,1]+750-500
                        #         seperate = True
                        # else:
                        #     y.append(FRONTLOCS_TIME[j,i,1])
                        #     newy = FRONTLOCS_TIME[j-1,i,1]

                        if(FRONTLOCS_TIME[j,i,1] > 750): 
                            if(FRONTLOCS_TIME[j-1,i,1] < 750):
                                y.append(FRONTLOCS_TIME[j,i,1])
                                seperate = True
                                newy.append(FRONTLOCS_TIME[j-1,i,1]-750+500)
                                newy.append(FRONTLOCS_TIME[j,i,1]-750+500)
                            else:
                                y.append(FRONTLOCS_TIME[j,i,1]-750+500)
                                newy.append(FRONTLOCS_TIME[j-1,i,1]-750+500)
                                newy.append(FRONTLOCS_TIME[j,i,1]-750+500)
                        elif(FRONTLOCS_TIME[j,i,1] < 500): 
                            if(FRONTLOCS_TIME[j-1,i,1] > 500):
                                y.append(FRONTLOCS_TIME[j,i,1])
                                seperate = True
                                newy.append(FRONTLOCS_TIME[j-1,i,1]+750-500)
                                newy.append(FRONTLOCS_TIME[j,i,1]+750-500)
                            else:
                                y.append(FRONTLOCS_TIME[j,i,1]+750-500)
                                newy.append(FRONTLOCS_TIME[j-1,i,1]+750-500)
                                newy.append(FRONTLOCS_TIME[j,i,1]+750-500)
                        elif((FRONTLOCS_TIME[j,i,1] < 750) and (FRONTLOCS_TIME[j-1,i,1] > 750)):
                            seperate = True
                            y.append(FRONTLOCS_TIME[j,i,1]+500-750)
                            newy.append(FRONTLOCS_TIME[j-1,i,1])
                            newy.append(FRONTLOCS_TIME[j,i,1])
                        elif((FRONTLOCS_TIME[j,i,1] > 500) and (FRONTLOCS_TIME[j-1,i,1] < 500)):
                            seperate = True
                            y.append(FRONTLOCS_TIME[j,i,1]+750-500)
                            newy.append(FRONTLOCS_TIME[j-1,i,1])
                            newy.append(FRONTLOCS_TIME[j,i,1])
                        else:
                            y.append(FRONTLOCS_TIME[j,i,1])
                            newy.append(FRONTLOCS_TIME[j-1,i,1])
                            newy.append(FRONTLOCS_TIME[j,i,1])

                        
                        if(seperate):
                            # print("In seperate")
                            ax[1].plot(x,y, 'r-', lw=2, alpha=0.5)

                            x.clear()
                            y.clear()

                            x = newx
                            y = newy

                    ax[1].plot(x,y, 'r-', lw=2, alpha=0.5)
                            

            ####L-Cell drawing
            if not noLcell and species == 'dros':
                for xx, yy in zip(Xint, Yint):
                    # print(xx, yy)
                    circle_temp = Point(xx, yy).buffer(
                        1)  # type(circle)=polygon
                    # ellipse = shapely.affinity.scale(circle_temp, 0.8*lcellsize_int[tind, 0]*CONVFAC/2, 0.8*lcellsize_int[tind, 1]*CONVFAC/2)  # type(ellipse)=polygon
                    ellipse = shapely.affinity.scale(
                        circle_temp, 1.0*lcell_int_large*CONVFAC/2, 1.0*lcell_int_small*CONVFAC/2)  # type(ellipse)=polygon
                    circle = Point(
                        xx - 0.333*v2[0], yy - 0.333*v2[1]).buffer(0.6*lcell_int_small*CONVFAC)
                    uni = circle.union(ellipse)
                    ax[1].add_patch(descartes.PolygonPatch(
                        uni, fc='g', ec='g', alpha=0.2))
                    
                    if(xx < 514):
                        xx = xx -500+1000
                        circle_temp = Point(xx, yy).buffer(
                        1)  # type(circle)=polygon
                    # ellipse = shapely.affinity.scale(circle_temp, 0.8*lcellsize_int[tind, 0]*CONVFAC/2, 0.8*lcellsize_int[tind, 1]*CONVFAC/2)  # type(ellipse)=polygon
                        ellipse = shapely.affinity.scale(
                            circle_temp, 1.0*lcell_int_large*CONVFAC/2, 1.0*lcell_int_small*CONVFAC/2)  # type(ellipse)=polygon
                        circle = Point(
                            xx - 0.333*v2[0], yy - 0.333*v2[1]).buffer(0.6*lcell_int_small*CONVFAC)
                        uni = circle.union(ellipse)
                        ax[1].add_patch(descartes.PolygonPatch(
                            uni, fc='g', ec='g', alpha=0.2))

            #### black circle around start
            for xx, yy in zip(np.ravel(Xrec), np.ravel(Yrec)):
                # if xx not in REC_Xorig and yy not in REC_Yorig:
                ax[1].add_patch(Circle((xx, yy), 0.5*CONVFAC,
                                ls='--', edgecolor='k', facecolor='none'))

            #ax[1].axis(axis_lim)
            ax[1].axis([500, 1000, 500, 750])
            remove_ticks_and_box(ax[1])

            fig.tight_layout()
            inset_axes(ax[1], width='20%', height='20%', loc=2)
            if type(c_st) == int:
                plt.plot(np.arange(20, hours[-1]+1, 0.01), 1-tanh_cust(
                    np.arange(20, hours[-1]+1, 0.01), c_st, sl_st), 'k-', lw=3)
                plt.scatter(tt, 1-tanh_cust(tt, c_st, sl_st),
                            c='r', s=50, zorder=5)
            else:
                for c in c_st:
                    plt.plot(np.arange(20, hours[-1]+1, 0.01), 1-tanh_cust(
                        np.arange(20, hours[-1]+1, 0.01), c, sl_st), 'k-', lw=1)
                    # xval = np.mean(np.sqrt((FRONTLOC[:, 0] - heel_langle_at_p26ocs_x[:, 0])**2 + (FRONTLOC[:, 1] - heel_locs_y[:, 0])**2))/CONVFAC
                    plt.scatter(tt, 1-tanh_cust(tt, c, sl_st),
                                c='r', s=10, zorder=5)
            # plt.plot(np.arange(0, 10, 0.01), 1-tanh_cust(np.arange(0, 10, 0.01)*CONVFAC, c_st, sl_st), 'k-', lw=3)
            # xval = np.mean(np.sqrt((FRONTLOC[:, 0] - heel_locs_x[:, 0])**2 + (FRONTLOC[:, 1] - heel_locs_y[:, 0])**2))/CONVFAC
            # plt.scatter(xval, 1-tanh_cust(xval*CONVFAC, c_st, sl_st), c='r', s=40, zorder=5)
            plt.ylim([-0.1, 1.1])
            plt.xlabel('Time (hour)', color='w')
            # plt.xlabel('Length of axon (um)', color='w')
            plt.ylabel('Stiffness', color='w')
            plt.gca().yaxis.set_label_position("right")
            plt.gca().yaxis.tick_right()
            plt.gca().tick_params(axis='x', colors='w')
            plt.gca().tick_params(axis='y', colors='w')

            inset_axes(ax[1], width='20%', height='20%', loc=1)
            plt.plot(np.arange(20, hours[-1]+1, 0.01), tanh_cust(
                np.arange(20, hours[-1]+1, 0.01), c_magnet, sl_magnet), 'k-', lw=3)
            plt.scatter(tt, tanh_cust(tt, c_magnet, sl_magnet),
                        c='r', s=40, zorder=5)
            plt.ylim([-0.1, 1.1])
            plt.xlabel('Time (hour)', color='w')
            plt.ylabel('L-cell attraction', color='w')
            plt.gca().tick_params(axis='x', colors='w')
            plt.gca().tick_params(axis='y', colors='w')

            filename = str('%04d' % ii + '.png')
            plt.savefig(filename, dpi=100)
            # plt.close(fig)
            
        ##########################################
        # calculate increment of the growth cones
        stiff = 1 - tanh_cust(tt, c_st, sl_st)  # sl_magnet)

        stiffpart = stiff*meanvec_old
        # print("Meanvec_old: ",meanvec_old[interest])
        sp = np.ravel(np.tile(speeds, (len(Xrec[0]), 1)).T)

        # print("Sp: ", sp, "        Stiff: ", stiff)
        if(not gil):
            filopart = (1-stiff)*sp*np.mean(radi*np.exp(1j*angl), axis=1)
        else:
            
            filopart = (1-stiff)*sp*average_vec(angs,lengths)

        # print("Filopart: ", filopart[interest])
        # print("Testfilo: ", )

        magnet_weight = tanh_cust(tt, c_magnet, sl_magnet)

        # check whether front location is in one of the L-cell regions
        inside_goal_area = np.zeros(len(FRONTLOC), dtype='bool')
        if np.ndim(FRONTLOCS_TIME) == 3:
            a, b = lcell_int_large/2*CONVFAC, lcell_int_small/2*CONVFAC
            # inside_goal_area = distance_to_exp(FRONTLOCS_TIME[0, :, :], FRONTLOCS_TIME[-1, :, :], a, b, Xint, Yint, v2inter[tind], v1inter[tind], mode='target_area')
            for kk in range(len(FRONTLOC)):
                for xx, yy in zip(Xint, Yint):
                    check1 = (np.sqrt((FRONTLOC[kk, 0] - xx + 0.333*v1[0])**2 + (FRONTLOC[kk, 1] - yy + 0.333*v1[1])**2) < (1.2*b + CONVFAC/4))
                    check2 = ((FRONTLOC[kk, 0] - xx)**2/(a+CONVFAC/4)**2 + (FRONTLOC[kk, 1] - yy)**2/(b+CONVFAC/4)**2 <= 1**2) #0.8
                    if check1 or check2:
                        inside_goal_area[kk] = True
        
        magnetpart = A_magnet * magnet_weight * np.exp(1j*np.arctan2(goal_loc[:, 1]-FRONTLOC[:, 1], goal_loc[:, 0]-FRONTLOC[:, 0])) #-lcellsize_int[tind, 0]/2*CONVFAC

        meanvec = np.where(inside_goal_area, (1-magnet_weight)*(stiffpart + filopart) + magnetpart, stiffpart + filopart)
        # meanvec = stiffpart + filopart + magnetpart

        # check for hard turns and soften them
        # meanvec = np.where(np.cos(np.angle(meanvec) - np.angle(meanvec_old)) < 0.9396926207859084, (meanvec + meanvec_old/norm) / 2, meanvec)

        # check if axons are close to magnet goal location, if so: stop movement
        meanvec = np.where(np.sqrt((goal_loc[:, 0]-FRONTLOC[:, 0])**2 + (goal_loc[:, 1]-FRONTLOC[:, 1])**2) <= 5, 0+0*1j, meanvec)

        change = np.array([meanvec.real, meanvec.imag]).T
        # print("Meanvec: ", meanvec)

        # if tind > 0:
        #     currmat = np.array([v1inter[tind], v2inter[tind]]).T
        #     prevmat = np.array([v1inter[tind-1], v2inter[tind-1]]).T
        #     change = (currmat@np.linalg.inv(prevmat)@change.T).T

        # print("Change: ", change[interest])
        FRONTLOC += change *dt

        # print("Nach chnage: ", FRONTLOC[interest])
        # print("Vorher: ", lengths[interest])
        # print("Vorher: ", angs[interest])

        if(update):
            update_parameters(FRONTLOC, angs, lengths, xfil, yfil)

        # print("Nachher: ", lengths[interest])

        ###### Hier [500, 1000, 500, 750]

        # teleported = []
        # for i in range(len(FRONTLOC)):
        #     if(FRONTLOC[i][0] < 500):
        #         FRONTLOC[i][0] = 1000 + FRONTLOC[i][0] - 500
        #         teleported.append(i)
        #     elif(FRONTLOC[i][0] > 1000):
        #         FRONTLOC[i][0] = 500 + FRONTLOC[i][0] - 1000
        #         teleported.append(i)
        #     if(FRONTLOC[i][1] < 500):
        #         FRONTLOC[i][1] = 750 + FRONTLOC[i][1] - 500
        #         teleported.append(i)
        #     elif(FRONTLOC[i][1] > 750):
        #         FRONTLOC[i][1] = 500 + FRONTLOC[i][1] - 750
        #         teleported.append(i)

        # for v in teleported:
        #     for i in xfil[v]:
        #         if(i < 500):
        #             i = 1000 + i - 500
        #         elif(i> 1000): 
        #             i = 500 + i - 1000

        # for v in teleported:
        #     for i in yfil[v]:
        #         if(i < 500):
        #             i = 750 + i - 500
        #         elif(i > 750): 
        #             i = 500 + i - 750

        # meanvec = change[:, 0]*dt + 1j*change[:, 1]*dt
        if tt == 20.:
            FRONTLOCS_TIME = np.vstack((FRONTLOCS_TIME, FRONTLOC))
        else:
            FRONTLOCS_TIME = np.vstack((FRONTLOCS_TIME, FRONTLOC.reshape(
                (1, np.shape(FRONTLOC)[0], np.shape(FRONTLOC)[1]))))
        FRONTLOCS_TIME = np.reshape(FRONTLOCS_TIME, (-1, len(FRONTLOC), 2))
        # (FRONTLOCS_TIME[-1, :, 0] - FRONTLOCS_TIME[0, :, 0] + 1j*(FRONTLOCS_TIME[-1, :, 1] - FRONTLOCS_TIME[0, :, 1]))/dt # meanvec.copy()
        meanvec_old = meanvec.copy()

        # print(FRONTLOCS_TIME)

        # how many of the filopodia are in the target area?
        ### Doesn't work with new x_fil, but not necessary for now

        if(not gil):
            indi = (xfil > axis_lim[0])*(xfil < axis_lim[1]) * \
                (yfil > axis_lim[2])*(yfil < axis_lim[3])
            xfil_in, yfil_in = xfil[indi], yfil[indi]
            counts = (density_lcell > 0)[np.array(
                yfil_in, dtype='int'), np.array(xfil_in, dtype='int')]
            outside_fil = 1 - np.sum(counts)/len(counts)
            # print("Outside fill: ", outside_fil)
            outside_fils.append(outside_fil)


        #######for accurate density, filo teleport

        for v in range(len(xfil)):
            for i in range(len(xfil[v])):
                if(xfil[v][i]  < 500):
                    xfil[v][i] = 1000 + xfil[v][i] - 500
                elif(xfil[v][i]> 1000): 
                    xfil[v][i] = 500 + xfil[v][i] - 1000

        for v in range(len(yfil)):
            for i in range(len(yfil[v])):
                if(yfil[v][i] < 500):
                    yfil[v][i] = 750 + yfil[v][i] - 500
                elif(yfil[v][i] > 750): 
                    # print("Y")
                    yfil[v][i] = 500 + yfil[v][i] - 750


        ###########################
        # create new density
        if(not gil):
            hist_fil = np.histogram2d(np.ravel(xfil), np.ravel(yfil),
                                            bins=[np.arange(2*Nx), np.arange(2*Ny)])[0].T

            if density_type != 'constant':
                hist_fil = np.histogram2d(np.ravel(np.linspace(xfil, FRONTLOC[:, 0, None], 50)), np.ravel(np.linspace(yfil, FRONTLOC[:, 1, None], 50)),
                                        bins=[np.arange(2*Nx), np.arange(2*Ny)])[0].T
                kernel = potential_parabola(POS2, np.array(
                    [[Nx], [Ny]]), np.array([15]), 1).reshape((2*Ny-1, 2*Nx-1))
        
                density_fil = fftconvolve(hist_fil, kernel, mode='same')
                density_fil /= np.max(density_fil)
                
                match species:
                    case 'dros':
                        if noLcell:
                            lcell_ell = potential_parabola(
                                POS2, np.c_[Xint, Yint].T, np.array([lcell_int_large, lcell_int_small])/2*CONVFAC, 0)
                            density_lcell = np.sum(lcell_ell, axis=1).reshape(
                                (2*Ny-1, 2*Nx-1))
                        else:
                            lcell_ell = potential_parabola(
                                POS2, np.c_[Xint, Yint].T, np.array([lcell_int_large, lcell_int_small])/2*CONVFAC, 5*np.max(density_fil))
                            lcell_circ = potential_parabola(
                                POS2, np.c_[Xint - 0.333*v2[0], Yint - 0.333*v2[1]].T, 1.2*np.ones(2)*lcell_int_small/2*CONVFAC, 5*np.max(density_fil))
                            density_lcell = np.sum(lcell_ell, axis=1).reshape(
                                (2*Ny-1, 2*Nx-1)) + np.sum(lcell_circ, axis=1).reshape((2*Ny-1, 2*Nx-1))
                    case 'bib':
                        lcell_ell = potential_parabola(
                            POS2, np.c_[Xint, Yint].T, np.array([30, 30]), 10)
                        density_lcell = np.sum(lcell_ell, axis=1).reshape(
                            (2*Ny-1, 2*Nx-1))

            """
            if species == 'dros':
                if noLcell:
                    lcell_ell = potential_parabola(
                        POS2, np.c_[Xint, Yint].T, LCELL_SIZE[0, :]/2*CONVFAC, 0)
                else:
                    lcell_ell = potential_parabola(
                        POS2, np.c_[Xint, Yint].T, LCELL_SIZE[0, :]/2*CONVFAC, 5 * np.max(density_fil))
                    lcell_circ = potential_parabola(
                    POS2, np.c_[Xint - 0.333*v2[0], Yint - 0.333*v2[1]].T, 1.2*np.ones(2)*LCELL_SIZE[0, 1]/2*CONVFAC, 5*np.max(density_heels))
                # print(np.min(density_fil))
                # print(np.max(density_fil))
            elif species == 'bib':
                lcell_ell = potential_parabola(
                    POS2, np.c_[Xint, Yint].T, np.array([30, 30]), 100)
            density_lcell = np.sum(lcell_ell, axis=1).reshape((2*Ny-1, 2*Nx-1))
            """
            
            match density_type:
                case 'full':
                    density = density_heels + density_lcell + density_nogo_heels + density_fil
                case 'fil_attract':
                    if include_nogozones:
                        density = density_fil - density_lcell - density_nogo_heels
                    else:
                        density = density_fil
                case 'fil_repuls':
                    if include_nogozones:
                        density = density_fil + density_lcell + density_nogo_heels
                    else:
                        density = density_fil
        
        
            #### Density teleport 

            # print(density[500:502, 500:502])

            density[751:852] = density[500:601]
            density[399:500] = density[650:751]

            density = np.transpose(density)

            density[1001:1102] = density[500:601]
            density[399:500] = density[980:1081]

            density = np.transpose(density)
            ###########################
            # draw new filopodia from that density

            xfil_old, yfil_old = xfil.copy(), yfil.copy()
            xfil, yfil, allinds = sample_roi(density, FRONTLOC, np.arctan2(FRONTLOCS_TIME[-1, :, 1] - FRONTLOCS_TIME[0, :, 1], FRONTLOCS_TIME[-1, :, 0] -
                                            FRONTLOCS_TIME[0, :, 0]), region, xind, yind, POS, Xrec, Yrec,
                                            np.tile(radii, (len(Xrec[0]), 1)).T,
                                            np.tile(circ_width, (len(Xrec[0]), 1)).T, 
                                            sampfac = sampling_factor)
            if sampling_factor == 1:
                xfil = np.reshape(xfil, (-1, n_fil))
                yfil = np.reshape(yfil, (-1, n_fil))
            else:
                xfil_new = np.reshape(xfil, (-1, int(n_fil/sampling_factor)))
                yfil_new = np.reshape(yfil, (-1, int(n_fil/sampling_factor)))
                ind = np.random.choice(range(n_fil), n_fil - len(xfil_new[0]), replace=False)
                xfil = np.c_[xfil_new, xfil_old[:, ind]]
                yfil = np.c_[yfil_new, yfil_old[:, ind]]

            teleported = []
            for i in range(len(FRONTLOC)):
                if(FRONTLOC[i][0] < 500):
                    FRONTLOC[i][0] = 1000 + FRONTLOC[i][0] - 500
                    teleported.append(i)
                elif(FRONTLOC[i][0] > 1000):
                    FRONTLOC[i][0] = 500 + FRONTLOC[i][0] - 1000
                    teleported.append(i)
                if(FRONTLOC[i][1] < 500):
                    FRONTLOC[i][1] = 750 + FRONTLOC[i][1] - 500
                    teleported.append(i)
                elif(FRONTLOC[i][1] > 750):
                    FRONTLOC[i][1] = 500 + FRONTLOC[i][1] - 750
                    teleported.append(i)
        
            # for v in teleported:
            #     for i in range(len(xfil[v])):
            #         if(xfil[v][i]  < 500):
            #             xfil[v][i] = 1000 + xfil[v][i] - 500
            #         elif(xfil[v][i]> 1000): 
            #             xfil[v][i] = 500 + xfil[v][i] - 1000

            # for v in teleported:
            #     for i in range(len(yfil[v])):
            #         if(yfil[v][i] < 500):
            #             yfil[v][i] = 750 + yfil[v][i] - 500
            #         elif(yfil[v][i] > 750): 
            #             # print("Y")
            #             yfil[v][i] = 500 + yfil[v][i] - 750


            angl = np.arctan2(
                yfil-FRONTLOC[:, 1][:, None], xfil-FRONTLOC[:, 0][:, None])
            radi = np.sqrt((yfil-FRONTLOC[:, 1][:, None])
                        ** 2 + (xfil-FRONTLOC[:, 0][:, None])**2)


            for v in teleported:
                for i in range(len(xfil[v])):
                    if(xfil[v][i]  < 500):
                        xfil[v][i] = 1000 + xfil[v][i] - 500
                    elif(xfil[v][i]> 1000): 
                        xfil[v][i] = 500 + xfil[v][i] - 1000
            print("Vorher: ", yfil[interest])
            for v in teleported:
                for i in range(len(yfil[v])):
                    if(yfil[v][i] < 500):
                        yfil[v][i] = 750 + yfil[v][i] - 500
                    elif(yfil[v][i] > 750): 
                        # print("Y")
                        yfil[v][i] = 500 + yfil[v][i] - 750
            print("Nachher: ", yfil[interest])

            print(xfil[interest])
            print(yfil[interest])
        else:
            xfil, yfil, angs, lengths, inhib, times = sample_gil(FRONTLOC, types, angs, lengths, inhib, times)
        # # change of bundle grid
        # # (1/2) kick
        # vel += acc * dt/2.0
        # # drift
        # bundle_pos += vel * dt
        # # update accelerations
        # acc = getAcc(bundle_pos)
        # # (1/2) kick
        # vel += acc * dt/2.0

        # # shift the heels, fronts, and filopodia
        # bundle_change_x = bundle_pos[:, 0] - Xint
        # bundle_change_y = bundle_pos[:, 1] - Yint
        # print(np.mean(np.abs(bundle_change_x)))
        # Xrec += bundle_change_x[None, :]
        # Yrec += bundle_change_y[None, :]
        # FRONTLOC[:, 0] += np.concatenate([bundle_change_x]*6)
        # FRONTLOC[:, 1] += np.concatenate([bundle_change_y]*6)
        # FRONTLOCS_TIME[:, :, 0] += np.concatenate([bundle_change_x]*6)[None, :]
        # FRONTLOCS_TIME[:, :, 1] += np.concatenate([bundle_change_y]*6)[None, :]
        # xfil += np.concatenate([bundle_change_x]*6)[:, None]
        # yfil += np.concatenate([bundle_change_y]*6)[:, None]

        # Xint, Yint = bundle_pos[:, 0].copy(), bundle_pos[:, 1].copy()
    
    vector = []
    for i in range(len(lengths)):
        vector.append(len(lengths[i]))
    
    # print(np.sort(collect))
    figs2, ax2 = plt.subplots()
    plt.close(fig)
    # ax2.hist(vector, align="left", bins=[0,1,2,3,4,5,6,7,8,9,10])
    ax2.hist(collect, align='left')
    # ax2.scatter(list(range(len(collect))), (collect))
    plt.savefig("barplot.png", dpi=100)

    # check if the axons went to the correct target
    correct_filtered, correct_full = distance_to_exp(
        FRONTLOCS_TIME[0, :, :], FRONTLOCS_TIME[-1, :, :], Xc, Yc, v1, v2, remove, goal_loc)
    return correct_filtered, correct_full, FRONTLOCS_TIME, outside_fils


if __name__ == "__main__":
    # mpl.use('Qt5Agg')  # turn on GUI for plotting

    create_movie = True
    if create_movie:
        mpl.use('Agg')  # do not show figures while creating movie
        
    density_type = 'full' # 'fil_attract' or 'fil_repuls' or 'full' or 'constant'
    include_nogozones = True

    A_magnet = 10
    CONVFAC = 25.2
    if density_type == 'fil_repuls':# or density_type == 'full':
        sampling_factor = 2
    else:
        sampling_factor = 1
    # dt = 1/3/sampling_factor
    dt = 1/30
    Nx, Ny = 1024, 512

    species = 'dros'  # 'dros' or 'bib'
    single_ablation = False
    noLcell = False

    radii = 1.0 * np.array([4.41, 3.60, 5.15, 1.0 * 3.58, 4.48, 4.74])*CONVFAC
    circ_width = 1.0 * np.array([[76.93677250002409, 66.1581562071056, 58.64359788352946, 68.19374152821266],
                                 [67.10341096706019, 63.6030367287868,
                                  64.58480307212197, 64.68382969250712],
                                 [59.99951401507621, 49.70632397768759,
                                  59.76089748808765, 62.63689465660343],
                                 [65.18759340181808, 64.6663385598332,
                                  54.128563634672744, 58.3338061658033],
                                 [73.7705462792746, 73.75319815476855,
                                  70.71260029847282, 67.54024198457094],
                                 [82.29011450458808, 70.63559641125377, 69.51033837975177, 58.19852457593135]]).mean(axis=1)
    # circ_width = 90 * np.ones(6)
    # print("Circ widht: ", circ_width)

    match species:
        case 'bib':
            speeds = 1.4*np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
        case 'dros':
            speeds = np.array([0.093, 0.053, 0.148, 0.09, 0.052, 0.077])
        case _:
            raise ValueError('Species not implemented')

    
    # measured size of heels
    rm_heels = np.array([[3.10231051831059, 2.91492280980449, 2.04300797855269, 2.16799853146643],
                         [3.36613563287586, 3.6474344212202,
                             2.59092160570733, 2.49818435400807],
                         [3.26, 3.45, 2.15, 2.1],
                         [2.5838454966756, 3.26122733330694,
                             2.43382971588824, 1.79884449690457],
                         [2.73134911063699, 3.95170014933712,
                             3.89010177959303, 3.75649529229474],
                         [3.97548396049507, 3.37776310930197, 3.23239134019713, 3.02306813633871]])

    rm_fronts = np.array([[3.39768913, 4.6265753, 4.61201485, 4.17743499],
                          [3.11686657, 3.37643334, 3.71268422, 2.80596312],
                          [3.63754606, 4.92912784, 3.42151975, 3.4512155],
                          [3.01362226, 3.0008909, 3.62965469, 3.61263049],
                          [3.4521001, 3.91277842, 4.38606931, 4.57895422],
                          [3.99243196, 4.05693107, 4.49508072, 4.3103579]]).mean(axis=1)

    # COL = np.array(['c', 'g', 'r', 'y', 'tab:purple', 'tab:orange'])

    # R2: 0.03, R5: 0.04, R3: 0.35
    # speeds = 0.12*np.ones(6)
    # speeds[2] = 0.12 # different speed for R3
    n_fil = 10

    region = 'circ_segment'
    # can be 'angle', 'magnet_strength', 'magnet_timing', 'stiffness', or 'paths'
    investigate = 'angle'  # 'perc_ablate'

    hours = np.array([20, 26, 30, 35, 40])

    LCELL_SIZE = np.array([[3.96, 2.23], [3.96, 2.23], [3.9, 2.28], [4.71, 2.61], [4.57, 2.63]])

    # front distances and angles from Marion's excel sheet (Cell 2015 paper)
    fronts_dist_old = np.array([[2.25454545454545, 2.6247619047619, 3.91962962962963, 4.17666666666667],
                                [1.65, 1.84166666666667,
                                    2.24222222222222, 2.27866666666667],
                                [2.76105263157895, 4.48,
                                    6.64208333333333, 6.92636363636363],
                                [2.1776, 2.64666666666667,
                                    3.77363636363636, 4.29125],
                                [2.26, 2.22333333333333, 2.356, 2.4675],
                                [2.02, 1.99666666666667, 3.43259259259259, 4.74222222222222]])

    front_ang_old = np.array([[28.2566666666667, 27.017619047619, 26.5007407407407, 30.0891666666667],
                              [139.81, 138.655833333333,
                                  141.716111111111, 143.435333333333],
                              [156.339473684211, 156.907333333333,
                                  159.620416666667, 163.584242424242],
                              [179.5208, 180.847272727273,
                                  181.234848484849, 179.0175],
                              [-154.426666666667, -146.302222222222, -
                                  145.667333333333, -144.179166666667],
                              [-26.0233333333333, -27.447619047619, -25.7514814814815, -31.4866666666667]])
    if species == 'dros':
        # load stuff and incorporate the old front locations
        # vec_opt, FRONTLOCS_ALL, HEELLOCS_AVG = load_grid_heels_fronts()
        # # FRONTLOCS_ALL += np.array([512, 256])[None, None, :]
        # # FRONTLOCS_ALL[:, :, 1] = -FRONTLOCS_ALL[:, :, 1]
        # HEELLOCS_AVG += np.array([1024, 512])[None, None, :]
        # FRONTLOCS_ALL = HEELLOCS_AVG + FRONTLOCS_ALL*25.2

        if noLcell:
            vec_opt = np.array(
                [[7.77392528e+01,  8.06643148e+01,  1.40255712e+02, -5.90464557e-02]])
            rot_ang = -np.arctan2(vec_opt[0, 3], vec_opt[0, 2])
            vec_opt[0, ::2], vec_opt[0, 1::2] = rotate_coord(
                vec_opt[0, ::2], vec_opt[0, 1::2], rot_ang)

        else:
            vec_opt = np.array([[8.83445973e+01, 7.67582595e+01, 1.70344600e+02, -6.24174303e+00],
                                [9.33712481e+01, 9.56606427e+01,
                                    1.77259729e+02, 2.26227648e-02],
                                [1.00776263e+02, 9.85328622e+01,
                                    1.98598157e+02, 1.44052931e+00],
                                [9.51118621e+01, 9.68803288e+01, 1.87098297e+02, 2.72226171e+00]])

            # rotate everything to have the grid perfectly horizontal
            rot_ang = []
            for hh in range(4):
                rot_ang.append(-np.arctan2(vec_opt[hh, 3], vec_opt[hh, 2]))
                # FRONTLOCS_ALL[:, hh, 0], FRONTLOCS_ALL[:, hh, 1] = rotate_coord(FRONTLOCS_ALL[:, hh, 0], FRONTLOCS_ALL[:, hh, 1], rot_ang[-1])
                # HEELLOCS_AVG[:, hh, 0], HEELLOCS_AVG[:, hh, 1] = rotate_coord(HEELLOCS_AVG[:, hh, 0], HEELLOCS_AVG[:, hh, 1], rot_ang[-1])
                vec_opt[hh, ::2], vec_opt[hh, 1::2] = rotate_coord(
                    vec_opt[hh, ::2], vec_opt[hh, 1::2], rot_ang[-1])

        startangs = np.pi/180 * \
            np.array([-140.6786, -64.3245, -17.25796667,
                     13.26706/2, 63.2865, 135.0751667])
        # startangs_std = np.array(
        #    [8.758102091, 10.27462811, 5.683716097, 7.21548511, 14.56159809, 5.942082166])
    elif species == 'bib':
        vec_opt = np.array(
            [[100*np.cos(np.pi/3), 100*np.sin(np.pi/3), 100, 0]])
        startangs = np.pi/180 * np.arange(90, -270, -60)

    print("Startangs: ", startangs)
    nrep = 1
    if nrep > 1:
        print('------------------------------')
        print('nrep is larger than 1')
        print('------------------------------')

    if create_movie:
        cent_stiff = 22 # 28 #np.array([25])
        sl_stiff = 0.45  # np.array([0.45])
        std_ang = 0 #15 * np.pi/180
        perc_ablate = 0
        jitter_heels = 0 #0.25
        if noLcell:
            A_magnet = 0
            c_magnet, sl_magnet = np.array([100]), np.array([0.45])
        else:
            c_magnet, sl_magnet = np.array([35]), np.array([0.45])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        # np.save('correct_full_WT_full11_magnetson', correct_full)

        COL = np.array(['c', 'g', 'r', 'y', 'tab:purple', 'tab:orange'])
        plt.figure(figsize=(8, 6))
        plt.bar(range(6), np.ones(6), width=0.8, facecolor='w', edgecolor=COL)
        plt.bar(range(6), np.mean(correct_full, axis=1), width=0.8, color=COL)
        plt.xticks(range(6), ['R'+str(i+1) for i in range(6)], fontsize=20)
        plt.ylabel('Wiring performance', fontsize=20)
        plt.yticks(fontsize=20)

        # COL = np.array(['c', 'g', 'r', 'y', 'tab:purple', 'tab:orange'])
        # plt.figure(figsize=(8, 6))
        # plt.bar(range(6), np.ones(6), width=0.8, facecolor='w', edgecolor=COL)
        # plt.bar(range(6), np.ones(6), width=0.8, color=COL)
        # plt.xticks(range(6), ['R'+str(i+1) for i in range(6)], fontsize=20)
        # plt.ylabel('Wiring performance', fontsize=20)
        # plt.yticks(fontsize=20)


        """
        #######
        density_type = 'full'
        include_nogozones = False
        cent_stiff = 0 # 28 #np.array([25])
        sl_stiff = 0.45  # np.array([0.45])
        std_ang = 0 #15 * np.pi/180
        perc_ablate = 0
        jitter_heels = 0
        c_magnet, sl_magnet = np.array([100]), np.array([0.45])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_full00', correct_full)
        
        include_nogozones = True
        cent_stiff = 0 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_full10', correct_full)

        include_nogozones = False
        cent_stiff = 22 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_full01', correct_full)

        include_nogozones = True
        cent_stiff = 22 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_full11', correct_full)
        
        ##############
        density_type = 'constant'
        include_nogozones = False
        cent_stiff = 0  # 28 #np.array([25])
        sl_stiff = 0.45  # np.array([0.45])
        std_ang = 0 #15 * np.pi/180
        jitter_heels = 0
        c_magnet, sl_magnet = np.array([100]), np.array([0.45])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_const00', correct_full)
        
        include_nogozones = True
        cent_stiff = 0 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_const10', correct_full)

        include_nogozones = False
        cent_stiff = 22 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_const01', correct_full)

        include_nogozones = True
        cent_stiff = 22 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_const11', correct_full)

        #########
        density_type = 'fil_repuls'
        include_nogozones = False
        cent_stiff = 0  # 28 #np.array([25])
        sl_stiff = 0.45  # np.array([0.45])
        std_ang = 0 #15 * np.pi/180
        jitter_heels = 0
        c_magnet, sl_magnet = np.array([100]), np.array([0.45])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_filrep00', correct_full)
        
        include_nogozones = True
        cent_stiff = 0 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_filrep10', correct_full)

        include_nogozones = False
        cent_stiff = 22 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_filrep01', correct_full)

        include_nogozones = True
        cent_stiff = 22 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_filrep11', correct_full)

        #########
        density_type = 'fil_attract'
        include_nogozones = False
        cent_stiff = 0  # 28 #np.array([25])
        sl_stiff = 0.45  # np.array([0.45])
        std_ang = 0 #15 * np.pi/180
        jitter_heels = 0
        c_magnet, sl_magnet = np.array([100]), np.array([0.45])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_filatt00', correct_full)
        
        include_nogozones = True
        cent_stiff = 0 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_filatt10', correct_full)

        include_nogozones = False
        cent_stiff = 22 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_filatt01', correct_full)

        include_nogozones = True
        cent_stiff = 22 # 28 #np.array([25])
        correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
            cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
        np.save('correct_full_noLcell_filatt11', correct_full)
        """
        # startpos = FRONTLOCS_TIME[0, :, :]
        # keep = np.where((startpos[:, 0] > 160) * (startpos[:, 0] < 1850)
        #                 * (startpos[:, 1] > 120) * (startpos[:, 1] < 925), True, False)

        # vec = np.repeat(
        #     np.array([CONVFAC * np.cos(startangs), np.sin(startangs)]).T, 162, axis=0)
        # proj = np.sum((FRONTLOCS_TIME[1:, :, :] - startpos[None, :, :]) /
        #               CONVFAC * vec[None, :, :], axis=2) / np.sqrt(np.sum(vec*vec, axis=1))

        # proj = proj[:, keep]
        # angles = 180/np.pi * np.arctan2(FRONTLOCS_TIME[1:, :, 1] - FRONTLOCS_TIME[:-1, :, 1],
        #                                 FRONTLOCS_TIME[1:, :, 0] - FRONTLOCS_TIME[:-1, :, 0])[:, keep]
        # rec = np.repeat(range(6), 162)[keep]

        # meanang = np.zeros((6, len(angles)))
        # stdang = np.zeros((6, len(angles)))
        # angu = np.zeros((6, len(angles)))
        # angd = np.zeros((6, len(angles)))
        # for i in range(6):
        #     angles2 = angles[:, rec == i]
        #     if i == 0 or i == 5:
        #         angles2 = np.where(angles2 < 0, angles2 + 360, angles2)
        #     meanang[i] = np.mean(angles2, axis=1)
        #     stdang[i] = np.std(angles2, axis=1)
        #     angu[i] = np.percentile(angles2, 75, axis=1)
        #     angd[i] = np.percentile(angles2, 25, axis=1)

        # COL = np.array(['c', 'g', 'r', 'y', 'tab:purple', 'tab:orange'])
        # plt.figure()
        # time_sim = np.linspace(20+dt, 35+dt, len(angles))
        # for i in range(6):
        #     if i == 0:
        #         meanang[i, :] -= 360
        #         angu[i, :] -= 360
        #         angd[i, :] -= 360
        #     plt.plot(time_sim, meanang[i, :], c=COL[i], lw=2)
        #     # plt.plot(time_sim, angu[i, :], c=COL[i])
        #     # plt.plot(time_sim, angd[i, :], c=COL[i])
        #     # plt.fill_between(time_sim, angu[i, :], angd[i, :], color=COL[i], alpha = 0.3)
        #     plt.fill_between(time_sim, meanang[i, :]-stdang[i, :],
        #                      meanang[i, :] + stdang[i, :], color=COL[i], alpha=0.3)
        # plt.xlabel('Developmental time (hAPF)')
        # plt.ylabel('Instantaneous gc angle (deg)')

        # plt.figure()
        # for i in range(6):
        #     plt.plot(time_sim, angu[i, :] - angd[i, :], c=COL[i], lw=2)
        # plt.xlabel('Developmental time (hAPF)')
        # plt.ylabel('IQR of gc angle (deg)')

        # plt.figure()
        # for i in range(6):
        #     plt.plot((FRONTLOCS_TIME[:, keep, 0][:, rec == i] - startpos[keep, 0][rec == i])/CONVFAC,
        #              (FRONTLOCS_TIME[:, keep, 1][:, rec == i] - startpos[keep, 1][rec == i])/CONVFAC, c=COL[i], lw=2, alpha=0.5)
        # plt.xlabel('x [um]')
        # plt.ylabel('y [um]')

        # FRONTLOCS_TIME_const = np.load('frontlocs_time_constantdensity.npy')
        # FRONTLOCS_TIME_constplus = np.load(
        #     'frontlocs_time_const+heels+lcells.npy')
        # FRONTLOCS_TIME_attractfil = np.load('frontlocs_time_attractfil.npy')
        # FRONTLOCS_TIME_attractfilstiff = np.load(
        #     'frontlocs_time_attractfil+stiff23.npy')
        # FRONTLOCS_TIME_repulsfil = np.load('frontlocs_time_repulsfil.npy')
        # FRONTLOCS_TIME_repulsfilstiff = np.load(
        #     'frontlocs_time_repulsfil+stiff23.npy')
        # FRONTLOCS_TIME_fulldens = np.load('frontlocs_time_fulldensity.npy')
        # FRONTLOCS_TIME_fulldensstiff = np.load(
        #     'frontlocs_time_fulldensity+stiff23.npy')

        # vec = np.repeat(
        #     np.array([CONVFAC * np.cos(startangs), np.sin(startangs)]).T, 162, axis=0)
        # proj_const = np.sum((FRONTLOCS_TIME_const[1:, :, :] - startpos[None, :, :]) /
        #                     CONVFAC * vec[None, :, :], axis=2) / np.sqrt(np.sum(vec*vec, axis=1))
        # proj_constplus = np.sum((FRONTLOCS_TIME_constplus[1:, :, :] - startpos[None, :, :]) /
        #                         CONVFAC * vec[None, :, :], axis=2) / np.sqrt(np.sum(vec*vec, axis=1))
        # proj_attractfil = np.sum((FRONTLOCS_TIME_attractfil[1:, :, :] - startpos[None, :, :]) /
        #                          CONVFAC * vec[None, :, :], axis=2) / np.sqrt(np.sum(vec*vec, axis=1))
        # proj_attractfilstiff = np.sum(
        #     (FRONTLOCS_TIME_attractfilstiff[1:, :, :] - startpos[None, :, :])/CONVFAC * vec[None, :, :], axis=2) / np.sqrt(np.sum(vec*vec, axis=1))
        # proj_repulsfil = np.sum((FRONTLOCS_TIME_repulsfil[1:, :, :] - startpos[None, :, :]) /
        #                         CONVFAC * vec[None, :, :], axis=2) / np.sqrt(np.sum(vec*vec, axis=1))
        # proj_repulsfilstiff = np.sum(
        #     (FRONTLOCS_TIME_repulsfilstiff[1:, :, :] - startpos[None, :, :])/CONVFAC * vec[None, :, :], axis=2) / np.sqrt(np.sum(vec*vec, axis=1))
        # proj_fulldens = np.sum((FRONTLOCS_TIME_fulldens[1:, :, :] - startpos[None, :, :]) /
        #                        CONVFAC * vec[None, :, :], axis=2) / np.sqrt(np.sum(vec*vec, axis=1))
        # proj_fulldensstiff = np.sum(
        #     (FRONTLOCS_TIME_fulldensstiff[1:, :, :] - startpos[None, :, :])/CONVFAC * vec[None, :, :], axis=2) / np.sqrt(np.sum(vec*vec, axis=1))

        # plt.figure()
        # for i in [0]:  # range(6):
        #     plt.plot(time_sim, np.std(
        #         proj_const[:, keep][:, rec == i], axis=1), c=COL[i], lw=2)
        #     plt.plot(time_sim, np.std(
        #         proj_constplus[:, keep][:, rec == i], axis=1), '--', c=COL[i], lw=2)
        #     plt.plot(time_sim, np.std(
        #         proj_attractfil[:, keep][:, rec == i], axis=1), '+', c=COL[i], lw=2)
        #     plt.plot(time_sim, np.std(
        #         proj_attractfilstiff[:, keep][:, rec == i], axis=1), 'd', c=COL[i], lw=2)
        #     plt.plot(time_sim, np.std(
        #         proj_repulsfil[:, keep][:, rec == i], axis=1), 's', c=COL[i], lw=2)
        #     plt.plot(time_sim, np.std(
        #         proj_repulsfilstiff[:, keep][:, rec == i], axis=1), 'h', c=COL[i], lw=2)
        #     plt.plot(time_sim, np.std(
        #         proj_fulldens[:, keep][:, rec == i], axis=1), 'o', c=COL[i], lw=2)
        #     plt.plot(time_sim, np.std(
        #         proj_fulldensstiff[:, keep][:, rec == i], axis=1), '1', c=COL[i], lw=2)

        # plt.xlabel('Developmental time (hAPF)')
        # plt.ylabel('Standard-deviation of paths (um)')
        # plt.legend(['Constant density', 'Constant density + no-go zones', 'Attractive filopodia',
        #             'Attractive filopodia with stiffness', 'Repulsive filopodia',
        #             'Repulsive filopodia with stiffness', 'Full density', 'Full density + stiffness'])

        ##############################################################
        # np.save('frontlocs_time_constantdensity', FRONTLOCS_TIME)
        # np.save('frontlocs_time_const+heels+lcells', FRONTLOCS_TIME)
        # np.save('frontlocs_time_attractfil', FRONTLOCS_TIME)
        # np.save('frontlocs_time_attractfil+stiff23', FRONTLOCS_TIME)
        # np.save('frontlocs_time_repulsfil+stiff23', FRONTLOCS_TIME)
        # np.save('frontlocs_time_fulldensity', FRONTLOCS_TIME)
        # np.save('frontlocs_time_fulldensity+stiff23', FRONTLOCS_TIME)
        # FRONTLOCS_TIME = np.load('frontlocs_time_constantdensity.npy')
        # np.save('len_avg_R'+str(receptor), np.array(len_avg_all).reshape((nrep, -1)))
        # print(str(np.round(np.mean(np.array(outside_fil_all)) * 100, 1)) + '% of filopodial tips are in the target area')
    else:
        if investigate == 'heel_jitter':
            print('heel_jitter')
            cent_stiff = 22
            sl_stiff = 0.45
            std_ang = 0
            perc_ablate = 0
            c_magnet, sl_magnet = np.array([100]), np.array([0.45])
            jitter = np.linspace(0, 0.4, 11)
            corrects = []
            lengths = []
            all_corrects = []
            for jitter_heels in jitter:
                print(jitter_heels)
                correct_filtered, correct_full, frontlocs, out = run_gc(
                    cent_stiff, sl_stiff, 0, perc_ablate, c_magnet, sl_magnet)
                corrects.append(np.mean(correct_filtered))
                lengths.append(np.size(correct_filtered))
                all_corrects.append(correct_full)

            np.save('jitter_heeljitter', jitter)
            np.save('correct_filtered_heeljitter', corrects)
            np.save('correct_full_heeljitter', all_corrects)
            np.save('lengths_heeljitter', lengths)

            #####################################
            jitter = np.load('jitter_heeljitter.npy')
            corrects = np.load('correct_filtered_heeljitter.npy')
            all_corrects = np.load('correct_full_heeljitter.npy')
            lengths = np.load('lengths_heeljitter.npy')
            
            jitter_const = np.load('jitter_heeljitter_constant11.npy')
            corrects_const = np.load('correct_filtered_heeljitter_constant11.npy')
            all_corrects_const = np.load('correct_full_heeljitter_constant11.npy')
            lengths_const = np.load('lengths_heeljitter_constant11.npy')
            

            plt.figure()
            plt.plot(jitter, corrects, 'ko--')
            plt.plot(jitter_const, corrects_const, 'bo--')
            plt.ylim([0, 1.05])
            plt.xlabel('Heel jitter (um)')
            plt.ylabel('Wiring performance')
            plt.grid()
            # for p, c, l in zip(jitter, corrects, lengths):
            #     plt.text(p, c - 0.05, 'n='+str(l),
            #               horizontalalignment='center', fontsize=8)

            # subtype specific performance curves
            perf_subt = np.zeros((len(jitter), 6))
            perf_subt_cont = np.zeros((len(jitter), 6))
            for i, ent in enumerate(all_corrects):
                for j, line in enumerate(ent):
                    perf_subt[i, j] = np.mean(line[line != -1])
            for i, ent in enumerate(all_corrects_const):
                for j, line in enumerate(ent):
                    perf_subt_cont[i, j] = np.mean(line[line != -1])

            COL = np.array(['c', 'g', 'r', 'y', 'tab:purple', 'tab:orange'])
            plt.figure()
            for i in range(6):
                plt.plot(jitter, perf_subt[:, i], 'o-', lw=2, c=COL[i])
                plt.plot(jitter_const, perf_subt_cont[:, i], 'o--', c=COL[i])
            plt.ylim([0, 1.05])
            plt.xlabel('Heel jitter (um)')
            plt.ylabel('Wiring performance')
            plt.grid()

        elif investigate == 'perc_ablate':
            print('perc_ablate')
            cent_stiff = 22
            sl_stiff = 0.45
            std_ang = 0
            jitter_heels = 0
            c_magnet, sl_magnet = np.array([100]), np.array([0.45])
            percs = np.arange(0, 51, 5)
            corrects = []
            lengths = []
            all_corrects = []
            for p in percs:
                print(p)
                correct_filtered, correct_full, frontlocs, out = run_gc(
                    cent_stiff, sl_stiff, 0, p, c_magnet, sl_magnet)
                corrects.append(np.mean(correct_filtered))
                lengths.append(np.size(correct_filtered))
                all_corrects.append(correct_full)

            np.save('percs_percablate', percs)
            np.save('correct_filtered_percablate', corrects)
            np.save('correct_full_percablate', all_corrects)
            np.save('lengths_percablate', lengths)

            #####################################
            percs = np.load('percs_percablate.npy')
            corrects = np.load('correct_filtered_percablate.npy')
            all_corrects = np.load('correct_full_percablate.npy')
            lengths = np.load('lengths_percablate.npy')
            
            percs_const = np.load('percs_percablate_constant11.npy')
            corrects_const = np.load('correct_filtered_percablate_constant11.npy')
            all_corrects_const = np.load('correct_full_percablate_constant11.npy')
            lengths_const = np.load('lengths_percablate_constant11.npy')
            

            plt.figure()
            plt.plot(percs, corrects, 'ko--')
            plt.plot(percs_const, corrects_const, 'bo--')
            plt.ylim([0, 1.05])
            plt.xlabel('Ablation percentage')
            plt.ylabel('Wiring performance')
            plt.grid()
            # for p, c, l in zip(jitter, corrects, lengths):
            #     plt.text(p, c - 0.05, 'n='+str(l),
            #               horizontalalignment='center', fontsize=8)

            # subtype specific performance curves
            perf_subt = np.zeros((len(percs), 6))
            perf_subt_cont = np.zeros((len(percs), 6))
            for i, ent in enumerate(all_corrects):
                for j, line in enumerate(ent):
                    perf_subt[i, j] = np.mean(line[line != -1])
            for i, ent in enumerate(all_corrects_const):
                for j, line in enumerate(ent):
                    perf_subt_cont[i, j] = np.mean(line[line != -1])

            COL = np.array(['c', 'g', 'r', 'y', 'tab:purple', 'tab:orange'])
            plt.figure()
            for i in range(6):
                #plt.plot(percs, perf_subt[:, i], 'o-', lw=2, c=COL[i])
                plt.plot(percs_const, perf_subt_cont[:, i], 'o--', c=COL[i])
            plt.ylim([0, 1.05])
            plt.xlabel('Ablation percentage')
            plt.ylabel('Wiring performance')
            plt.grid()


        elif investigate == 'angle':
            print('angle')
            cent_stiff = 0  # 22
            sl_stiff = 0.45
            c_magnet, sl_magnet = 100, 0.45
            perc_ablate = 0
            jitter_heels = 0

            samp = 11

            std_angs = np.linspace(0, 30, samp) * np.pi/180
            # data = Parallel(n_jobs=-1, verbose=100)(delayed(run_gc)(cent_stiff, sl_stiff,
            #                                                         startang, A_magnet, circ_width, r) for startang in np.repeat(startangs, rep))

            corrmat = np.zeros(samp)
            for i, std_ang in enumerate(std_angs):
                correct_filtered, correct_full, FRONTLOCS_TIME, outside_fil = run_gc(
                    cent_stiff, sl_stiff, A_magnet, perc_ablate, c_magnet, sl_magnet)
                corrmat[i] = np.mean(correct_filtered)
                np.save('corr_full_angles_' + density_type +'_'+str(int(include_nogozones))+str(int(cent_stiff>0)) + '_' +str(i), correct_full)

            np.save('angles_sweep_filokde_'+density_type+'_'+str(int(include_nogozones))+str(int(cent_stiff>0)), corrmat)

            # # load and plot
            cor = np.load(
                '/home/eric/axon_guidance/angles_sweep_filokde_constant_00.npy')
            plt.figure()
            plt.plot(std_angs/np.pi*180, cor, 'o-')
            plt.xlabel('Angular jitter (deg)')
            plt.ylabel('Wiring performance')
            plt.ylim([-0.05, 1.05])
            plt.grid()

            cor2 = np.zeros((6, samp))
            cor2_nogo = np.zeros((6, samp))
            for i in range(samp):
                fullcorr = np.load(
                    '/home/eric/axon_guidance/corr_full_angles_'+density_type+'_'+str(int(False))+str(int(cent_stiff>0)) + '_' +str(i) + '.npy')
                cor2[:, i] = np.mean(fullcorr, axis=1)
                fullcorr_nogo = np.load(
                    '/home/eric/axon_guidance/corr_full_angles_'+density_type+'_'+str(int(True))+str(int(cent_stiff>0)) + '_' +str(i) + '.npy')
                cor2_nogo[:, i] = np.mean(fullcorr_nogo, axis=1)
                
            COL = np.array(['c', 'g', 'r', 'y', 'tab:purple', 'tab:orange'])
            plt.figure()
            for i in range(6):
                plt.plot(std_angs/np.pi*180, cor2[i, :], 'o-', c=COL[i])
                plt.plot(std_angs/np.pi*180, cor2_nogo[i, :], '+-', c=COL[i])
            plt.xlabel('Angular jitter (deg)')
            plt.ylabel('Wiring performance')
            plt.ylim([-0.05, 1.05])
            plt.grid()

        elif investigate == 'stiffness':
            print('stiffness')
            # 2D parameter sweep: central point and slope
            # nsamp = 11
            # cent_stiff = np.linspace(25, 32, nsamp)
            # sl_stiff = np.linspace(0.1, 2, nsamp)
            # std_ang = 0

            # corrmat = np.zeros((nsamp, nsamp))
            # for i, c_st in enumerate(cent_stiff):
            #     for j, sl_st in enumerate(sl_stiff):
            #         cor, cor2, fr, out = run_gc(c_st, sl_st, A_magnet, 0)
            #         corrmat[i, j] = np.mean(cor)
            #         np.save('corr_full_'+str(i)+'_'+str(j), cor2)

            # np.save('stiffness_sweep_filokde_WT', corrmat)

            # 1D parameter sweep: fixed slope
            nsamp = 11
            cent_stiff = np.linspace(20, 32, nsamp)
            sl_stiff = 0.45
            std_ang = 0

            corrmat = np.zeros(nsamp)
            for i, c_st in enumerate(cent_stiff):
                cor, cor2, fr, out = run_gc(c_st, sl_stiff, A_magnet, 0)
                corrmat[i] = np.mean(cor)
                np.save('corr_full_1dstiff_'+str(i), cor2)

            np.save('stiffness_1dsweep_filokde_WT', corrmat)

            # 2D load and plot
            nsamp = 11
            cent_stiff = np.linspace(25, 32, nsamp)
            sl_stiff = np.linspace(0.1, 2, nsamp)
            cor2 = np.zeros((6, nsamp, nsamp))
            for i in range(nsamp):
                for j in range(nsamp):
                    try:
                        fullcorr = np.load(
                            '/home/eric/axon_guidance/remote_computing/corr_full_'+str(i)+'_'+str(j)+'.npy')
                        cor2[:, i, j] = np.mean(fullcorr, axis=1)
                    except:
                        None

            fig35, ax35 = plt.subplots(1, 1, figsize=(8, 8))
            heatmap(cor2[5, :, :], np.round(cent_stiff, 1), np.where(sl_stiff < 1, np.round(sl_stiff, 2), np.round(sl_stiff, 2)), extralines=False,
                    cbarlabel='Mean performance', xlabel='Slope of stiffness', ylabel='Central point (hAPF)', ax1=ax35, vmin=0, vmax=1, origin='lower')
            plt.tight_layout()

            # 1D load and plot
            cor2 = np.zeros((6, nsamp))
            for i in range(nsamp):
                try:
                    fullcorr = np.load(
                        '/home/eric/axon_guidance/remote_computing/corr_full_1dstiff_'+str(i)+'.npy')
                    cor2[:, i] = np.mean(fullcorr, axis=1)
                except:
                    None

            COL = np.array(['c', 'g', 'r', 'y', 'tab:purple', 'tab:orange'])
            plt.figure()
            for i in range(6):
                plt.plot(cent_stiff, cor2[i, :], 'o-', c=COL[i])
            plt.xlabel('Stiffness transition (hAPF)')
            plt.ylabel('Wiring performance')
            plt.ylim([-0.05, 1.05])
            plt.grid()
            plt.title('Slope = 0.45/hour')

        elif investigate == 'magnet_strength':
            cent_stiff = np.array([25])
            sl_stiff = np.array([0.4])
            # cent_stiff = 25.2*np.array([5])
            # sl_stiff = np.array([0.02])
            samp = 11
            As = np.linspace(0, 10, samp)

            novalleydyn = False
            rep = 1
            data = Parallel(n_jobs=-1, verbose=100)(delayed(run_gc)(cent_stiff, sl_stiff,
                                                                    startang, A_magnet, circ_width, r) for A_magnet in np.repeat(As, rep))
            cumerr = np.array([d[1] for d in data])
            cumerr = np.reshape(np.array(cumerr), (samp, rep))
            err = np.mean(cumerr, axis=1)/nr_of_rec

            # novalleydyn = True
            # rep = 1
            # data = Parallel(n_jobs=-1, verbose=100)(delayed(run_gc)(cent_stiff, sl_stiff, startang, A_magnet) for A_magnet in np.repeat(As, rep))
            # cumerr_noval = np.array([d[1] for d in data])
            # cumerr_noval = np.reshape(np.array(cumerr_noval), (samp, rep))
            # err_noval = np.mean(cumerr_noval, axis=1)/nr_of_rec

            plt.figure()
            plt.plot(As, err, lw=3, label='Full model')
            # plt.plot(As, err_noval, lw=3, label='No valley dynamics')
            # plt.errorbar(As, err, yerr=np.std(cumerr/nr_of_rec, axis=1)/np.sqrt(rep), lw=3)
            # plt.errorbar(As, err_noval, yerr=np.std(cumerr_noval/nr_of_rec, axis=1)/np.sqrt(rep), lw=3)
            plt.xlabel('R' + str(receptor)+' magnet strength')
            plt.ylabel('Mean performance')
            # plt.legend()
            plt.ylim([0, 1.05])

        elif investigate == 'magnet_timing':
            nsamp = 10
            cent_magnet = np.linspace(25, 38, nsamp)
            slopes_magnet = np.linspace(0.1, 2, nsamp)

            # c_st = 1.0
            # sl_st = 0.02
            c_st = 25
            sl_st = 0.45
            data = Parallel(n_jobs=-1, verbose=100)(delayed(run_gc)(c_st, sl_st, startang, A_magnet, circ_width,
                                                                    r, c_magnet, sl_magnet) for c_magnet in cent_magnet for sl_magnet in slopes_magnet)
            cumerr_35 = np.array([d[0] for d in data])
            cumerr_35 = np.reshape(np.array(cumerr_35), (nsamp, nsamp))
            cumerr_40 = np.array([d[1] for d in data])
            cumerr_40 = np.reshape(np.array(cumerr_40), (nsamp, nsamp))
            cumerr_35_v = np.array([d[2] for d in data])
            cumerr_35_v = np.reshape(np.array(cumerr_35_v), (nsamp, nsamp))
            cumerr_40_v = np.array([d[3] for d in data])
            cumerr_40_v = np.reshape(np.array(cumerr_40_v), (nsamp, nsamp))
            len_at_times = np.array([d[4] for d in data])
            len_at_times = np.mean(len_at_times, axis=0)

            fig_35, ax_35 = plt.subplots(1, 1, figsize=(10, 6))
            # ax_h[0].set_xlim(0, 1.2)
            heatmap(cumerr_35/nr_of_rec, np.round(cent_magnet, 1), np.where(slopes_magnet < 2, np.round(slopes_magnet, 3), np.array(slopes_magnet, dtype='int')),
                    extralines=False, cbarlabel='Mean performance', xlabel='Slope', ylabel='Central point (hAPF)', ax1=ax_35, vmin=0, vmax=1, origin='lower')
            plt.tight_layout()

            fig_40, ax_40 = plt.subplots(1, 1, figsize=(10, 6))
            # ax_h[0].set_xlim(0, 1.2)
            heatmap(cumerr_40/nr_of_rec, np.round(cent_magnet, 1), np.where(slopes_magnet < 2, np.round(slopes_magnet, 3), np.array(slopes_magnet, dtype='int')),
                    extralines=False, cbarlabel='Mean performance', xlabel='Slope', ylabel='Central point (hAPF)', ax1=ax_40, vmin=0, vmax=1, origin='lower')
            plt.tight_layout()

            fig_35_v, ax_35_v = plt.subplots(1, 1, figsize=(10, 6))
            # ax_h[0].set_xlim(0, 1.2)
            heatmap(cumerr_35_v/nr_of_rec, np.round(cent_magnet, 1), np.where(slopes_magnet < 2, np.round(slopes_magnet, 3), np.array(slopes_magnet, dtype='int')),
                    extralines=False, cbarlabel='Mean performance', xlabel='Slope', ylabel='Central point (hAPF)', ax1=ax_35_v, vmin=0, vmax=1, origin='lower')
            plt.tight_layout()

            fig_40_v, ax_40_v = plt.subplots(1, 1, figsize=(10, 6))
            # ax_h[0].set_xlim(0, 1.2)
            heatmap(cumerr_40_v/nr_of_rec, np.round(cent_magnet, 1), np.where(slopes_magnet < 2, np.round(slopes_magnet, 3), np.array(slopes_magnet, dtype='int')),
                    extralines=False, cbarlabel='Mean performance', xlabel='Slope', ylabel='Central point (hAPF)', ax1=ax_40_v, vmin=0, vmax=1, origin='lower')
            plt.tight_layout()

        elif investigate == 'paths':
            include_equator = False
            c_st = 25
            sl_st = 0.45
            nsamp = 16
            nrep = 5
            # circ_widths_all = circ_width * np.linspace(0.5, 2, nsamp)
            # data = Parallel(n_jobs=-1, verbose=100)(delayed(run_gc)(c_st, sl_st, startang, A_magnet, circw) for circw in np.repeat(circ_widths_all, nrep))
            radius_all = r * np.linspace(0.5, 2, nsamp)
            data = Parallel(n_jobs=-1, verbose=100)(delayed(run_gc)(c_st, sl_st, startang,
                                                                    A_magnet, circ_width, radius) for radius in np.repeat(radius_all, nrep))

            FRONTLOCS_TIME = np.array([d[0] for d in data]).reshape(
                (nsamp, nrep, -1, 1, 2))
            paths = np.mean(FRONTLOCS_TIME, axis=1)[:, :, 0, :]
            REC_X = np.array([d[1] for d in data])[0, :, :]
            REC_Y = np.array([d[2] for d in data])[0, :, :]
            Xint = np.array([d[3] for d in data])[0, :]
            Yint = np.array([d[4] for d in data])[0, :]

            plt.figure()
            cmap = plt.get_cmap('Purples', nsamp)
            for i in range(nsamp):
                plt.plot(paths[i, :, 0], paths[i, :, 1], lw=2, c=cmap(i))
            norm = mpl.colors.Normalize(vmin=0.5, vmax=2)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cb = plt.colorbar(sm, ticks=np.linspace(0.5, 2, nsamp))
            cb.ax.tick_params(labelsize=10)
            # cb.set_label('Relative size of ROI')
            cb.set_label('Relative radius of ROI')

            plt.scatter(REC_X, REC_Y)

            for xx, yy in zip(Xint, Yint):
                circle_temp = Point(xx, yy).buffer(1)  # type(circle)=polygon
                # ellipse = shapely.affinity.scale(circle_temp, 0.8*lcellsize_int[tind, 0]*CONVFAC/2, 0.8*lcellsize_int[tind, 1]*CONVFAC/2)  # type(ellipse)=polygon
                ellipse = shapely.affinity.scale(
                    circle_temp, 1.0*LCELL_SIZE[-1, 0]*CONVFAC/2, 1.0*LCELL_SIZE[-1, 1]*CONVFAC/2)  # type(ellipse)=polygon
                circle = Point(
                    xx - 0.333*vec_opt[-1, 2], yy).buffer(0.6*LCELL_SIZE[-1, 1]*CONVFAC)
                uni = circle.union(ellipse)
                plt.gca().add_patch(descartes.PolygonPatch(uni, fc='g', ec='g', alpha=0.2))

            plt.axis([1000, 1300, 400, 600])

    if create_movie:
        command = ('mencoder',
                   'mf://00*.png',
                   '-mf',
                   'type=png:w=800:h=600:fps=3',
                   '-ovc',
                   'lavc',
                   '-lavcopts',
                   'vcodec=mpeg4:vbitrate=4000',
                   '-oac',
                   'copy',
                   '-o',
                   'front_movement_model.avi')

        print('Starting to stitch images together')
        os.spawnvp(os.P_WAIT, 'mencoder', command)
        # for file in glob.glob("*.png"):
        #     os.remove(file)
        mpl.use('Qt5Agg')  # turn on GUI for plotting

    print('\nTotal time: ' + str(np.round(time.time()-start, 2)) + ' seconds')

#%%
# fitting a binomial GLM to the results, to identify the most influential parameters for good wiring performance
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm
import xarray as xr

xarray_3d = xr.Dataset(
    {"R1": (("model", "test"), np.array([[0.96, 0.57, 1.0, 1.0, 0.99, 0.75, 1.0, 1.0, 0.33, 0.45, 0.88, 0.96, 0.52, 0.51, 1.0, 1.0],
                                         [0.92, 0.43, 0.97, 0.89, 0.94, 0.66, 0.98, 0.88, 0.37, 0.44, 0.67, 0.68, 0.5, 0.5, 0.88, 0.87],
                                         [0.92, 0.61, 1.0, 1.0, 0.99, 0.76, 1.0, 1.0, 0.39, 0.47, 0.78, 0.77, 0.43, 0.5, 0.99, 1.0],
                                         [0.92, 0.59, 1.0, 1.0, 0.96, 0.79, 1.0, 1.0, 0.38, 0.41, 0.78, 0.94, 0.45, 0.52, 1.0, 1.0],
                                         [0.92, 0.82, 1.0, 1.0, 0.98, 0.97, 1.0, 1.0, 0.55, 0.62, 1.0, 1.0, 0.63, 0.5, 1.0, 1.0]]).T),

     "R2": (("model", "test"), np.array([[1.0, 0.91, 1.0, 1.0, 1.0, 0.95,  1.0, 1.0, 0.8,  0.87, 1.0, 1.0, 0.79, 0.76, 1.0, 1.0],
                                         [1.0, 0.77, 1.0, 0.99, 1.0, 0.89, 1.0, 1.0, 0.76, 0.85, 0.89, 1.0, 0.81, 0.81, 1.0, 0.99],
                                         [1.0, 0.95, 1.0, 1.0, 1.0, 0.94,  1.0, 1.0, 0.82, 0.84, 1.0, 0.99, 0.86, 0.81, 1.0, 0.99],
                                         [1.0, 0.97, 1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 0.77, 0.96, 0.99, 1.0, 0.8, 0.73, 1.0, 1.0],
                                         [1.0, 1.0,  1.0, 1.0, 1.0, 0.98,  1.0, 1.0, 0.93, 0.93, 1.0, 1.0, 0.95, 0.95, 1.0, 1.0]]).T),

    "R3": (("model", "test"), np.array([[0.69, 0.66, 0.99, 0.99, 0.77, 0.76, 1.0,  1.0,  0.51, 0.7,  0.84, 0.98, 0.42, 0.48, 1.0,  0.98],
                                        [0.56, 0.68, 0.64, 0.72, 0.56, 0.65, 0.69, 0.81, 0.47, 0.6,  0.68, 0.68, 0.46, 0.4,  0.67, 0.77],
                                        [0.68, 0.7,  0.98, 1.0,  0.76, 0.74, 0.98, 1.0,  0.55, 0.68, 0.78, 0.75, 0.53, 0.48, 0.96, 0.95],
                                        [0.68, 0.8,  0.99, 1.0,  0.69, 0.68, 0.92, 1.0,  0.52, 0.59, 0.75, 0.99, 0.36, 0.49, 0.95, 0.97],
                                        [0.63, 0.53, 1.0, 0.98,  0.53, 0.6,  0.97, 0.92, 0.47, 0.53, 0.9,   1.0, 0.43, 0.37, 0.9,  0.9]]).T),

    "R4": (("model", "test"), np.array([[0.96, 0.95, 1.0,  1.0,  0.6,  0.99, 1.0,  1.0,  0.93, 0.99, 1.0,  1.0,  0.93, 0.9,  1.0,  1.0],
                                        [0.88, 0.95, 0.94, 0.98, 0.62, 0.99, 0.78, 0.97, 0.92, 1.0,  0.99, 1.0,  0.88, 0.88, 0.95, 0.92],
                                        [0.95, 0.97, 1.0,  1.0,  0.57, 0.94, 0.98, 0.94, 0.88, 0.96, 0.98, 0.94, 0.84, 0.83, 0.98, 0.98],
                                        [0.92, 0.99, 1.0,  1.0,  0.53, 1.0,  0.91, 1.0,  0.89, 1.0,  0.96, 1.0,  0.93, 0.95, 1.0,  1.0],
                                        [0.83, 0.77, 1.0,  1.0,  0.2,  0.25, 0.28, 0.47, 0.83, 0.95, 1.0,  0.87, 0.9,  0.9,  0.45, 0.38]]).T),

    "R5": (("model", "test"), np.array([[1.0, 0.97, 1.0,  1.0, 0.99, 0.99, 1.0, 1.0, 0.88, 0.99, 1.0,   1.0,  0.95, 0.99, 1.0, 1.0],
                                        [0.99, 1.0, 1.0, 0.99, 1.0,  0.98, 1.0, 0.99, 0.86, 0.94, 0.96, 0.99, 0.94, 0.96, 1.0, 1.0],
                                        [0.99, 1.0, 1.0,  1.0, 0.98, 1.0,  1.0, 1.0, 0.87, 0.97, 0.99,  0.99, 0.97, 0.96, 1.0, 1.0],
                                        [1.0, 0.99, 1.0,  1.0, 0.97, 1.0,  1.0, 1.0, 0.82, 0.94, 0.98,  1.0,  0.98, 0.95, 1.0, 1.0],
                                        [0.98, 1.0, 1.0,  1.0, 0.98, 1.0,  1.0, 1.0, 0.78, 0.83, 1.0,   1.0,  0.93,  0.9, 1.0, 1.0]]).T),

    "R6": (("model", "test"), np.array([[0.93, 0.74, 1.0,  1.0,  0.99, 0.9,  1.0,  1.0,  0.62, 0.61, 1.0,  1.0,  0.69, 0.69, 1.0,  1.0],
                                        [0.9,  0.62, 0.98, 0.95, 0.95, 0.69, 0.99, 0.89, 0.51, 0.64, 0.75, 0.93, 0.56, 0.56, 0.93, 0.95],
                                        [0.88, 0.76, 1.0,  1.0,  0.98, 0.81, 1.0,  1.0,  0.54, 0.6,  0.91, 0.9,  0.6,  0.53, 1.0,  1.0],
                                        [0.93, 0.81, 1.0,  1.0,  0.9,  0.85, 0.97, 1.0,  0.47, 0.68, 0.75, 1.0,  0.6,  0.59, 1.0,  0.99],
                                        [1.0,  0.97, 1.0,  1.0,  0.97, 0.98, 1.0,  1.0,  0.83, 0.92, 1.0,  1.0,  0.8,  0.82, 1.0,  0.98]]).T)},

    coords={
        "model": ['c00', 'c10', 'c01', 'c11', 'r00', 'r10', 'r01', 'r11', 'a00', 'a10', 'a01', 'a11', 'f00', 'f10', 'f01', 'f11'],
        "test": ['WT', 'ang_var', 'heel_jitter', 'rand_ablate', 'no_Lcell']}
)

df_3d = xarray_3d.to_dataframe()
print(df_3d)

# R3 data
# data = np.array([['', 'success', 'stiffness', 'nogozones', 'repuls', 'attract', 'const', 'full'],
#                  ['c',       69,           0,           0,        0,         0,       1],
#                  ['c_n',     66,           0,           1,        0,         0,       1],
#                  ['c_n_s',   99,           1,           1,        0,         0,       1],
#                  ['c_s',     99,           1,           0,        0,         0,       1],
#                  ['a_n',     74,           0,           1,        0,         1,       0],
#                  ['a',       51,           0,           0,        0,         1,       0],
#                  ['a_n_s',  100,           1,           1,        0,         1,       0],
#                  ['r_n',     76,           0,           1,        1,         0,       0],
#                  ['r',       77,           0,           0,        1,         0,       0],
#                  ['r_n_s',  100,           1,           1,        1,         0,       0],
#                  ['r_s',    100,           1,           0,        1,         0,       0]])

# pdf = pd.DataFrame(data=data[1:,1:],
#                    index=data[1:,0],LCELL_SIZE = np.array([[3.96, 2.23], [3.9, 2.28], [4.71, 2.61], [4.57, 2.63]])

#                    columns=data[0,1:])

pdf = df_3d['R3'].unstack()

X = pdf.assign(intercept=1).astype('float')
X = X.assign(fail = 100 - X['success'])
y = pd.concat([X.pop(x) for x in ['success', 'fail']], axis=1)

#%%
# try a model without any indicators
model_no_indicators = sm.GLM(y, X['intercept'], family=sm.families.Binomial())
result_no_indicators = model_no_indicators.fit()
print(result_no_indicators.summary())
print(result_no_indicators.aic)

yplot = y['success']/(y['success'] + y['fail'])

y_fitted_plot = np.array(result_no_indicators.fittedvalues, dtype='float')
plt.figure()
plt.plot(yplot, y_fitted_plot, 'o')
plt.plot(yplot, yplot, '--', label='y = x')

plt.ylabel("fitted value")
plt.xlabel("observed value")
plt.legend()

model_nogo = sm.GLM(y, X[['intercept', 'nogozones']],
                    family=sm.families.Binomial())
result_nogo = model_nogo.fit()
print(result_nogo.summary())
print(result_nogo.aic)

model_stiff = sm.GLM(y, X[['intercept', 'stiffness']],
                    family=sm.families.Binomial())
result_stiff = model_stiff.fit()
print(result_stiff.summary())
print(result_stiff.aic)
plt.figure()
plt.plot(yplot, result_stiff.fittedvalues, 'o')
plt.plot(yplot, yplot, '--', label='y = x')
plt.ylabel("fitted value")
plt.xlabel("observed value")
plt.legend()

model_attract = sm.GLM(y, X[['intercept', 'attract']],
                    family=sm.families.Binomial())
result_attract = model_attract.fit()
print(result_attract.summary())

print(result_attract.aic)

plt.figure()
plt.plot(yplot, result_attract.fittedvalues, 'o')
plt.plot(yplot, yplot, '--', label='y = x')
plt.ylabel("fitted value")
plt.xlabel("observed value")
plt.legend()

model_full = sm.GLM(y.astype('float'),
                    X[['intercept', 'stiffness', 'nogozones', 'repuls']].astype('float'),
                    family=sm.families.Binomial(link=sm.families.links.Logit()))
result_full = model_full.fit()
print(result_full.aic)

print(result_full.summary())


plt.figure()
plt.plot(yplot, result_full.fittedvalues, 'o')

plt.plot(yplot, yplot, '--', label='y = x')
plt.ylabel("fitted value")
plt.xlabel("observed value")
plt.legend()

#%%
# visualizations for heel jitter and random ablation
import numpy as np
import matplotlib.pyplot as plt

CONVFAC = 25.2
vec_opt = np.array([[8.83445973e+01, 7.67582595e+01, 1.70344600e+02, -6.24174303e+00],
                    [9.33712481e+01, 9.56606427e+01,
                        1.77259729e+02, 2.26227648e-02],

                    [1.00776263e+02, 9.85328622e+01,
                        1.98598157e+02, 1.44052931e+00],
                    [9.51118621e+01, 9.68803288e+01, 1.87098297e+02, 2.72226171e+00]])

species = 'dros'
v1, v2 = vec_opt[0, :2], vec_opt[0, 2:]
mini, maxi = -20, 20
n1, n2 = np.meshgrid(np.arange(mini, maxi), np.arange(mini, maxi))
Xint = (n1*v1[0] + n2*v2[0]).flatten()
Yint = (n1*v1[1] + n2*v2[1]).flatten()
Xrec, Yrec, Xc, Yc = create_starting_grid2(Xint, Yint, 20)
Xrec = np.ravel(Xrec)
Yrec = np.ravel(Yrec)

plt.figure()

plt.scatter(Xrec, Yrec, s=100, color=[[0.5, 0.5, 0.5]])
plt.scatter(Xrec + 0.3*CONVFAC*np.random.randn(*np.shape(Xrec)), Yrec + 0.3*CONVFAC*np.random.randn(*np.shape(Yrec)), s=100, c='k')
plt.axis([-200, 200, -200, 200])
plt.xticks([])
plt.yticks([])
plt.gca().set_aspect('equal')

l = len(Xrec)
Xrec_abl, Yrec_abl = Xrec.copy(), Yrec.copy()
perc_ablation = 25
rem = np.random.choice(l, int(l * perc_ablation/100), replace=False)
Xrec_abl[rem] = 1e4
Yrec_abl[rem] = 1e4

plt.figure()
plt.scatter(Xrec, Yrec, s=100, color=[[0.5, 0.5, 0.5]])
plt.scatter(Xrec_abl, Yrec_abl, s=100, c='k')
plt.axis([-200, 200, -200, 200])
plt.xticks([])
plt.yticks([])
plt.gca().set_aspect('equal')

#%%
# look at the growth of the lamina
def rotate_coord(x, y, angle):
    """
    angle in rad
    """
    return x*np.cos(angle) - y*np.sin(angle), x*np.sin(angle) + y*np.cos(angle)

LCELL_SIZE = np.array([[3.96, 2.23], [3.9, 2.28], [4.71, 2.61], [4.57, 2.63]])

CONVFAC = 25.2
vec_opt = np.array([[8.83445973e+01, 7.67582595e+01, 1.70344600e+02, -6.24174303e+00],
                    [9.33712481e+01, 9.56606427e+01,
                        1.77259729e+02, 2.26227648e-02],
                    [1.00776263e+02, 9.85328622e+01,
                        1.98598157e+02, 1.44052931e+00],
                    [9.51118621e+01, 9.68803288e+01, 1.87098297e+02, 2.72226171e+00]])
# rotate everything to have the grid perfectly horizontal
rot_ang = []
for hh in range(4):
    rot_ang.append(-np.arctan2(vec_opt[hh, 3], vec_opt[hh, 2]))
    # FRONTLOCS_ALL[:, hh, 0], FRONTLOCS_ALL[:, hh, 1] = rotate_coord(FRONTLOCS_ALL[:, hh, 0], FRONTLOCS_ALL[:, hh, 1], rot_ang[-1])
    # HEELLOCS_AVG[:, hh, 0], HEELLOCS_AVG[:, hh, 1] = rotate_coord(HEELLOCS_AVG[:, hh, 0], HEELLOCS_AVG[:, hh, 1], rot_ang[-1])
    vec_opt[hh, ::2], vec_opt[hh, 1::2] = rotate_coord(
        vec_opt[hh, ::2], vec_opt[hh, 1::2], rot_ang[-1])

v1, v2 = vec_opt[:, :2]/CONVFAC, vec_opt[:, 2:]/CONVFAC

# inter-bundle growth vs intra-bundle growth
devtime = [25, 30, 35, 40]
plt.figure()
plt.plot(devtime, (v2[:, 0] - LCELL_SIZE[:, 0]), label='inter-bundle distance (edge-edge)')
plt.plot(devtime, LCELL_SIZE[:, 0], label='L-cell size')
plt.plot(devtime, v2[:, 0], label='inter-bundle distance (center-center)')
plt.legend()
plt.ylim([0, 8])
plt.xlabel('Developmental time (hAPF)')
plt.ylabel('D-V Distance (um)')


plt.figure()
plt.plot(devtime, (v2[:, 0] - LCELL_SIZE[:, 0]), label='inter-bundle distance (edge-edge)')
plt.plot(devtime, LCELL_SIZE[:, 0], label='L-cell size')
plt.plot(devtime, v2[:, 0], label='inter-bundle distance (center-center)')
plt.legend()
plt.ylim([0, 8])
plt.xlabel('Developmental time (hAPF)')
plt.ylabel('Distance (um)')



plt.figure()
plt.plot(np.sqrt(v1[:, 0]**2 + v1[:, 1]**2)/CONVFAC, 'o--', color=[0.2, 0.2, 0.2])
plt.plot(np.sqrt(v2[:, 0]**2 + v2[:, 1]**2)/CONVFAC, 'o--', color=[0.2, 0.2, 0.2])
plt.plot(LCELL_SIZE[:, 0], 'o--', color=[0.8, 0.8, 0.8])
plt.plot(LCELL_SIZE[:, 1], 'o--', color=[0.9, 0.9, 0.9])
