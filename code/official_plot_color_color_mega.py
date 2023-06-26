# !/bin/python
# Plotting colour_colour plot

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import iqr
import scipy.stats as st

plt.style.use('style.mplstyle')

# Load in the table
hdulist = fits.open(
    '/net/vdesk/data2/bach1/ballieux/master_project_1/data/official_mega_master_clean.fits')
sample = 'VLASS'  # This van be either VLASS for the sample LoTSS NVSS VLASS, or LoLSS for LoLSS LoTSS NVSS
tbdata = hdulist[1].data
orig_cols = hdulist[1].columns
# print(orig_cols)
hdulist.close()

flux_sample = tbdata[sample + '_flux']


def hist_norm_height(n, bins, const):
    ''' Function to normalise bin height by a constant. 
        Needs n and bins from np.histogram.'''

    n = np.repeat(n, 2)
    n = np.float32(n) / const
    new_bins = [bins[0]]
    new_bins.extend(np.repeat(bins[1:], 2))

    return n, new_bins[:-1]


# Only use sources brighter than a certain limit
brightness_limit_list = [0., 0.05, 0.05, 0.1, 0.2, 0.5, 1., 2.]
bins_list = [50, 36, 30, 30, 25, 15, 5]
i_bins = 0

# range(len(bins_list)):, change range to plot for all brightness cuts or just the 0 mJy cut
for i_bins in range(1):
    brightness_limit = brightness_limit_list[i_bins]
    bins = bins_list[i_bins]
    # SNR = tbdata['VLASS_flux'] / (tbdata['VLASS_isl_rms'] / 1000)
    source_bright_ind = np.where((flux_sample > brightness_limit))

    alpha_low = tbdata['alpha_low_'+sample][source_bright_ind]
    alpha_high = tbdata['alpha_high_'+sample][source_bright_ind]


    err_alpha_low = np.median(tbdata['e_alpha_low_'+sample][source_bright_ind])
    err_alpha_high = np.median(tbdata['e_alpha_high_'+sample][source_bright_ind])
    print('err alpha_low:', err_alpha_low, 'err alpha_high:', err_alpha_high)

    # plt.clf()
    fig = plt.figure(0, figsize=(15, 15))
    gs = plt.GridSpec(2, 2, wspace=0, hspace=0, width_ratios=[
                      7, 1], height_ratios=[1, 7])
    ax2 = plt.subplot(gs[0])  # top histo
    ax1 = plt.subplot(gs[1:, -1])  # right histo
    ax = plt.subplot(gs[2])  # main

    # Limits of plot
    xmin = -2.5
    xmax = 3
    ymin = -2.5
    ymax = 1.5

    cmap = cm.get_cmap("gray")
    cmap._init()
    cmap._lut[:-3, :-1] = 0.
    cmap._lut[:-3, -1] = np.linspace(1, 0, cmap.N)

    X = np.linspace(xmin, xmax, bins + 1)
    Y = np.linspace(ymin, ymax, bins + 1)

    try:
        H, X, Y = np.histogram2d(
            alpha_low.flatten(), alpha_high.flatten(), bins=(X, Y))
        # weights=kwargs.get('weights', None))

    except ValueError:
        raise ValueError("It looks like at least one of your sample columns "
                         "have no dynamic range. You could try using the "
                         "`extent` argument.")

    V = 1.0 - np.exp(-0.475 * np.arange(0.5, 2.1, 0.5) ** 2)
    Hflat = H.flatten()
    inds = np.argsort(Hflat)[::-1]
    Hflat = Hflat[inds]
    sm = np.cumsum(Hflat)
    sm /= sm[-1]

    for i, v0 in enumerate(V):
        try:
            V[i] = Hflat[sm <= v0][-1]
        except:
            V[i] = Hflat[0]

    X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])
    X, Y = X[:-1], Y[:-1]


    # cs_scatter = ax.scatter(alpha_low_mps, alpha_high_mps,s=10, color='k', zorder=-1, c = duff_curve[high_limit_duff_ind],lw=0, cmap=plt.cm.copper,
    #     rasterized=True, norm=matplotlib.colors.LogNorm())
    # cs_scatter = ax.scatter(alpha_low_mps, alpha_high_mps,s=10, color='r',zorder=-1)

    # cs_scatter.cmap.set_under('k')
    # cs_scatter.set_clim(0.8,10.) # set limits of colurbar.
    # cax = fig.add_axes([0.15, 0.65, 0.2, 0.02]) # adding axes for colourbar to put in the plot
    # cbar = fig.colorbar(cs_scatter,cax = cax,orientation='horizontal',ticks=[0.8,0.9,1,2,3,4,5,6,7,8,9,10])
    # cbar.set_ticklabels([0.8,'',1,'','','',5,'','','','',10])

    if i_bins < 2:
        # ax.contourf(X1, Y1, H.T, [16,255], cmap=LinearSegmentedColormap.from_list("cmap",([1] * 3,[1] * 3),N=2), antialiased=False)
        if sample == 'VLASS':
            ax.contourf(X1, Y1, H.T, [104, 296, 916, 1975], cmap=LinearSegmentedColormap.from_list(
                    "cmap", ([1] * 3, [1] * 3), N=2), antialiased=False)
        elif sample == 'LoLSS':
            ax.contourf(X1, Y1, H.T, [24, 80, 234, 421], cmap=LinearSegmentedColormap.from_list(
                "cmap", ([1] * 3, [1] * 3), N=2), antialiased=False)
        #ax.contourf(X1, Y1, H.T, [15,75,145,265], cmap=LinearSegmentedColormap.from_list("cmap",([1] * 3,[1] * 3),N=2), antialiased=False)
        print('Use these values for contours', V[::-1])
        # [V[-1], H.max()]
        # if you don't want the grey scale beind, comment this line out.
        ax.pcolor(X, Y, H.max() - H.T, cmap=cmap, shading='auto', zorder=10)
        CS = ax.contour(X1, Y1, H.T, V[::-1], colors='k', linewidths=1.5)

    fake_cs_contour = ax.scatter(np.linspace(1000, H.max()), np.linspace(1000, H.max(
    )), c=np.linspace(0, H.max()), cmap=plt.cm.Greys)  # fake data to get the correct colourbar
    # adding axes for colourbar to put in the plot
    cax_contour = fig.add_axes([0.15, 0.68, 0.2, 0.02])
    if sample == 'VLASS':
        cbar_contour = fig.colorbar(fake_cs_contour, cax=cax_contour,
                                    orientation='horizontal', ticks=[0, 500, 1000, 1500, 2000])
        cbar_contour.ax.tick_params(size=4, labelsize=10)
        # because you need to subtract the value from the max above, the tick labels are inverted.
        cbar_contour.set_ticklabels([0, 500, 1000, 1500, 2000])
    else:
        cbar_contour = fig.colorbar(fake_cs_contour, cax=cax_contour,
                                    orientation='horizontal', ticks=[0, 100, 200, 300, 400])
        cbar_contour.ax.tick_params(size=4, labelsize=10)
        # because you need to subtract the value from the max above, the tick labels are inverted.
        cbar_contour.set_ticklabels([0, 100, 200, 300, 400])
    cbar_contour.solids.set_rasterized(True)

    ax.set_xlabel(r'$\alpha_{\mathrm{low}}$', fontsize=25)
    ax.set_ylabel(r'$\alpha_{\mathrm{high}}$', fontsize=25)
    ax.tick_params(left=True, labelleft=True, right=True,
                   bottom=True, labelbottom=True, top=True)

    # cs_scatter = ax.scatter(alpha_low_mps, alpha_high_mps,s=10, color='r')

    # Plotting lines
    compx = np.linspace(xmin + 0.1*xmin, xmax + 0.1*xmax)
    compy = np.zeros(len(compx))
    ax.plot(compx, compy, linestyle='-', color='k', linewidth=1., alpha=0.4)
    compy_2 = np.linspace(ymin + 0.1*ymin, ymax + 0.1*ymax)
    compx_2 = np.zeros(len(compy))
    ax.plot(compx_2, compy_2, linestyle='-',
            color='k', linewidth=1., alpha=0.4)

    y_diag = compx
    ax.plot(compx, y_diag, linestyle='--', color='red')

    if sample == 'VLASS':
        x_limit = err_alpha_low
        y_limit = - err_alpha_high
    else:
        x_limit = err_alpha_low
        y_limit = - err_alpha_high

    ax.vlines(x_limit, ymin, ymax, color='dodgerblue', zorder=10)
    ax.hlines(y_limit, xmin, xmax, color='dodgerblue', zorder=10)

    # GPS sel line
    #gps_x_sel = np.arange(-2.6,2.1,0.01)
    #gps_y_sel = np.ones(len(gps_x_sel))*-0.5
    #ax.plot(gps_x_sel,gps_y_sel,linestyle = '-',color = 'dodgerblue',linewidth = 2)

    # Tick labels for main plot
    ax.yaxis.set_ticks([-3, -2.5, -2, -1.5, -1., -0.5, 0, 0.5, 1, 1.5, 2])
    ax.yaxis.set_ticklabels([-3., '', -2., '', -1., '', 0., '', 1., '', 2.])
    ax.yaxis.set_ticks_position('both')

    # Plotting mock spectra

    # MPS Quarter
    x1_mock = np.linspace(1.2, 1.45)
    y1_mock = x1_mock - 3.55
    x2_mock = np.linspace(1.45, 1.7)
    y2_mock = -x1_mock - 0.9
    ax.plot(x1_mock, y1_mock, color='0.45')
    ax.plot(x2_mock, y2_mock, color='0.45')

    # Negative powerlaw
    x3_mock = np.linspace(-1.9, -1.65)
    y3_mock = y2_mock
    ax.plot(x3_mock, y3_mock, color='0.45')

    # Postive powerlaw
    x4_mock = np.linspace(1.35, 1.6)
    y4_mock = x4_mock + -0.25
    ax.plot(x4_mock, y4_mock, color='0.45')

    # Top left corner
    x5_mock = np.linspace(-1.75, -1.5)
    y5_mock = x5_mock + 2.85
    x6_mock = np.linspace(-2, -1.75)
    y6_mock = -x6_mock - 0.65
    ax.plot(x5_mock, y5_mock, color='0.45')
    ax.plot(x6_mock, y6_mock, color='0.45')

    (_, caps, _) = ax.errorbar(-2.15, .5, xerr=err_alpha_low,
                               yerr=err_alpha_high, color='k', linestyle='none', elinewidth=2)
    for cap in caps:
        cap.set_color('k')
        cap.set_markeredgewidth(2)

    if sample == "VLASS":
        ind_peaked = np.where((alpha_low > (tbdata['e_alpha_low_'+sample][source_bright_ind]))
                              & (alpha_high < - (tbdata['e_alpha_high_'+sample][source_bright_ind])))

    else:
        ind_peaked = np.where((alpha_low > (tbdata['e_alpha_low_'+sample][source_bright_ind]))
                              & (alpha_high < - (tbdata['e_alpha_high_'+sample][source_bright_ind])))
        # ind_peaked = np.where((alpha_low >= 0.) & (alpha_high <= 0.0))
    perc_peaked_source = (np.shape(ind_peaked)[
                          1] / np.shape(source_bright_ind)[1]) * 100
    ax.plot(alpha_low, alpha_high, "o", color='k',
            ms=3, zorder=-1, alpha=0.2, rasterized=True)
    #Do PS sources in different colors
    ax.plot(alpha_low[ind_peaked], alpha_high[ind_peaked], "o", color='rebeccapurple',
            ms=3, zorder=-1, alpha=0.2, rasterized=True)

    # ax.annotate(r'% peaked sources: ' + str(np.round(perc_peaked_source,2)), xy=(-2.4,1.70), xycoords='data', fontsize = 15)
    # ax.annotate(r'#N all sources: ' + str(np.round(np.shape(source_bright_ind)[1])), xy=(-2.4,1.50), xycoords='data', fontsize = 15)

    # Plot histograms on the side of the plot
    n_high, bins_high = np.histogram(alpha_high.flatten(), bins=Y)
    n_high, bins_high = hist_norm_height(n_high, bins_high, n_high.max())

    # Fit a gaussian to the histograms
    g_init = models.Gaussian1D(amplitude=1., mean=-0.82, stddev=1.)
    fit_g = fitting.LevMarLSQFitter()
    g_high = fit_g(g_init, bins_high, n_high)

    ax1.step(n_high, bins_high, color='k')
    ax1.plot(np.linspace(0, 1, 50), np.ones(50)*np.median(alpha_high), 'k--')
    ax1.plot(g_high(np.linspace(ymin, ymax, 100)),
             np.linspace(ymin, ymax, 100), 'r-')
    ax1.tick_params(left=True, right=True, bottom=True, top=True, labelsize=10)
    ax1.yaxis.set_ticks([-3, -2.5, -2, -1.5, -1., -0.5, 0, 0.5, 1, 1.5, 2])
    ax1.axes.set_yticklabels([])
    ax1.xaxis.set_ticks_position('top')
    ax1.set_ylim(ymin=ymin, ymax=ymax)
    ax1.set_xlim(xmin=0, xmax=1)
    ax1.xaxis.set_ticks([0., 0.25, 0.5, 0.75, 1])
    ax1.axes.set_xticklabels(['', 0.25, 0.5, 0.75, 1.0])

    n_low, bins_low = np.histogram(alpha_low.flatten(), bins=X)
    n_low, bins_low = hist_norm_height(n_low, bins_low, n_low.max())
    g_low = fit_g(g_init, bins_low, n_low)
    ax2.step(bins_low, n_low, color='k', zorder=0)
    ax2.plot(np.ones(50)*np.median(alpha_low), np.linspace(0, 1, 50), 'k--')
    ax2.plot(np.linspace(xmin, xmax, 100), g_low(
        np.linspace(xmin, xmax, 100)), 'r-', zorder=1)
    ax2.tick_params(left=True, right=True, bottom=True, top=True, labelsize=10)
    ax2.axes.set_xticklabels([])
    ax2.set_ylim(ymin=0, ymax=1)
    ax2.set_xlim(xmin, xmax)
    ax2.yaxis.set_ticks([0., 0.25, 0.5, 0.75, 1])
    ax2.axes.set_yticklabels(['', 0.25, 0.5, 0.75, 1.0])

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.grid(True)

    #ax2.set_title(r'Color-color diagram $\alpha_{low}, \alpha_{high}$. Brightness limit of ' + str(brightness_limit) + 'Jy', fontsize = 30)
    plt.savefig(
        '/net/vdesk/data2/bach1/ballieux/master_project_1/colour_colour_plots/official_colour_colour_'+sample+'_mega.png')
    # plt.show()

    #Something with interquartile range
    a_low_siqr = iqr(alpha_low)/2
    a_high_siqr = iqr(alpha_high)/2

    print('Brightness limit of ' + str(brightness_limit) +
          ' Jy has a total of ' + str(np.shape(source_bright_ind)[1]) + ' sources')
    print(str(np.round(perc_peaked_source, 2)) +
          ' % of sources are PS (', np.shape(ind_peaked)[1], ')')
    print('alpha_low median + std', np.median(alpha_low),
          np.std(alpha_low))
    print('alpha_high median + std', np.median(alpha_high),
          np.std(alpha_high))
    print('outliers<-2', len(alpha_high[alpha_high < -2]))
    print('Source with alpha_low > 2.5:', len(alpha_low[alpha_low > 2.5]), 'with coordinates', tbdata['RA'][source_bright_ind][alpha_low > 2.5],
      tbdata['DEC'][source_bright_ind][alpha_low>2.5], 'and alpha_low', alpha_low[alpha_low>2.5] , \
          '+/-' , tbdata['e_alpha_low_VLASS'][source_bright_ind][alpha_low>2.5], 'and alpha_high', alpha_high[alpha_low>2.5] , \
              '+/-' , tbdata['e_alpha_high_VLASS'][source_bright_ind][alpha_low>2.5] )
        
    #     #For LF
    # limit_95 = iqr(alpha_low, rng=(5,50)) 
    # print("95% distribution of alpha_low", np.median(alpha_low)-limit_95)
        #For HF
    limit_95 = iqr(alpha_low, rng=(50,95)) 
    print("95% distribution of alpha_low", np.median(alpha_low)+limit_95)