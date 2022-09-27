import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from collections import OrderedDict
import os
import math
from matplotlib import rc # To plot labels in serif rather than default.

#run on vdesk

cdict_models = {'powlaw': 'DarkOrange',
        'curve':'k'
        } 

plt.style.use('style.mplstyle')
rc('text', usetex=False)
rc('font', family='sans-serif')

# Defining plotting routine

def sed(models,paras,freq,flux,flux_err, name, name_title,#alpha_low, alpha_high,
        grid = False, freq_labels = False, log = True, bayes = False, resid = True, savefig=False, error = None):
    
    # Ensuring that freq and flux are approability matched up.
    ind = np.argsort(freq)
    freq = freq[ind]
    flux = flux[ind]
    flux_err = flux_err[ind]
    
    if resid == True:
        gs = plt.GridSpec(2,1, height_ratios = [3,1], hspace = 0)
        ax = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
        ax1.set_xlabel('Frequency (MHz)', fontsize = 20)
        ax1.set_ylabel(r'$\chi$', fontsize = 20)
        ax1.xaxis.labelpad = 15
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(2)
            ax1.spines[axis].set_linewidth(2)
    else:
        fig = plt.figure(1,figsize=(10, 8))#(12, 8))
        gs = plt.GridSpec(1,1)
        ax = plt.subplot(gs[0]) 
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    freq_cont = np.array(list(range(40,1700)))
    ax.set_xlim(40., 1500.)

    # Take highest flux in curved model as y-lim (if curved model is fitted):
    if math.isnan(paras[2][1]) == False:
        flux_max_curv = max(models[2](freq_cont, *paras[2]))
        flux_max = np.maximum(flux_max_curv, max(flux))
    else:
        flux_max = max(flux)

    ax.set_ylim(min(flux)-0.2*min(flux), flux_max + 0.2*flux_max)
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Flux Density (Jy)')
    ax.set_title(name_title, size = 22)
    try:
        tt = len(models)
    except TypeError:
        tt = 1
        models = [models]
        paras = [paras]

    for i in range(tt):

        # Defining colours for models to make it easy to tell the difference
        try:
            color = cdict_models[models[i].__name__] # In case none of the models are listed here.        
        except KeyError:
            color = 'darkgreen'

        ax.plot(freq_cont, models[i](freq_cont, *paras[i]), color = color)
        
        if resid == True:
            model_points = models[i](freq,*paras[i])
            residuals = flux-model_points
            chi_sing = residuals/flux_err
            chi_sing_err = np.ones(len(freq)) # Chi errors are 1.
            compx = np.linspace(min(freq)-0.1*min(freq),max(freq)+0.1*max(freq))
            compy = np.zeros(len(compx))
            ax1.plot(compx,compy,linestyle = '--',color = 'gray',linewidth = 2)
            ax1.set_xlim(min(freq_cont), max(freq_cont))
            # ax1.errorbar(freq,residuals,flux_err,color = color, linestyle='none',marker = '.')
            for i in range(len(freq)):
                if freq[i] in [54]:
                    (_, caps, _) = ax1.errorbar(freq[i], chi_sing[i], chi_sing_err[i], marker = 's',elinewidth=2, markersize=8, color = 'dodgerblue', linestyle='none', label = 'LoLSS')
                    for cap in caps:
                        cap.set_color('dodgerblue') 
                        cap.set_markeredgewidth(2)              
                elif freq[i] in [74]:
                    (_, caps, _) = ax1.errorbar(freq[i], chi_sing[i], chi_sing_err[i], marker = '^',elinewidth=2,  markersize=8, color = 'crimson', linestyle='none', label = 'VLSSr')
                    for cap in caps:
                        cap.set_color('crimson') 
                        cap.set_markeredgewidth(2)
                elif freq[i] in [150]:
                    (_, caps, _) = ax1.errorbar(freq[i], chi_sing[i], chi_sing_err[i], marker = 'o',elinewidth=2,  markersize=8, color = 'DarkGoldenRod', linestyle='none', label = 'TGSS')
                    for cap in caps:
                        cap.set_color('DarkGoldenRod') 
                        cap.set_markeredgewidth(2)        
                elif freq[i] in [144]:
                    (_, caps, _) = ax1.errorbar(freq[i], chi_sing[i], chi_sing_err[i], marker = '<',elinewidth=2,  markersize=8, color = 'forestgreen', linestyle='none', label = 'LoTSS')
                    for cap in caps:
                        cap.set_color('forestgreen') 
                        cap.set_markeredgewidth(2)
                elif freq[i] in [1400]:
                    (_, caps, _) = ax1.errorbar(freq[i], chi_sing[i], chi_sing_err[i], marker = 'v',elinewidth=2,  markersize=8, color = 'navy', linestyle='none', label = 'NVSS')
                    for cap in caps:
                        cap.set_color('navy') 
                        cap.set_markeredgewidth(2)

        if log == True:
            ax.set_xscale('log')
            ax.set_yscale('log')
            
            # Making sure minor ticks are marked properly.

            def ticks_format_x(value, index):
                """
                get the value and returns the value as:
                   integer: [0,99]
                   1 digit float: [0.1, 0.99]
                   n*10^m: otherwise
                To have all the number of the same size they are all returned as latex strings
                """
                exp = np.floor(np.log10(value))
                base = value/10**exp
                if value in np.array([50,70,90,400,600,800,900]): # This will remove 90 and 900 MHz, replace number for anything you don't want to appear.
                    return ''
                if exp == 0 or exp == 1 or exp == 2 or exp == 3:   
                    return '${0:d}$'.format(int(value))
                if exp == -1:
                    return '${0:.1f}$'.format(value)
                else:
                    return '${0:d}\\times10^{{{1:d}}}$'.format(int(base), int(exp))

            def ticks_format_y(value, index):
                """
                get the value and returns the value as:
                   integer: [0,99]
                   1 digit float: [0.1, 0.99]
                   n*10^m: otherwise
                To have all the number of the same size they are all returned as latex strings
                """
                if value in np.array([0.03,0.05,0.07,0.09,0.3,0.5,0.7,0.9,3,5,7,9,13,15,17,19,23,25,27,29]): # This will remove 90 and 900 MHz, replace number for anything you don't want to appear.
                    return ''
                exp = np.floor(np.log10(value))
                base = value/10**exp
                if exp == 0 or exp == 1 or exp == 2 or exp ==3:
                    return '${0:d}$'.format(int(value))
                if exp == -1:
                    return '${0:.1f}$'.format(value)
                if exp == -2:
                    return '${0:.2f}$'.format(value)
                else:
                    return '${0:d}\\times10^{{{1:d}}}$'.format(int(base), int(exp))

            subsx = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] # ticks to show per decade
            subsy = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] 
            ax.xaxis.set_minor_locator(ticker.LogLocator(subs=subsx)) #set the ticks position
            
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(ticks_format_x))   # remove the major ticks
            ax.xaxis.set_minor_formatter(ticker.FuncFormatter(ticks_format_x))

            ax.yaxis.set_minor_locator(ticker.LogLocator(subs=subsy))
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(ticks_format_y))   # remove the major ticks

            ax.yaxis.set_minor_formatter(ticker.FuncFormatter(ticks_format_y))


            ax.tick_params(axis='both',which='both', top=True,right=True)
            ax.tick_params(axis='both',which='major', width = 1.5)
            ax.tick_params(axis='both',which='minor', length = 5,width = 1.5)

            if resid == True:
                ax.set_xticklabels('',minor = True)
                ax.set_xlabel('')
                ax1.set_xscale('log')
                ax1.xaxis.set_minor_locator(ticker.LogLocator(subs=subsx)) #set the ticks position
                ax1.xaxis.set_major_formatter(ticker.NullFormatter())   # remove the major ticks
                ax1.xaxis.set_minor_formatter(ticker.FuncFormatter(ticks_format_x))
                ax1.tick_params(axis='both',which='both')
                ax1.tick_params(axis='both',which='major')
                ax1.tick_params(axis='both',which='minor')

    if freq_labels == True:
        
        for i in range(len(freq)):
            LoLSS_range = np.array([[54-42, 66-54]]).T# Range of observational frequencies LoLSS
            LoTSS_range = np.array([[144-120, 168-144]]).T
            FIRST_range = np.array([[1400-1365,1435-1400]]).T
            
            if freq[i] in [54]:
                (_, caps, _) = ax.errorbar(freq[i], flux[i], xerr = LoLSS_range, yerr = flux_err[i], marker = 's',elinewidth=2, markersize=8,  color = 'dodgerblue', linestyle='none', label = 'LoLSS',markeredgecolor='none')
                for cap in caps:
                    cap.set_color('dodgerblue') 
                    cap.set_markeredgewidth(2)
            elif freq[i] in [74]:
                (_, caps, _) = ax.errorbar(freq[i], flux[i], flux_err[i],marker = 'o',elinewidth=2, markersize=10,  color = 'crimson', linestyle='none', label = 'VLSSr',markeredgecolor='none')
                for cap in caps:
                    cap.set_color('crimson') 
                    cap.set_markeredgewidth(2)                
            elif freq[i] in [150]:
                (_, caps, _) = ax.errorbar(freq[i], flux[i], flux_err[i], marker = '>',elinewidth=2,  markersize=10, color = 'DarkGoldenRod', linestyle='none', label = 'TGSS',markeredgecolor='none')
                for cap in caps:
                    cap.set_color('DarkGoldenRod') 
                    cap.set_markeredgewidth(2)    
            elif freq[i] in [144]:
                (_, caps, _) = ax.errorbar(freq[i], flux[i], xerr = LoTSS_range, yerr = flux_err[i], marker = '<',elinewidth=2, markersize=10,  color = 'forestgreen', linestyle='none', label = 'LoTSS',markeredgecolor='none')
                for cap in caps:
                    cap.set_color('forestgreen') 
                    cap.set_markeredgewidth(2)
            elif freq[i] in [1400]:
                (_, caps, _) = ax.errorbar(freq[i], flux[i], flux_err[i], marker = 'v',elinewidth=2, markersize=10,  color = 'navy', linestyle='none', label = 'NVSS',markeredgecolor='none')
                for cap in caps:
                    cap.set_color('navy') 
                    cap.set_markeredgewidth(2)
            elif freq[i] in [1400.001]:
                (_, caps, _) = ax.errorbar(freq[i], flux[i], xerr = FIRST_range, yerr = flux_err[i], marker = '^',elinewidth=2, markersize=10,  color = 'hotpink', linestyle='none', label = 'FIRST',markeredgecolor='none')
                for cap in caps:
                    cap.set_color('hotpink') 
                    cap.set_markeredgewidth(2)
            elif freq[i] in [128., 144.00001, 160.]:
                (_, caps, _) = ax.errorbar(freq[i], flux[i], xerr = np.array([[8, 8]]).T, yerr = flux_err[i], marker = 'D',elinewidth=2, markersize=10,  color = 'indigo', linestyle='none', label = 'LoTSS in-band',markeredgecolor='none')
                for cap in caps:
                    cap.set_color('indigo') 
                    cap.set_markeredgewidth(2)
            else:
                ax.errorbar(freq[i], flux[i], flux_err[i], marker = 'o',elinewidth=2, markersize=10, color = 'darkgreen', linestyle='none', label = 'Data',markeredgecolor='none')
                    # Elimanating doubled up legend values.
            # handles, labels = ax.get_legend_handles_labels()
            # by_label = OrderedDict(zip(labels, handles))
            # ax.legend(by_label.values(), by_label.keys(),loc='upper right')

    else:   
        ax.errorbar(freq, flux, flux_err, marker = '.', color = 'darkgreen', linestyle='none', label = 'Data')
    
    if grid == True:
        ax.grid(which='both')
        if resid == True:
            ax1.grid(axis = 'x',which = 'both')

    if savefig == True:
        if not os.path.exists(os.getcwd()+'/kross_seds'):
            os.makedirs(os.getcwd()+'/kross_seds')
            print('Creating directory ', os.getcwd()+'/kross_seds/ and saving figures in png format with title names.')
        plt.savefig('kross_seds/'+name+'.png')
    # plt.show()
