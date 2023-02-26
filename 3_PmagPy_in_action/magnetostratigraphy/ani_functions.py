import pmagpy.pmag as pmag
import pmagpy.pmagplotlib as pmagplotlib
from pmagpy import contribution_builder as cb
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def ani_depthplot(spec_file='specimens.txt', samp_file='samples.txt',
                  meas_file='measurements.txt', site_file='sites.txt',
                  age_file="", sum_file="", fmt='svg', dmin=-1, dmax=-1,
                  depth_scale='core_depth', dir_path='.', contribution=None):
    """
    returns matplotlib figure with anisotropy data plotted against depth
    available depth scales: 'composite_depth', 'core_depth' or 'age' (you must provide an age file to use this option).
    You must provide valid specimens and sites files, and either a samples or an ages file.
    You may additionally provide measurements and a summary file (csv).

    Parameters:
        spec_file (str): default "specimens.txt"
        samp_file (str): default "samples.txt"
        meas_file (str): default "measurements.txt"
        site_file (str): default "sites.txt"
        age_file (str): default ""
        sum_file (str): default ""
        fmt (str): str, default "svg"
            format for figures, ["svg", "jpg", "pdf", "png"]
        dmin (number): default -1
            minimum depth to plot (if -1, default to plotting all)
        dmax (number): default -1
            maximum depth to plot (if -1, default to plotting all)
        depth_scale (str): default "core_depth"
            scale to plot, ['composite_depth', 'core_depth', 'age'].
            if 'age' is selected, you must provide an ages file.
        dir_path (str): default "."
            directory for input files
        contribution : cb.Contribution, default None
            if provided, use Contribution object instead of reading in
            data from files

    Returns:
        plot 
            matplotlib plot, or False if no plot could be created
        name 
            figure name, or error message if no plot could be created
    """
    if depth_scale == 'sample_core_depth':
        depth_scale = 'core_depth'
    if depth_scale == 'sample_composite_depth':
        depth_scale = 'composite_depth'

    pcol = 4
    tint = 9
    plots = 0

    dmin, dmax = float(dmin), float(dmax)

    # if contribution object is not provided, read in data from files
    if isinstance(contribution, cb.Contribution):
        con = contribution
    else:
        # format files to use full path
        meas_file = pmag.resolve_file_name(meas_file, dir_path)
        spec_file = pmag.resolve_file_name(spec_file, dir_path)
        samp_file = pmag.resolve_file_name(samp_file, dir_path)
        site_file = pmag.resolve_file_name(site_file, dir_path)

        if age_file:
            age_file = pmag.resolve_file_name(age_file, dir_path)
            if not os.path.isfile(age_file):
                print(
                    'Warning: you have provided an invalid age file.  Attempting to use sample file instead')
                age_file = None
                depth_scale = 'core_depth'
            else:
                samp_file = age_file
                depth_scale = 'age'
                print(
                    'Warning: you have provided an ages format file, which will take precedence over samples')

        samp_file = pmag.resolve_file_name(samp_file, dir_path)

        label = 1

        if sum_file:
            sum_file = pmag.resolve_file_name(sum_file, dir_path)
            core_df=pd.read_csv(sum_file)
            depths=core_df['Top depth cored CSF (m)'].values

        # contribution

        dir_path = os.path.split(spec_file)[0]
        tables = ['measurements', 'specimens', 'samples', 'sites']
        con = cb.Contribution(dir_path, read_tables=tables,
                              custom_filenames={'measurements': meas_file, 'specimens': spec_file,
                                                'samples': samp_file, 'sites': site_file})



    for ftype in ['specimens', 'samples', 'sites']:
        if not con.tables.get(ftype):
            if ftype == 'samples':
                if con.tables.get('ages'):
                    depth_scale = 'age'
                    continue
            print("-W- This function requires a {} file to run.".format(ftype))
            print("    Make sure you include one in your working directory")
            return False, "missing required file type: {}".format(ftype)


    # propagate needed values
    con.propagate_cols(['core_depth'], 'samples', 'sites')
    con.propagate_location_to_specimens()

    # get data read in
    isbulk = 0  # tests if there are bulk susceptibility measurements
    ani_file = spec_file

    SampData = con.tables['samples'].df
    AniData = con.tables['specimens'].df
    # add sample into specimens (AniData)
    AniData = pd.merge(
        AniData, SampData[['sample', depth_scale]], how='inner', on='sample')
    # trim down AniData
    cond = AniData[depth_scale].astype(bool)
    AniData = AniData[cond]
    if dmin != -1:
        AniData = AniData[AniData[depth_scale] < dmax]
    if dmax != -1:
        AniData = AniData[AniData[depth_scale] > dmin]
    AniData['core_depth'] = AniData[depth_scale]

    if not age_file:
        Samps = con.tables['samples'].convert_to_pmag_data_list()
    else:
        con.add_magic_table(dtype='ages', fname=age_file)
        Samps = con.tables['ages'].convert_to_pmag_data_list()
        # get age unit
        age_unit = con.tables['ages'].df['age_unit'][0]
        # propagate ages down to sample level

    for s in Samps:
        # change to upper case for every sample name
        s['sample'] = s['sample'].upper()

    if 'measurements' in con.tables:
        isbulk = 1
        Meas = con.tables['measurements'].df  # convert_to_pmag_data_list()

    if isbulk:
        Meas = Meas[Meas['specimen'].astype('bool')]
        Meas = Meas[Meas['susc_chi_volume'].astype(bool)]
        # add core_depth into Measurements dataframe
        Meas = pd.merge(Meas[['susc_chi_volume', 'specimen']], AniData[[
                        'specimen', 'core_depth']], how='inner', on='specimen')
        Bulks = list(Meas['susc_chi_volume'] * 1e6)
        BulkDepths = list(Meas['core_depth'])
    else:
        Bulks, BulkDepths = [], []

    # now turn Data from pandas dataframe to a list of dicts
    Data = list(AniData.T.apply(dict))

    if len(Bulks) > 0:  # set min and max bulk values
        bmin = min(Bulks)
        bmax = max(Bulks)
    xlab = "Depth (m)"

    #
    if len(Data) > 0:
        location = Data[0].get('location', 'unknown')
        if cb.is_null(location):
            location = 'unknown'
            try:
                location = con.tables['sites'].df['location'][0]
            except KeyError:
                pass
    else:
        return False, 'no data to plot'

    # collect the data for plotting tau  V3_inc and V1_dec
    Depths, Tau1, Tau2, Tau3, V3Incs, P, V1Decs = [], [], [], [], [], [], []
    F23s = []
    Axs = []  # collect the plot ids
    if len(Bulks) > 0:
        pcol += 1

    Data = pmag.get_dictitem(Data, 'aniso_s', '', 'not_null')
    # get all the s1 values from Data as floats
    aniso_s = pmag.get_dictkey(Data, 'aniso_s', '')
    aniso_s = [a.split(':') for a in aniso_s if a is not None]
    #print('aniso_s', aniso_s)
    s1 = [float(a[0]) for a in aniso_s]
    s2 = [float(a[1]) for a in aniso_s]
    s3 = [float(a[2]) for a in aniso_s]
    s4 = [float(a[3]) for a in aniso_s]
    s5 = [float(a[4]) for a in aniso_s]
    s6 = [float(a[5]) for a in aniso_s]
    # we are good with s1 - s2
    nmeas = pmag.get_dictkey(Data, 'aniso_s_n_measurements', 'int')
    sigma = pmag.get_dictkey(Data, 'aniso_s_sigma', 'f')
    Depths = pmag.get_dictkey(Data, 'core_depth', 'f')
    # Ss=np.array([s1,s4,s5,s4,s2,s6,s5,s6,s3]).transpose() # make an array
    Ss = np.array([s1, s2, s3, s4, s5, s6]).transpose()  # make an array
    # Ts=np.reshape(Ss,(len(Ss),3,-1)) # and re-shape to be n-length array of
    # 3x3 sub-arrays
    for k in range(len(Depths)):
        # tau,Evecs= pmag.tauV(Ts[k]) # get the sorted eigenvalues and eigenvectors
        # v3=pmag.cart2dir(Evecs[2])[1] # convert to inclination of the minimum
        # eigenvector
        fpars = pmag.dohext(nmeas[k] - 6, sigma[k], Ss[k])
        V3Incs.append(fpars['v3_inc'])
        V1Decs.append(fpars['v1_dec'])
        Tau1.append(fpars['t1'])
        Tau2.append(fpars['t2'])
        Tau3.append(fpars['t3'])
        P.append(Tau1[-1]/Tau3[-1])
        F23s.append(fpars['F23'])
    if len(Depths) > 0:
        if dmax == -1:
            dmax = max(Depths)
            dmin = min(Depths)
        tau_min = 1
        for t in Tau3:
            if t > 0 and t < tau_min:
                tau_min = t
        tau_max = max(Tau1)
        # tau_min=min(Tau3)
        P_max = max(P)
        P_min = min(P)
        # dmax=dmax+.05*dmax
        # dmin=dmin-.05*dmax

        main_plot = plt.figure(1, figsize=(11, 7))  # make the figure
        # main_plot = plt.figure(1, figsize=(10, 8))  # make the figure

        version_num = pmag.get_version()
        plt.figtext(.02, .01, version_num)  # attach the pmagpy version number
        ax = plt.subplot(1, pcol, 1)  # make the first column
        Axs.append(ax)
        ax.plot(Tau1, Depths, 'rs')
        ax.plot(Tau2, Depths, 'b^')
        ax.plot(Tau3, Depths, 'ko')
        if sum_file:
            for depth in depths:
                if depth >= dmin and depth < dmax:
                    plt.axhline(depth,color='blue',linestyle='dotted')
        if tau_min>.3: tau_min=.3
        if tau_max<.36: tau_max=.36
        ax.axis([tau_min, tau_max, dmax, dmin])
        ax.set_xlabel('Eigenvalues')
        if depth_scale == 'core_depth':
            ax.set_ylabel('Depth (mbsf)')
        elif depth_scale == 'age':
            ax.set_ylabel('Age (' + age_unit + ')')
        else:
            ax.set_ylabel('Depth (mcd)')
        ax2 = plt.subplot(1, pcol, 2)  # make the second column
        ax2.yaxis.set_major_locator(plt.NullLocator())
        ax2.plot(P, Depths, 'rs')
        ax2.axis([P_min, P_max, dmax, dmin])
        ax2.set_xlabel('P')
        ax2.set_title(location)
        if sum_file:
            for depth in depths:
                if depth >= dmin and depth < dmax:
                    plt.axhline(depth,color='blue',linestyle='dotted')
        Axs.append(ax2)
        ax3 = plt.subplot(1, pcol, 3)
        Axs.append(ax3)
        ax3.plot(V3Incs, Depths, 'ko')
        ax3.axis([0, 90, dmax, dmin])
        ax3.set_xlabel('V3 Inclination')
        ax3.yaxis.set_major_locator(plt.NullLocator())
        if sum_file:
            for depth in depths:
                if depth >= dmin and depth < dmax:
                    plt.axhline(depth,color='blue',linestyle='dotted')
        ax4 = plt.subplot(1, np.abs(pcol), 4)
        Axs.append(ax4)
        ax4.plot(V1Decs, Depths, 'rs')
        ax4.axis([0, 360, dmax, dmin])
        ax4.set_xlabel('V1 Declination')
        ax4.yaxis.set_major_locator(plt.NullLocator())
        if sum_file:
            for depth in depths:
                if depth >= dmin and depth < dmax:
                    plt.axhline(depth,color='blue',linestyle='dotted')
        # ax5=plt.subplot(1,np.abs(pcol),5)
        # Axs.append(ax5)
        # ax5.plot(F23s,Depths,'rs')
        # bounds=ax5.axis()
        # ax5.axis([bounds[0],bounds[1],dmax,dmin])
        # ax5.set_xlabel('F_23')
        # ax5.semilogx()
        # if sum_file:
        #    for core in Cores:
        #         depth=float(core[core_depth_key])
        #         if depth>=dmin and depth<=dmax:
        #            plt.plot([bounds[0],bounds[1]],[depth,depth],'b--')
        #            if pcol==5 and label==1:plt.text(bounds[1],depth+tint,core[core_label_key])
        # if pcol==6:
        if pcol == 5:
            # ax6=plt.subplot(1,pcol,6)
            ax6 = plt.subplot(1, pcol, 5)
            Axs.append(ax6)
            ax6.plot(Bulks, BulkDepths, 'bo')
            #ax6.axis([bmin - 1, 1.1 * bmax, dmax, dmin])
            ax6.set_ylim(dmax, dmin)
            ax6.set_xlabel('Bulk Susc. (uSI)')
            ax6.yaxis.set_major_locator(plt.NullLocator())
            if sum_file:
                for depth in depths:
                    if depth >= dmin and depth < dmax:
                        plt.axhline(depth,color='blue',linestyle='dotted')
        for x in Axs:
            # this makes the x-tick labels more reasonable - they were
            # overcrowded using the defaults
            pmagplotlib.delticks(x)
        fig_name = location + '_ani_depthplot.' + fmt
        return main_plot, [fig_name]
    else:
        return False, "No data to plot"