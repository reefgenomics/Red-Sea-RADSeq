#!/usr/bin/env python3

from sputils.spbars import SPBars
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class Buitrago:
    def __init__(self):
        with open('/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/pver.14reef.277ind.ordered.strata.txt',
                  'r') as f:
            self.pver_sample_list = [line.split('\t')[0] for line in list(f)[1:]]

        with open('/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/spis.18reef.368ind.ordered.strata.txt',
                  'r') as f:
            self.spis_sample_list = [line.split('\t')[0] for line in list(f)[1:]]
        self.spis_sample_list.remove('SWAJ-R1-43')
        self.spis_sample_list.remove('SMAQ-R1-30')

        foo = 'bar'
        # not found spis {'SWAJ-R1-43', 'SMAQ-R1-30'}
        self.fig = plt.figure(figsize=self._mm2inch((250, 320)))
        # bars to legends at ratio of 4:1
        gs = gridspec.GridSpec(nrows=5, ncols=6)
        # TODO we are here, setting up the axes that we will then pass into the sputils.
        # TODO we still need to write the legend making for the utils module and use it here.
        # Then we need to look at ordinations and hierarchical clustering and annotating accorrding to the site
        # and reef.
        # There will be 6 axes for plotting and 3 legend axes
        # 2 sets of 3 one for each
        self.pver_genera_ax = plt.subplot(gs[:4, :1])
        self.pver_seq_ax = plt.subplot(gs[:4, 1:2])
        self.pver_profile_ax = plt.subplot(gs[:4, 2:3])
        self.spis_genera_ax = plt.subplot(gs[:4, 3:4])
        self.spis_seq_ax = plt.subplot(gs[:4, 4:5])
        self.spis_profile_ax = plt.subplot(gs[:4, 5:6])
        self.genera_leg_ax = plt.subplot(gs[4:5, :2])
        self.seq_leg_ax = plt.subplot(gs[4:5, 2:4])
        self.profile_leg_ax = plt.subplot(gs[4:5, 4:6])


        spb = SPBars(seq_count_table_path='/Users/benjaminhume/Documents/projects/20210113_buitrago/sp_output/post_med_seqs/131_20201203_DBV_20201207T095144.seqs.absolute.abund_and_meta.txt',
            profile_count_table_path='/Users/benjaminhume/Documents/projects/20210113_buitrago/sp_output/its2_type_profiles/131_20201203_DBV_20201207T095144.profiles.absolute.abund_and_meta.txt',
            plot_type='seq_and_profile', orientation='v', legend=False, relative_abundnce=True)
        self.seq_color_dict = spb.seq_color_dict
        self.profile_color_dict = spb.profile_color_dict
        # create an instance of SPBars just to generate a seq and profile dict for the whole dataset
        # then use this dictionary for plotting the actual plots.

        sp_bars = SPBars(
            seq_count_table_path='/Users/benjaminhume/Documents/projects/20210113_buitrago/sp_output/post_med_seqs/131_20201203_DBV_20201207T095144.seqs.absolute.abund_and_meta.txt',
            profile_count_table_path='/Users/benjaminhume/Documents/projects/20210113_buitrago/sp_output/its2_type_profiles/131_20201203_DBV_20201207T095144.profiles.absolute.abund_and_meta.txt',
            plot_type='seq_only', orientation='v', legend=False, relative_abundnce=True,
            color_by_genus=True, sample_outline=True, sample_names_included=self.pver_sample_list, bar_ax=self.ax_arr[0]
        )

        sp_bars.plot()

        foo = 'bar'

    def _mm2inch(self, *tupl):
        inch = 25.4
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

Buitrago()


