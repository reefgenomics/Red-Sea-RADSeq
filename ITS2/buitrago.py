#!/usr/bin/env python3

from sputils.spbars import SPBars
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import os

class Buitrago:
    def __init__(self):
        self.root_dir = os.path.dirname(os.path.abspath( __file__ ))
        self.pver_df = pd.read_csv(
            filepath_or_buffer=os.path.join(self.root_dir, 'pver.14reef.277ind.ordered.strata.txt'),
            index_col=0, sep='\t'
        )
        self.pver_df = self.pver_df.iloc[::-1]
        self.spis_df = pd.read_csv(
            filepath_or_buffer=os.path.join(self.root_dir, 'spis.18reef.368ind.ordered.strata.txt'),
            index_col=0, sep='\t'
        )
        self.spis_df = self.spis_df.iloc[::-1]
        self.spis_df.drop(labels=['SWAJ-R1-43', 'SMAQ-R1-30'], axis=0, inplace=True)

        self.fig = plt.figure(figsize=self._mm2inch((200, 320)))
        # bars to legends at ratio of 4:1
        gs = gridspec.GridSpec(nrows=4, ncols=6)
        # TODO we are here, setting up the axes that we will then pass into the sputils.
        # TODO we still need to write the legend making for the utils module and use it here.
        # Then we need to look at ordinations and hierarchical clustering and annotating accorrding to the site
        # and reef.
        # There will be 6 axes for plotting and 3 legend axes
        # 2 sets of 3 one for each
        self.titles = ['pver_genera', 'pver_seq', 'pver_profile', 'spis_genera', 'spis_seq', 'spis_profile']
        self.bar_ax_arr = [
            # pver_genera_ax
            plt.subplot(gs[:4, :1]),
            # pver_seq_ax
            plt.subplot(gs[:4, 1:2]),
            # pver_profile_ax
            plt.subplot(gs[:4, 2:3]),
            # spis_genera_ax
            plt.subplot(gs[:4, 3:4]),
            # spis_seq_ax
            plt.subplot(gs[:4, 4:5]),
            # spis_profile_ax
            plt.subplot(gs[:4, 5:6]),
        ]

        # self.genera_leg_ax = plt.subplot(gs[4:5, :2])
        # self.seq_leg_ax = plt.subplot(gs[4:5, 2:4])
        # self.profile_leg_ax = plt.subplot(gs[4:5, 4:6])

        # Count table relative paths
        self.seq_count_table_path = os.path.join(self.root_dir, 'sp_output/post_med_seqs/131_20201203_DBV_20201207T095144.seqs.absolute.abund_and_meta.txt')
        self.profile_count_table_path = os.path.join(self.root_dir, 'sp_output/its2_type_profiles/131_20201203_DBV_20201207T095144.profiles.absolute.abund_and_meta.txt')

        # create an instance of SPBars just to generate a seq and profile dict for the whole dataset
        # then use this dictionary for plotting the actual plots.
        spb = SPBars(
            seq_count_table_path=self.seq_count_table_path,
            profile_count_table_path=self.profile_count_table_path,
            plot_type='seq_and_profile', orientation='v', legend=False, relative_abundnce=True
        )
        self.seq_color_dict = spb.seq_color_dict
        self.profile_color_dict = spb.profile_color_dict

        config_tups = [
            ('seq_only', self.seq_color_dict, None, True),
            ('seq_only', self.seq_color_dict, None, False),
            ('profile_only', None, self.profile_color_dict, False)]
        # Now we can plot up each of the axes
        for i, species_sample_lists in enumerate([self.pver_df.index.values, self.spis_df.index.values]):
            for j, (plot_type, seq_color_dict, profile_color_dict, color_by_genus) in enumerate(config_tups):
                sp_bars = SPBars(
                    seq_count_table_path=self.seq_count_table_path,
                    profile_count_table_path=self.profile_count_table_path,
                    plot_type=plot_type, orientation='v', legend=False, relative_abundnce=True,
                    color_by_genus=color_by_genus, sample_outline=False, sample_names_included=species_sample_lists,
                    bar_ax=self.bar_ax_arr[(i*3) + j], seq_color_dict=seq_color_dict,
                    profile_color_dict=profile_color_dict
                )
                sp_bars.plot()

        for ax, title in zip(self.bar_ax_arr, self.titles):
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(title, fontsize='small')
        self.fig.show()
        plt.savefig('bar_plots.svg')
        plt.savefig('bar_plots.png', dpi=1200)
        # TODO we are here

        # put horizontal lines to show the reef boundaries
        # TODO but should probably prioritise getting some ordinations plotted up

        foo = 'bar'

    def _mm2inch(self, *tupl):
        inch = 25.4
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

Buitrago()


