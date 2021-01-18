#!/usr/bin/env python3

from sputils.spbars import SPBars
from sputils.sphierarchical import SPHierarchical
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import os
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import numpy as np

class Buitrago:
    """A base class that will give access to the basic meta info dfs"""
    def __init__(self):
        self.root_dir = os.path.dirname(os.path.abspath(__file__))
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
        self.sample_names = list(self.spis_df.index) + list(self.pver_df.index)

        # Count table relative paths
        self.seq_count_table_path = os.path.join(self.root_dir,
                                                 'sp_output/post_med_seqs/131_20201203_DBV_20201207T095144.seqs.absolute.abund_and_meta.txt')
        self.profile_count_table_path = os.path.join(self.root_dir,
                                                     'sp_output/its2_type_profiles/131_20201203_DBV_20201207T095144.profiles.absolute.abund_and_meta.txt')

    def _mm2inch(self, *tupl):
        inch = 25.4
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

class BuitragoHier(Buitrago):
    """
    Plot up a series of dendrograms
    """
    def __init__(self):
        super().__init__()

        # setup fig
        # 6 rows for the dendro and 1 for the coloring by species
        gs = gridspec.GridSpec(nrows=8, ncols=1)

        self.fig = plt.figure(figsize=(self._mm2inch(183, 80)))
        self.dendro_ax = plt.subplot(gs[:4, :])
        self.species_ax = plt.subplot(gs[4:5, :])
        self.seq_bars_ax = plt.subplot(gs[5:8, :])

        dist_path = 'sp_output/between_sample_distances/A/20201207T095144_braycurtis_sample_distances_A_sqrt.dist'

        # We have to work out which samples to plot.
        # I think that some samples were not used in the end for the study due to the
        # host data being no good.
        # Also, there will be samples not included in the dist matrix as they contained no A.

        # Run SPHier through blank to get the list of samples we have in the A matrix
        # THen find the interset of samples listed in the self.pver and self.spis dfs.
        sph = SPHierarchical(dist_output_path=dist_path, no_plotting=True)
        symbiodinium_names = sph.obj_name_to_obj_uid_dict.keys()
        symbiodinium_host_names = set(symbiodinium_names).intersection(set(self.sample_names))

        sph = SPHierarchical(
            dist_output_path=dist_path, ax=self.dendro_ax,
            sample_names_included=symbiodinium_host_names)
        sph.plot()
        self.dendro_ax.collections[0].set_linewidth(0.5)
        self.sample_uid_to_sample_name = sph.obj_uid_to_obj_name_dict
        self.sample_name_to_sample_uid = sph.obj_name_to_obj_uid_dict

        self._plot_species_ax(sph)

        # We want to plot the bars in the order of the hierarchical
        # we will use the sph.dendrogram['ivl'] uids converted to names
        dendrogram_sample_name_order = [self.sample_uid_to_sample_name[_] for _ in sph.dendrogram['ivl']]

        # Now plot the bars
        spb = SPBars(
            seq_count_table_path=self.seq_count_table_path,
            profile_count_table_path=self.profile_count_table_path,
            plot_type='seq_only', orientation='h', legend=False, relative_abundance=True,
            sample_names_included=dendrogram_sample_name_order, bar_ax=self.seq_bars_ax
        )
        spb.plot()
        self.seq_bars_ax.set_xticks([])
        self.seq_bars_ax.set_yticks([])
        self.seq_bars_ax.set_ylabel('sequences', rotation='vertical', fontsize='x-small')


        plt.savefig('dendro_bars.svg')
        plt.savefig('dendro_bars.png', dpi=1200)
        foo = 'bar'

    def _plot_species_ax(self, sph):
        # plot the species colors on the bottom axis
        # Rectangles should be centered around the x coordinate of the leaves. Width will be
        # the distance between the x coordinates.
        x_coords = range(5, (len(sph.dendrogram['ivl']) * 10) + 5, 10)
        sample_name_to_x_coord_dict = {
            sample_name: x_coord for
            sample_name, x_coord in
            zip(sph.dendrogram['ivl'], x_coords)
        }
        c_dict = {'P': 'brown', 'S': 'blue'}
        width = 10
        color_list = []
        rectangles = []
        for sample_uid, x_coord in sample_name_to_x_coord_dict.items():
            if self.sample_uid_to_sample_name[sample_uid][0] in ['S', 'P']:
                rectangles.append(Rectangle(
                    (x_coord - width / 2, 0),
                    width,
                    1, color=c_dict[self.sample_uid_to_sample_name[sample_uid][0]]))
                color_list.append(c_dict[self.sample_uid_to_sample_name[sample_uid][0]])
            else:
                # negative sample
                rectangles.append(Rectangle(
                    (x_coord - width / 2, 0),
                    width,
                    1, color='black'))
                color_list.append('black')
        # patches_collection = PatchCollection(self.bar_patches, match_original=True)
        listed_color_map = ListedColormap(color_list)
        patches_collection = PatchCollection(rectangles, cmap=listed_color_map)
        patches_collection.set_array(np.arange(len(rectangles)))
        self.species_ax.add_collection(patches_collection)
        self.species_ax.set_xlim((x_coords[0] - width, x_coords[-1] + width))
        # Remove the axis
        self.species_ax.set_xticks([])
        self.species_ax.set_yticks([])
        self.species_ax.set_ylabel('species', rotation='vertical', fontsize='x-small')


class BuitragoBars(Buitrago):
    def __init__(self):
        super().__init__()
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

BuitragoHier()
# BuitragoBars()


# # Get the x coordinates of the leaves
        # # Get the point coordinates of the branches that come from the leaftips.
        # # I have saved this code but its basically pointless.
        # x_coords = []
        # leaf_lines = []
        # for coord_list in sph.ax.collections[0].get_segments():
        #     # The x coordinates of leafs
        #     leaf_x_coords = [_[0] for _ in coord_list if _[1] == 0]
        #     x_coords.extend(leaf_x_coords)
        #     # The y coordinates of the first line originating from the leaf base
        #     leaf_y_coords = [_[1] for _ in coord_list if _[0] in leaf_x_coords and _[1] != 0]
        #     leaf_lines.extend([[(x, 0), (x, y)] for x, y in zip(leaf_x_coords, leaf_y_coords)])
        # x_coords = sorted(x_coords)
        # leaf_lines.sort(key=lambda x: x[0][0])