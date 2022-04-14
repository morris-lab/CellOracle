# -*- coding: utf-8 -*-


import os, sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from .config import CONFIG

from .development_module_visualization import (plot_cluster_whole,
                            plot_quiver, plot_quiver_random,
                            plot_simulation_flow_on_grid,
                            plot_simulation_flow_random_on_grid)

class Oracle_visualization:
    """docstring for ."""

    def __init__(self):
        pass

    def plot_cluster_whole(self, ax=None, s=CONFIG["s_scatter"], args=CONFIG["default_args"]):
        plot_cluster_whole(self=self, ax=ax, s=s, args=args)

    def plot_quiver(self, ax=None, scale=CONFIG["scale_simulation"], color=None, s=CONFIG["s_scatter"], show_background=True, args=CONFIG["default_args"]):
        plot_quiver(self=self, ax=ax, scale=scale, color=color, s=s, show_background=show_background, args=args)

    def plot_quiver_random(self, ax=None, scale=CONFIG["scale_simulation"], color=None, s=CONFIG["s_scatter"], show_background=True, args=CONFIG["default_args"]):
        plot_quiver_random(self=self, ax=ax, scale=scale, color=color, s=s, show_background=show_background, args=args)

    def plot_simulation_flow_on_grid(self, ax=None, scale=CONFIG["scale_simulation"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
        plot_simulation_flow_on_grid(self=self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)

    def plot_simulation_flow_random_on_grid(self, ax=None, scale=CONFIG["scale_simulation"], show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
        plot_simulation_flow_random_on_grid(self=self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)
