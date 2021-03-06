#!/usr/bin/env python

# ==============================================================================
# plots.py
#
# by Adam Freese <maxwellsdemon137@gmail.com> copyright 2017
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This file is part of dglap90.
#
# dglap90 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dglap90 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with dglap90.  If not, see <http://www.gnu.org/licenses/>.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# This script makes a plot of the u, d and s distributions at leading order
# at two Q**2 values, and compares the result of evolving from one to the other.
# The high and low Q**2 PDFs are from CJ15.
# The user should run the tester program to produce the evolved PDF before
# running this script.

# Modules used
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import pandas

# Fix font
mpl.rc('font',size=16,family='cmr10',weight='normal')
mpl.rc('text',usetex=True)

# Dicts for making the plots
parton_colors = {
        'u' : '#ff7f00',
        'd' : '#377eb8',
        's' : '#984ea3',
        }
lines = {
        'low'     : ':',
        'evolved' : '--',
        'high'    : '-'
        }

# Files containing dataframes
files = {
        'low'     : 'pdf_grids/cj15_Q2_5_GeV2.csv',
        'evolved' : 'pdf_grids/evolved_Q2_1000_GeV2.csv',
        'high'    : 'pdf_grids/cj15_Q2_1000_GeV2.csv'
        }
labels = {
        'low'     : r'LO CJ15 at $Q^2=5$~GeV$^2$',
        'evolved' : r'LO DGLAP to $Q^2=1000$~GeV$^2$',
        'high'    : r'LO CJ15 at $Q^2=1000$~GeV$^2$'
        }

# Method to plot PDFs from a single file
def plot_pdfs(key):
    # Prepare data frame
    df = pandas.read_csv(files[key])
    # Loop over partons
    for parton in parton_colors:
        pdf = df['x'] * df[parton]
        # Apply labels on up quark
        label = ""
        if parton=='u':
            label = labels[key]
        plt.plot( df['x'], pdf, lines[key], color=parton_colors[parton], label=label )

# Prepare canvas
plt.semilogx()
plt.xlabel(r'$x$')
plt.ylabel(r'$xq(x)$')

# Plot the PDFs
for key in files:
    plot_pdfs( key )

# Legend
plt.legend( loc=1, prop={'size':12} )

# Save the file
plt.tight_layout()
plt.savefig("test_plot.pdf")

# Destroy the plot to free up memory
plt.close()
