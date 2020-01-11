"""
Compute read-length CDF distributions and create plots.
"""

import numpy as np
import os
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt


#
# Make Plot
#

def get_size_list(sample, cells=None, read_type='subread', per_cell=False):
    """
    Return a list of tuples where the first element is a label and the second
    is an array of read sizes associated with that label.

    If `per_cell` is `True`, then the list has one entry per cell, and the label
    is the cell name. If `per_cell` is `False`, then the list has one entry for
    all reads in the sample, and the label is the sample name.

    :param sample: Sample name.
    :param cells: List of cells to process. If `None`, get all cells for the sample.
    :param read_type: "subread" to plot all subreads or "longest" to plot the
        longest subread per read.
    :param per_cell: Plot one distribution per cell.
    """

    # Get cell list
    if cells is None:
        cells = os.listdir(os.path.join(sample, 'cells'))

    # Check all cell files before reading (takes a long time to read, fail early)
    cells = [cell.strip() if cell is not None else '' for cell in cells]

    for cell in cells:
        if not cell:
            raise RuntimeError('Empty cell name in cell list')

        file_name = os.path.join(sample, 'cells', cell, 'zmw_summary.tab.gz')

        if not os.path.isfile(file_name):
            raise RuntimeError('Missing cell file: {}'.format(file_name))

    # Get size array list
    size_array_list = get_cell_sizes(sample, cells, read_type)

    if per_cell:
        size_list = list(zip(cells, size_array_list))
    else:
        size_list = [(sample, np.concatenate(size_array_list))]

    # Return CDF list
    return size_list


def get_cdf_plot(
        size_list, width=7, height=7, dpi=300, z_cut=2.5, legend=True,
        grid_spec = {'color': '#b8b8c8', 'which': 'both', 'lw': 0.5}
    ):
    """
    Get a CDF plot.

    :param size_list: A list of 2-element tuples, [0] sample name, [1] a numpy array of sizes. Sizes do not need to
        be sorted.
    :param width: Figure width.
    :param height: Figure height.
    :param dpi: Figure resolution.
    :param z_cut: Cut high values at this z-score cutoff. This prevents large outliers from expanding the x-axis
        and compressing the useful visualazation to the left side of x.
    :param legend: Display a legend with line colors and sample names.
    :param grid_spec: A dictionary of grid keywords (for axes.grid) to specify how the grid should be constructed. Set
        to `None` to disable the grid and show a blank background.
    """

    # Make figure
    #fig, ax = plt.figure(1, figsize=(width, height), dpi=dpi)
    # ax = fig.add_subplot(1, 1, 1)

    fig, ax = plt.subplots(1, 1, figsize=(width, height), dpi=dpi)

    for label, sizes in size_list:

        sizes = np.sort(sizes) / 1e3  # Sort sizes and convert to kbp
        cdist = np.flip(np.cumsum(np.flip(sizes))) / 1e6  # Cumulative distribution in Gbp (sizes already in kbp)

        # Trim distribution (remove skew from ultra-longreads)
        if z_cut is not None:
            max_index = np.sum(((sizes - np.mean(sizes)) / np.std(sizes)) > z_cut)

            sizes = sizes[:-max_index]
            cdist = cdist[:-max_index]

        # Add to plot
        ax.plot(sizes, cdist, '-', label=label)

    # Add labels
    ax.set_ylabel('Cumulative throughput (Gbp)', fontsize=12)
    ax.set_xlabel('Read size (kbp)', fontsize=12)

    ax.get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    ax.get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',.2f'))
    )

    # Add legend
    if legend:
        ax.legend(prop={'size': 10})

    # Add grid
    if grid_spec is not None:
        ax.grid(**grid_spec)

    # Return plot
    return fig


def get_cell_sizes(sample, cells, read_type='subread'):
    """
    Get a list of Numpy arrays where each array is the read sizes from one cell.

    :param sample: Sample name.
    :param cells: List of cells to process. If `None`, get all cells for the sample.
    :param read_type: "subread" to plot all subreads or "longest" to plot the
        longest subread per read.
    :param cells: List of cells to process. If `None`, get all cells for the sample.
    """

    # Check sample
    if sample is None:
        raise RuntimeError('Sample is None')

    sample = sample.strip()

    if not sample:
        raise RuntimeError()

    # Check type
    read_type = read_type.lower().strip()

    if read_type not in {'subread', 'longest'}:
        raise RuntimeError('Read type must be "subread" or "longest": {}'.format(read_type))

    # Read
    size_array_list = list()

    for cell in cells:
        file_name = os.path.join(sample, 'cells', cell, 'zmw_summary.tab.gz')

        if not os.path.isfile(file_name):
            raise RuntimeError('Missing cell file in get_cell_sizes(): {}'.format(file_name))

        df = pd.read_csv(file_name, sep='\t', header=0)

        # Get size array
        if read_type == 'subread':
            size_array_list.append(np.asarray(
                [val for val_list in df['LIST'] for val in val_list.split(',')],
                np.int32
            ))

        elif read_type == 'longest':
            size_array_list.append(np.asarray(df['MAX'], np.int32))

        else:
            raise RuntimeError('BUG: Unknown read type {}'.format(read_type))

    # Return list of sizes
    return size_array_list
