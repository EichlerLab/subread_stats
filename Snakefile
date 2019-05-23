"""
Get stats for raw sequencing data in each sample.
"""

import collections
import os
import pysam

import pandas as pd
import numpy as np

localrules: subread_stats_get_tables


###############
# Definitions #
###############

# Get config
FOFN_FILE_NAME = config.get('fofn', None)

SAMPLE_NAME = config.get('sample', None)

if FOFN_FILE_NAME is None:
    raise RuntimeError('Missing FOFN file in config (list of input PacBio BAM or BAX files). Example: snakemake ... --config fofn=/path/to/sample.fofn')

if SAMPLE_NAME is None:
    if FOFN_FILE_NAME.lower().endswith('.fofn'):
        SAMPLE_NAME = os.path.basename(FOFN_FILE_NAME).rsplit('.', 1)[0]
    else:
        raise RuntimeError(
            'Sample name must be set (e.g. snakemake ... --config fofn=/path/to/sample.fofn sample=NAME). If the '
            'sample name is not set, the FOFN file name must end in .fofn, and the sample name will be the FOFN file '
            'name without the .fofn extension.'
        )


#######################
### Read FOFN Table ###
#######################


def get_cell_dict():
    """
    Get a list of input files per cell. Returns a dictionary where keys are cell names and values are a list of
    input BAM or BAX files for that cell.
    """

    with open(FOFN_FILE_NAME, 'r') as in_file:

        cell_dict = collections.defaultdict(list)

        line_count = 0

        for line in in_file:

            line_count += 1

            line = line.strip()

            if not line or line.startswith('#'):
                continue

            if not line.lower().endswith('.subreads.bam'):
                raise RuntimeError(
                    'Unrecognized file type in FOFN on line {}: Must end with ".subreads.bam": {}'.format(
                        line_count,
                        line
                    )
                )

            cell_name = os.path.basename(line).rsplit('.', 2)[0]

            if cell_name in cell_dict:
                raise RuntimeError('Duplicate cell on line {}: {} ("{}" conflicts with previous entry "{}")'.format(
                    line_count,
                    cell_name,
                    os.path.basename(line),
                    os.path.basename(cell_dict[cell_name])
                ))

            cell_dict[cell_name] = line

    return cell_dict


#########
# Rules #
#########

# subread_stats_get_tables
#
# Get all tables.
rule subread_stats_get_tables:
    input:
        tab_sample='tables/samples/{0}/sample_summary.tab'.format(SAMPLE_NAME),
        tab_cell='tables/samples/{0}/cell_summary.tab'.format(SAMPLE_NAME)

# subread_stats_summarize_sample
#
# Summarize stats over a sample.
rule subread_stats_summarize_sample:
    input:
        tab=expand('tables/samples/{0}/cells/{{CELL}}/zmw_summary.tab.gz'.format(SAMPLE_NAME), CELL=get_cell_dict().keys())
    output:
        tab='tables/samples/{0}/sample_summary.tab'.format(SAMPLE_NAME),
        xlsx='tables/samples/{0}/sample_summary.xlsx'.format(SAMPLE_NAME)
    run:

        # Get a list of all subread lengths
        subread_list = [
            int(val) for tab_file in input.tab for sublist in pd.read_csv(tab_file, sep='\t')['LIST'].apply(
                lambda vals: vals.split(',')
            ) for val in sublist
        ]

        subread_list = sorted(subread_list)

        # Get a list of longest subreads
        roi_list = [
            val for tab_file in input.tab for val in pd.read_csv(tab_file, sep='\t')['MAX']
        ]

        roi_list = sorted(roi_list)

        # Get stats
        df_summary = pd.DataFrame(
            {
                'SAMPLE': SAMPLE_NAME,
                'CELLS': len(input.tab),
                'SUB_MEAN': pd.Series(np.mean(subread_list)),
                'SUB_MED': np.median(subread_list),
                'SUB_N50': subread_list[sum(np.cumsum(subread_list) < np.sum(subread_list) * 0.5)],
                'SUB_SD': np.std(subread_list),
                'SUB_SUM': np.sum(subread_list),
                'SUB_N': len(subread_list),
                'ROI_MEAN': pd.Series(np.mean(roi_list)),
                'ROI_MED': np.median(roi_list),
                'ROI_N50': roi_list[sum(np.cumsum(roi_list) < np.sum(roi_list) * 0.5)],
                'ROI_SD': np.std(roi_list),
                'ROI_SUM': np.sum(roi_list),
                'ROI_N': len(roi_list),
                'MAX': max(subread_list)
            }
        )

        # Write
        df_summary.to_csv(output.tab, sep='\t', index=False)
        df_summary.to_excel(output.xlsx, index=False)


# subread_stats_get_all_subread_bam
#
# Merge all all cell sumaries for this sample.
rule subread_stats_merge_subread_stats:
    input:
        tab=lambda wildcards: [
            'tables/samples/{0}/cells/{1}/cell_summary.tab'.format(SAMPLE_NAME, cell)
            for cell in get_cell_dict().keys()
        ]
    output:
        tab='tables/samples/{0}/cell_summary.tab'.format(SAMPLE_NAME),
        xlsx='tables/samples/{0}/cell_summary.xlsx'.format(SAMPLE_NAME)
    run:

        # Merge and write
        df = pd.concat(
            [pd.read_csv(tab_file, sep='\t', header=0) for tab_file in input.tab]
        )

        df.to_csv(output.tab, sep='\t', index=False)
        df.to_excel(output.xlsx, index=False)

# subread_stats_get_subread_stats
#
# Get stats per subread.
rule subread_stats_get_subread_stats:
    output:
        tab_summary=protected('tables/samples/{0}/cells/{{cell}}/cell_summary.tab'.format(SAMPLE_NAME)),
        tab_zmw=protected('tables/samples/{0}/cells/{{cell}}/zmw_summary.tab.gz'.format(SAMPLE_NAME))
    run:

        # Get subread file
        seq_file = get_cell_dict().get(wildcards.cell, None)

        if seq_file is None:
            raise RuntimeError('No sequence data file for cell {}'.format(wildcards.cell))

        # Init
        zmw_list = list()
        len_list = list()

        # Read each line in BAM
        with pysam.AlignmentFile(seq_file, 'rb', check_sq=False) as sam_file:
            for record in sam_file:

                # Record name format: cell/zmw/start_end (start/end 0-based, half-open, BED format)
                cell, zmw, coord = record.query_name.split('/')
                start, stop = [int(val) for val in coord.split('_')]

                zmw_list.append(int(zmw))
                len_list.append(stop - start)

        # Make into series (key = zmw, value = length)
        ser = pd.Series(len_list, zmw_list)
        ser.index.name = 'ZMW'
        ser.name = 'LEN'

        del(len_list, zmw_list)

        # Group by ZMW and get stats
        df = ser.groupby(ser.index).agg([len, np.max, np.sum, lambda l: ','.join(l.astype(str))])
        df.columns = ('N', 'MAX', 'SUM', 'LIST')

        # Write stats per ZMW
        df.to_csv(output.tab_zmw, sep='\t', compression='gzip', index=True)

        # Sort values for N50 calculations
        df.sort_values('MAX', inplace=True)
        ser.sort_values(inplace=True)

        # Summarize by cell
        df_summary = pd.DataFrame(
            {
                'SUB_MEAN': pd.Series(np.mean(ser)),
                'SUB_MED': pd.Series(np.median(ser)),
                'SUB_N50': np.min(ser[np.cumsum(ser) >= np.sum(ser) * 0.5]),
                'SUB_SD': np.std(ser),
                'SUB_SUM': np.sum(ser),
                'SUB_N': ser.shape[0],
                'ROI_MEAN': np.mean(df['MAX']),
                'ROI_MED': np.median(df['MAX']),
                'ROI_N50': np.min(df['MAX'][np.cumsum(df['MAX']) >= np.sum(df['MAX']) * 0.5]),
                'ROI_SD': np.std(df['MAX']),
                'ROI_SUM': np.sum(df['MAX']),
                'ROI_N': df.shape[0],
                'MAX': max(ser)
            }
        )

        df_summary = df_summary.loc[: , (
            'SUB_MEAN', 'SUB_MED', 'SUB_SD', 'SUB_SUM', 'SUB_N', 'SUB_N50',
            'ROI_MEAN', 'ROI_MED', 'ROI_SD', 'ROI_SUM', 'ROI_N', 'ROI_N50',
            'MAX'
        )]

        df_summary.index = [wildcards.cell]
        df_summary.index.name = 'CELL'

        # Write
        df_summary.to_csv(output.tab_summary, sep='\t', index=True)
