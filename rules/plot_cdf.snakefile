"""
Make cumulative read size distribution.
"""

# subread_stats_plot_cdf
#
# Make CDF plot.
rule subread_stats_plot_cdf:
    input:
        tab=lambda wildcards: [
            '{0}/cells/{1}/zmw_summary.tab.gz'.format(SAMPLE_NAME, cell)
            for cell in get_cell_dict().keys()
        ]
    output:
        pdf='{sample}/plot/cdf/cdf_{per_cell}_{read_type}.pdf',
        png='{sample}/plot/cdf/cdf_{per_cell}_{read_type}.png'
    wildcard_constraints:
        per_cell='cell|sample',
        read_type='subread|longest'
    run:

        # Get stat list
        size_list = substatlib.cdf.get_size_list(
            wildcards.sample,
            cells=None,
            read_type=wildcards.read_type,
            per_cell=(wildcards.per_cell == 'cell')
        )

        # Get figure
        fig = substatlib.cdf.get_cdf_plot(size_list)

        # Write
        fig.savefig(output.pdf, bbox_inches='tight')
        fig.savefig(output.png, bbox_inches='tight')
