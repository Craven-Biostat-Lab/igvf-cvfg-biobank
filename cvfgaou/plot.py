import matplotlib.pyplot as plt
from matplotlib import colors, ticker
import numpy as np

# Color codes
clinvar_colors = {
    'Benign': '#1D7AAB',
    'Likely benign': '#64A1C4',
    'Uncertain significance': '#A0A0A0',
    'Conflicting': 'grey',
    'Likely pathogenic': '#E6B1B8',
    'Pathogenic': '#CA7682',
}

cohort_colors = {
    'cases': '#CA7682',
    'overlap': 'grey',
    'controls': '#1D7AAB'
}

points_colors = {
    '≤ -16': '#176082',
    '≤ -12': '#176082',
    '≤ -11': '#1D7AAB',
    '≤ -10': '#4B91A6',
    '≤ -9': '#63A1C4',
    '≤ -8': '#7AB5D1',
    '≤ -7': '#99C8DC',
    '≤ -6': '#B8DCE8',
    '≤ -5': '#D0E8F0',
    '≤ -4': '#E4F1F6',
    '≤ -3': '#EDF6FA',
    '≤ -2': '#F4F9FC',
    '≤ -1': '#F9FCFE',
    '≥ +1': '#F7E4E7',
    '≥ +2': '#F0D0D5',
    '≥ +3': '#E6B1B8',
    '≥ +4': '#D68F99',
    '≥ +5': '#CA7682',
    '≥ +6': '#B85C6B',
    '≥ +7': '#A84957',
    '≥ +8': '#943744',
    '≥ +9': '#7F2936',
    '≥ +10': '#671B28',
    '≥ +11': '#520F1C',
    '≥ +12': '#3A060D',
    '≥ +16': '#3A060D'
}

# Plotters

def grouped_bar_plot(ax, groups, y_values, bar_colors=None, xlabel=None, log_scale=False, **style_args):
    # Grouped bar positioning:
    # If we have n bars in a group, we want the width of each bar to be
    # w=1/(n+1) (to leave a bar's width worth of space between groups.
    # The first bar is then positioned at tick-w*(n-1)/2; this way, the
    # n-th bar will be at tick-w*(n-1)/2 + w*(n-1) = tick+w*(n-1)/2;
    # The resulting expression for the i-th (0-indexed) bar is
    # tick+w*(i - (n - 1)/2)

    nbars=len(groups)
    bar_size=1/(nbars+1)
    
    if bar_colors is None:
        bar_colors = list(
            list(plt.rcParams['axes.prop_cycle'])[n]['color']
            for n in range(nbars)
        )
    else:
        bar_colors = list(bar_colors)

    for n, (label, series) in enumerate(groups):

        if log_scale:
            ax.set_axisbelow(True)
            ax.set_xscale('log')
            
        color = bar_colors[n]

        bc = ax.barh(
            y=y_values - bar_size * (n - (nbars - 1)/2),
            width=series,
            height=bar_size,
            label=label,
            facecolor=colors.to_rgba(color, 0.2),
            edgecolor=color,
            **style_args
        )
        #ax.bar_label(bc, label_type='edge')
        #left += series

    if xlabel is not None:
        ax.set_xlabel(xlabel)

    #ax.legend()

def populate_classifier_row(subfigs, row, classifier, plot_df):
    
    subfigs[row].suptitle(classifier)
    
    axs = subfigs[row].subplots(1, 4, sharey=True)
    
    for ax in axs:
        ax.grid(True, which='major', axis='x', linestyle='--')
        ax.grid(True, which='minor', axis='x', linestyle=':')
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    if row > 0:
        # Share x across OR plots
        axs[0].sharex(subfigs[row-1].axes[0])
        #for col in range(3):
        #    axs[col].sharex(subfigs[row-1].axes[col])

    # Establish Y positions
    y_labels = plot_df.Classification.cat.remove_unused_categories().cat.categories
    y_map = {
        label: index
        for index, label in enumerate(y_labels)
    }

    ##################
    # Log OR subplot #
    ##################

    axs[0].errorbar(
        x=plot_df.LogOR,
        y=plot_df.Classification.map(y_map),
        xerr=np.absolute(plot_df[['LogOR_LI', 'LogOR_UI']].to_numpy() - plot_df[['LogOR']].to_numpy()).transpose(),
        fmt='.',
        capsize=4,
        capthick=1.5,
        elinewidth=1,
        color='black'
    )

    axs[0].axvline(x=0, color='red', linestyle=':')

    axs[0].set_xlabel("Log OR")

    axs[0].set_yticks(range(len(y_labels)), y_labels)
    
    #########################
    # Variant count subplot #
    #########################

    grouped_bar_plot(
        axs[1],
        [
            ('cases', plot_df.case_only_variant_count),
            ('overlap', plot_df.overlap_variant_count),
            ('controls', plot_df.control_only_variant_count)
        ],
        plot_df.Classification.map(y_map),
        bar_colors=(cohort_colors[c] for c in ('cases', 'overlap', 'controls')),
        xlabel='Variant count'
    )

    cohort_handles, cohort_labels = axs[1].get_legend_handles_labels()
    subfigs[row].legend(
        cohort_handles,
        cohort_labels,
        loc='outside upper right',
        title='Cohort'
    )
    
    #############################
    # Patricipant count subplot #
    #############################

    grouped_bar_plot(
        axs[2],
        [
            ('cases', plot_df.cases_with_variants),
            ('controls', plot_df.controls_with_variants)
        ],
        plot_df.Classification.map(y_map),
        bar_colors=(cohort_colors[c] for c in ('cases', 'controls')),
        xlabel='Participant count',
        log_scale=True
    )
    
    #############################
    # Clinvar breakdown subplot #
    #############################
    
    clinvar_subplot_title = 'ClinVar significance'
    
    class_clinvar_bars, cohort_clinvar_bars = (
        [
            (significance, plot_df[f'{prefix} {significance}'])
            for significance in (
                'Benign',
                'Likely benign',
                'Uncertain significance',
                'Conflicting',
                'Likely pathogenic',
                'Pathogenic',
                #'Other / not in ClinVar'
            )
        ]
        for prefix in ('Class ClinVar', 'Cohort ClinVar')
    )
    clinvar_bar_colors = [
        clinvar_colors[significance]
        for significance in (
            'Benign',
            'Likely benign',
            'Uncertain significance',
            'Conflicting',
            'Likely pathogenic',
            'Pathogenic',
            #'Other / not in ClinVar'
        )
    ]
    #grouped_bar_plot(
    #    axs[3],
    #    class_clinvar_bars,
    #    plot_df.Classification.map(y_map),
    #    #log_scale=True,
    #    xlabel='ClinVar significance'
    #)
    grouped_bar_plot(
        axs[3],
        cohort_clinvar_bars,
        plot_df.Classification.map(y_map),
        bar_colors=clinvar_bar_colors,
        #hatch='//'
        xlabel='ClinVar significance'
    )
    clinvar_handles, clinvar_labels = axs[3].get_legend_handles_labels()
    subfigs[row].legend(
        clinvar_handles,
        clinvar_labels,
        loc='outside lower right',
        title='ClinVar significance'
    )

    #handles, labels = [], []
    #for ax in axs:
    #    ax_handles, ax_labels = ax.get_legend_handles_labels()
    #    handles += ax_handles
    #    labels += ax_labels
    #
    #subfigs[row].legend(handles, labels, loc='outside right')

def summary_fig(
    or_estimates_df,
    combined_fig = None,
    cols='Gene',
    saturation_threshold=0.7,
    frame_on=False,
    left_ticks=False,
    title_rotation=0
):

    if combined_fig is None:
        combined_fig = plt.figure(
            layout='constrained',
            figsize=(15, 10)
        )

    plot_df = or_estimates_df[
        (or_estimates_df.cases_with_variants > 0)
    ]

    col_dfs = [
        (col_id, col_df)
        for col_id, col_df in plot_df.groupby(cols)
        if not col_df.empty
    ]

    axs = combined_fig.subplots(
        1, #combined_or_estimates_df['Gene'].nunique(),
        len(col_dfs),
        sharex = True,
        sharey = True
    )

    class_tick_list = plot_df.Classification.cat.categories.to_list()
    class_tick_map = {val: ind for ind, val in enumerate(class_tick_list)}

    for col, (gene, col_df) in enumerate(col_dfs):

        # Figure out which bars will get a dark background
        dark_rows=(
            (
                col_df.LogOR_LI > 0
            ) & (
                col_df['Classification']
                .astype(str) # The categorical wrapper breaks things, so convert to str
                .map(points_colors)
                .map(colors.to_rgb)
                .map(colors.rgb_to_hsv)
                .apply(lambda c: c[2] < saturation_threshold)
            )
        )
        
        # Need to draw separately for dark and light background rows
        for dark in True, False:
            subset_df = col_df[dark_rows == dark]
            color = 'white' if dark else 'black'
            axs[col].errorbar(
                x=subset_df.LogOR,
                y=subset_df.Classification.map(class_tick_map),
                xerr=np.absolute(
                    subset_df[['LogOR_LI', 'LogOR_UI']].to_numpy()
                    - subset_df[['LogOR']].to_numpy()
                ).transpose(),
                fmt='.',
                capsize=4,
                capthick=1.5,
                elinewidth=1,
                color=color
            )

        axs[col].axvline(x=0, color='black', linestyle=':')
        axs[col].set_yticks(range(len(class_tick_list)), class_tick_list)
        axs[col].set_title(gene, rotation=title_rotation)

        # Highlight rows with significant ORs
        for classification in col_df.loc[col_df.LogOR_LI > 0, 'Classification']:
            y_val = class_tick_map[classification]
            axs[col].axhspan(y_val-0.45, y_val+0.45, color=points_colors[classification])

        axs[col].set_frame_on(frame_on)
        axs[col].tick_params(left=left_ticks)

    return combined_fig