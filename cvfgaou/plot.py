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