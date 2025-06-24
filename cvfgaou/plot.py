import matplotlib.pyplot as plt
import numpy as np

def grouped_bar_plot(ax, groups, y_values, xlabel=None, log_scale=False):
    # Grouped bar positioning:
    # If we have n bars in a group, we want the width of each bar to be
    # w=1/(n+1) (to leave a bar's width worth of space between groups.
    # The first bar is then positioned at tick-w*(n-1)/2; this way, the
    # n-th bar will be at tick-w*(n-1)/2 + w*(n-1) = tick+w*(n-1)/2;
    # The resulting expression for the i-th (0-indexed) bar is
    # tick+w*(i - (n - 1)/2)

    nbars=len(groups)
    bar_size=1/(nbars+1)

    for n, (label, series) in enumerate(groups):

        if log_scale:
            ax.set_xscale('log')

        bc = ax.barh(
            y=y_values - bar_size * (n - (nbars - 1)/2),
            width=series,
            height=bar_size,
            label=label
        )
        ax.bar_label(bc, label_type='edge')
        #left += series

    if xlabel is not None:
        ax.set_xlabel(xlabel)

    ax.legend()

def populate_classifier_row(subfigs, row, classifier, plot_df):
    
    subfigs[row].suptitle(classifier)
    
    axs = subfigs[row].subplots(1, 3, sharey=True)

    if row > 0:
        for col in range(3):
            axs[col].sharex(subfigs[row-1].axes[col])

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
        xlabel='Variant count'
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
        xlabel='Participant count',
        log_scale=True
    )