import os
from collections import defaultdict
from src.utils import constants as c
import pandas as pd
import datetime
import numpy as np
import math
import re
from jinja2 import Template
from statistics import median

import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as pltx
from matplotlib.colors import rgb2hex
from matplotlib.colors import Normalize

from version import __version__

from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.models.widgets import TableColumn, Slider, Div, Select, Panel, Tabs, CheckboxGroup
from bokeh.models import (DataTable,
                          CustomJS)
from bokeh.embed import components
from bokeh.resources import INLINE
from bokeh.plotting import ColumnDataSource


SUM_OF_SCORES = 'Sum of scores'


def create_heatmap_bar(output_dir):
    fig = pltx.figure(figsize=(2, 2))
    x = np.arange(25).reshape(5, 5)
    cm = sns.diverging_palette(12, 250, sep=30, l=50, n=15, s=90, as_cmap=True)
    ax = sns.heatmap(x, cmap=cm, cbar=False)
    cbar = pltx.colorbar(ax.get_children()[0], orientation='horizontal')

    # hide border
    cbar.outline.set_visible(False)

    # hide ticks
    cbar.set_ticks([])

    # hide heatmap
    pltx.gca().set_visible(False)

    fig.savefig(os.path.join(output_dir, 'heatmap_bar.png'), dpi=100, format='png', bbox_inches='tight', pad_inches=-.035, transparent=True)
    pltx.close(fig)


def create_title_div(id, name, info):
    div = Div(text="""<div style="text-align:left;font-size: 20pt;font-weight: bold;">{1}<span style="float: right;font-size: 10pt;font-weight:normal;">{2}</span>""".format(id, name, info), css_classes=['bk-width-auto']) # width=DIV_WIDTH, height=DIV_HEIGHT)
    return div


def get_columns(labels, current_columns):
        columns = [c.GS]
        for label in labels:
            if label in current_columns:
                columns.append(label)
        return columns


def get_rank_to_sample_pd(pd_metrics):
    rank_to_sample_pd = defaultdict(dict)
    pd_copy = pd_metrics.copy()

    # label ranks, so that the rows don't get lost
    pd_copy.loc[pd_copy[pd.isnull(pd_copy['rank'])].index, 'rank'] = 'rank independent'

    pd_copy = pd_copy[-pd_copy['metric'].isin(['time', 'memory'])]

    # transform table
    pd_grouped = pd_copy.pivot_table(index=['rank', 'tool', 'sample'], columns='metric', values='value')

    # add sample "(average over samples)"
    pd_mean_over_samples = pd_grouped.groupby(['rank', 'tool'], sort=False).mean()
    # convert means to string to store '<mean> (<standard error>)'
    pd_mean_over_samples_str = pd_mean_over_samples.astype(str)
    pd_sem_over_samples = pd_grouped.groupby(['rank', 'tool'], sort=False).sem()

    # add standard error to "(average over samples)"
    df_columns = pd_sem_over_samples.columns.values
    for index, row in pd_sem_over_samples.iterrows():
        for j, item in enumerate(row):
            sem_value = pd_sem_over_samples.at[index, df_columns[j]]
            mean_value = '{:>.{precision}g}'.format(pd_mean_over_samples.at[index, df_columns[j]], precision=3)
            # if standard error is nan because there is only one sample, then just copy mean
            if np.isnan(sem_value):
                pd_mean_over_samples_str.at[index, df_columns[j]] = mean_value
            else:
                sem_value = '{:>.{precision}g}'.format(sem_value, precision=3)
                sem_value = '<div class="tooltip tooltip-sem">{}<span class="tooltiptext">standard error: {}</span></div>'.format(sem_value, sem_value)
                pd_mean_over_samples_str.at[index, df_columns[j]] = "{} ({})".format(mean_value, sem_value)
    pd_mean_over_samples_str['sample'] = '(average over samples)'
    pd_mean_over_samples_str = pd_mean_over_samples_str.reset_index().set_index(['rank', 'tool', 'sample'])

    pd_grouped = pd.concat([pd_grouped, pd_mean_over_samples_str])

    # copy unifrac values to every taxonomic rank
    pd_grouped_copy = pd_grouped.copy()
    for index, row in pd_grouped_copy.iterrows():
        pd_grouped.loc[index][c.UNIFRAC] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNIFRAC]
        pd_grouped.loc[index][c.UNW_UNIFRAC] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNW_UNIFRAC]

    for (rank, sample), g in pd_grouped.groupby(['rank', 'sample']):
        rank_to_sample_pd[rank][sample] = g.reset_index().rename(columns={'tool': 'Tool'}).drop(['rank', 'sample'], axis=1).set_index('Tool').T
        # drop all metrics except unifrac for 'rank independent'
        if rank == 'rank independent':
            rank_to_sample_pd[rank][sample] = rank_to_sample_pd[rank][sample].drop(c.ALL_METRICS[2:])
    return rank_to_sample_pd


def get_formatted_pd_rankings(pd_rankings):
    df_list = []
    df_list_unsorted_pos = []
    metrics_list = []
    for metric, g in pd_rankings.groupby(level=0):
        metrics_list.append(metric)
        df = g.reset_index().sort_values('position')
        df2 = g.reset_index()
        df_list.append(pd.DataFrame({metric: df['tool'].tolist(), 'score' + metric: df['position'].tolist()}))
        df_list_unsorted_pos.append(pd.DataFrame({metric: df2['tool'].tolist(), 'score' + metric: df2['position'].tolist()}))

    df_sum = pd_rankings.groupby(['tool'])['position'].sum().reset_index().sort_values('position')
    df_sum_unsorted_pos = pd_rankings.groupby(['tool'])['position'].sum().reset_index()
    df_list.append(
        pd.DataFrame({SUM_OF_SCORES: df_sum['tool'].tolist(), 'score' + SUM_OF_SCORES: df_sum['position'].tolist()}))
    df_list_unsorted_pos.append(
        pd.DataFrame({SUM_OF_SCORES: df_sum_unsorted_pos['tool'].tolist(), 'score' + SUM_OF_SCORES: df_sum_unsorted_pos['position'].tolist()}))

    pd_show = pd.concat(df_list, axis=1)
    pd_show_unsorted_pos = pd.concat(df_list_unsorted_pos, axis=1)
    return pd_show, pd_show_unsorted_pos


def create_rankings_html(pd_rankings):
    pd_show, pd_show_unsorted_pos = get_formatted_pd_rankings(pd_rankings)

    table_source = ColumnDataSource(pd_show)

    columns = [
        TableColumn(field=SUM_OF_SCORES, title=SUM_OF_SCORES, sortable=False),
        TableColumn(field='score' + SUM_OF_SCORES, title='', width=50),
        TableColumn(field=c.RECALL, title=c.RECALL, sortable=False),
        TableColumn(field='score' + c.RECALL, title='', width=50),
        TableColumn(field=c.PRECISION, title=c.PRECISION, sortable=False),
        TableColumn(field='score' + c.PRECISION, title='', width=50),
        TableColumn(field=c.L1NORM, title=c.L1NORM, sortable=False),
        TableColumn(field='score' + c.L1NORM, title='', width=50),
        TableColumn(field=c.UNIFRAC, title=c.UNIFRAC, sortable=False),
        TableColumn(field='score' + c.UNIFRAC, title='', width=50),
    ]
    data_table = DataTable(source=table_source, columns=columns, width=800, height=25 + len(pd_show) * 25)

    top = [float(x) for x in pd_show_unsorted_pos['score' + SUM_OF_SCORES]]
    source = ColumnDataSource(data=dict(x=pd_show_unsorted_pos[SUM_OF_SCORES].tolist(),
                                        top=top,
                                        recall=pd_show_unsorted_pos['score' + c.RECALL],
                                        precision=pd_show_unsorted_pos['score' + c.PRECISION],
                                        l1norm=pd_show_unsorted_pos['score' + c.L1NORM],
                                        unifrac=pd_show_unsorted_pos['score' + c.UNIFRAC]))

    callback = CustomJS(args=dict(source=source), code="""
        var data = source.data;
        var wrecall = weight_recall.value;
        var wprecision = weight_precision.value;
        var wl1norm = weight_l1norm.value;
        var wunifrac = weight_unifrac.value;
        topx = data['top'];
        recall = data['recall'];
        precision = data['precision'];
        l1norm = data['l1norm'];
        unifrac = data['unifrac'];
        
        
        for (i = 0; i < topx.length; i++) {
            topx[i] = recall[i] * wrecall + precision[i] * wprecision + l1norm[i] * wl1norm + unifrac[i] * wunifrac;
        }
        
        source.change.emit();
    """)

    weight_recall = Slider(start=0, end=10, value=1, step=.1, title=c.RECALL + " weight", callback=callback)
    callback.args["weight_recall"] = weight_recall

    weight_precision = Slider(start=0, end=10, value=1, step=.1, title=c.PRECISION + " weight", callback=callback)
    callback.args["weight_precision"] = weight_precision

    weight_l1norm = Slider(start=0, end=10, value=1, step=.1, title=c.L1NORM + " weight", callback=callback)
    callback.args["weight_l1norm"] = weight_l1norm

    weight_unifrac = Slider(start=0, end=10, value=1, step=.1, title=c.UNIFRAC + " weight", callback=callback)
    callback.args["weight_unifrac"] = weight_unifrac

    p = figure(x_range=pd_show_unsorted_pos[SUM_OF_SCORES].tolist(), plot_width=800, plot_height=400, title=SUM_OF_SCORES + " - lower is better")
    p.vbar(x='x', top='top', source=source, width=0.5, bottom=0, color="firebrick")

    col_rankings = column([Div(text="<font color='navy'><u>Hint 1:</u> click on the columns of scores for sorting.</font>", style={"width": "500px", "margin-bottom": "10px"}),
                           data_table,
                           Div(text="<font color='navy'><u>Hint 2:</u> slide the bars to change the weight of the metrics.</font>", style={"width": "500px", "margin-top": "18px"}),
                           row(weight_recall, weight_precision),
                           row(weight_l1norm, weight_unifrac),
                           p], css_classes=['bk-padding-top'])
    return col_rankings


def get_colors_and_ranges(name, all_values, df_metrics):
    color1 = 'dodgerblue'
    color2 = 'red'
    hue1 = 12
    hue2 = 240
    if name == c.PRECISION or name == c.RECALL or name == c.F1_SCORE or name == c.JACCARD:
        return color1, color2, hue1, hue2, 0, 1
    if name == c.FP or name == c.UNIFRAC or name == c.UNW_UNIFRAC:
        return color2, color1, hue2, hue1, 0, max(all_values)
    if name == c.TP:
        return color1, color2, hue1, hue2, 0, max(all_values)
    if name == c.FN:
        fn_values = df_metrics.loc[c.TP, ].values
        # convert "<mean> (<standard error>)" to float of <mean>
        if len(fn_values) > 0 and isinstance(fn_values[0], str):
            fn_values = [float(x.split(' ')[0]) for x in fn_values]
        return color2, color1, hue2, hue1, 0, max(fn_values)
    if name == c.L1NORM:
        return color2, color1, hue2, hue1, 0, 2
    if name == c.BRAY_CURTIS:
        return color2, color1, hue2, hue1, 0, 1
    return color1, color2, hue1, hue2, max(all_values), min(all_values)


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


def get_heatmap_colors(pd_series, **args):
    values = pd_series.tolist()

    # convert "<mean> (<standard error>)" to float of <mean>
    if len(values) > 0 and isinstance(values[0], str):
        values = [float(x.split(' ')[0]) for x in values]
    all_values = values[:]

    dropped_gs = False
    if pd_series.index[0] == c.GS:
        pd_series.drop(c.GS)
        values = values[1:]
        dropped_gs = True
    if len(values) == 0:
        return ['']

    notnan_values = [x for x in values if isinstance(x, (float, int)) and not np.isnan(x)]
    notnan_all_values = [x for x in all_values if isinstance(x, (float, int)) and not np.isnan(x)]
    if len(notnan_values) == 0:
        red = 'background-color: red'
        return [red for x in values] if not dropped_gs else [''] + [red for x in values]

    color1, color2, hue1, hue2, min_value, max_value = get_colors_and_ranges(pd_series.name, all_values, args["df_metrics"])

    cm = sns.diverging_palette(hue1, hue2, sep=50, l=80, n=15, s=90, as_cmap=True)
    norm = MidpointNormalize(vmin=min_value, vmax=max_value, midpoint=median(notnan_all_values))

    normed = norm([round(x, 3) if not math.isnan(x) else max_value for x in values])
    heatmap_colors = [rgb2hex(x) for x in pltx.cm.get_cmap(cm)(normed)]
    return_colors = []
    for val, x in zip(values, heatmap_colors):
        if val == min_value:
            return_colors.append('background-color: {}'. format(color2))
        elif val == max_value:
            return_colors.append('background-color: {}'. format(color1))
        elif math.isnan(val):
            return_colors.append('background-color: red')
        else:
            return_colors.append('background-color: {}'. format(x))

    if dropped_gs:
        return [''] + return_colors
    else:
        return return_colors


def create_metrics_table(pd_metrics, labels, sample_ids_list):
    rank_to_sample_pd = get_rank_to_sample_pd(pd_metrics)

    all_sample_ids = sample_ids_list[:]
    all_sample_ids.insert(0, '(average over samples)')

    presence_metrics = [c.RECALL, c.PRECISION, c.F1_SCORE, c.TP, c.FP, c.FN, c.JACCARD]
    estimates_metrics = [c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.BRAY_CURTIS]
    alpha_diversity_metics = [c.OTUS, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT]
    rank_independent_metrics = [c.UNIFRAC, c.UNW_UNIFRAC]
    all_metrics = [presence_metrics, estimates_metrics, alpha_diversity_metics]

    presence_metrics_label = 'Presence/absence of taxa'
    estimates_metrics_label = 'Abundance estimates'
    alpha_diversity_metics = 'Alpha diversity'
    all_metrics_labels = [presence_metrics_label, estimates_metrics_label, alpha_diversity_metics]

    styles = [{'selector': 'td', 'props': [('width', '100pt')]},
              {'selector': 'th', 'props': [('width', '100pt'), ('text-align', 'left')]},
              {'selector': 'th:nth-child(1)', 'props': [('width', '120pt'), ('font-weight', 'normal')]},
              {'selector': '', 'props': [('width', 'max-content'), ('width', '-moz-max-content'), ('border-top', '1px solid lightgray'), ('border-spacing', '0px')]},
              {'selector': 'expand-toggle:checked ~ * .data', 'props': [('background-color', 'white !important')]}]
    styles_hidden_thead = styles + [{'selector': 'thead', 'props': [('display', 'none')]}]

    d = {c.RECALL: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.RECALL, c.RECALL, c.TOOLTIP_RECALL),
         c.PRECISION: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.PRECISION, c.PRECISION, c.TOOLTIP_PRECISION),
         c.F1_SCORE: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.F1_SCORE, c.F1_SCORE, c.TOOLTIP_F1_SCORE),
         c.TP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.TP, c.TP, c.TOOLTIP_TP),
         c.FP: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.FP, c.FP, c.TOOLTIP_FP),
         c.FN: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.FN, c.FN, c.TOOLTIP_FN),
         c.JACCARD: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.JACCARD, c.JACCARD, c.TOOLTIP_JACCARD),
         c.UNIFRAC: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.UNIFRAC, c.UNIFRAC, c.TOOLTIP_UNIFRAC),
         c.UNW_UNIFRAC: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.UNW_UNIFRAC, c.UNW_UNIFRAC, c.TOOLTIP_UNW_UNIFRAC),
         c.L1NORM: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.L1NORM, c.L1NORM, c.TOOLTIP_L1NORM),
         c.BRAY_CURTIS: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.BRAY_CURTIS, c.BRAY_CURTIS, c.TOOLTIP_BRAY_CURTIS),
         c.OTUS: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.OTUS, c.OTUS, c.TOOLTIP_OTUS),
         c.SHANNON_DIVERSITY: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.SHANNON_DIVERSITY, c.SHANNON_DIVERSITY, c.TOOLTIP_SHANNON_DIVERSITY),
         c.SHANNON_EQUIT: '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(c.SHANNON_EQUIT, c.SHANNON_EQUIT, c.TOOLTIP_SHANNON_EQUIT)}
    pattern = re.compile('|'.join(map(re.escape, d)))

    def translate(match):
        return d[match.group(0)]

    rank_to_sample_to_html = defaultdict(list)
    for rank in rank_to_sample_pd:
        for sample_id in all_sample_ids:
            if sample_id not in rank_to_sample_pd[rank]:
                rank_to_sample_to_html[rank].append("")
                continue
            mydf = rank_to_sample_pd[rank][sample_id]
            mydf.index.name = None
            html = ''
            first_metrics = True
            for metrics, metrics_label in zip(all_metrics, all_metrics_labels):
                if rank == 'rank independent':
                    metrics = rank_independent_metrics
                    metrics_label = estimates_metrics_label
                html += '<p style="margin-bottom: auto"><b>{}</b></p>'.format(metrics_label)
                mydf_metrics = mydf.loc[metrics]

                sorted_columns = get_columns(labels, mydf_metrics.columns.tolist())
                mydf_metrics = mydf_metrics.loc[:, sorted_columns]

                if first_metrics:
                    this_style = styles
                else:
                    this_style = styles_hidden_thead
                if metrics_label == presence_metrics_label or metrics_label == estimates_metrics_label:
                    html += mydf_metrics.style.apply(get_heatmap_colors, df_metrics=mydf_metrics, axis=1).set_precision(3).set_table_styles(this_style).render()
                else:
                    html += mydf_metrics.style.set_precision(3).set_table_styles(this_style).render()
                if rank == 'rank independent':
                    break
                first_metrics = False
            html = pattern.sub(translate, html)
            rank_to_sample_to_html[rank].append(html)

    mytable1 = Div(text="""<div>{}</div>""".format(rank_to_sample_to_html[c.ALL_RANKS[0]][0]), css_classes=['bk-width-auto'])

    select_rank = Select(title="Taxonomic rank:", value=c.ALL_RANKS[0], options=c.ALL_RANKS + ['rank independent'], css_classes=['bk-fit-content'])

    select_sample = Select(title="Sample:", value='0', options=list(zip(map(str, range(len(all_sample_ids))), all_sample_ids)), css_classes=['bk-fit-content'])

    source = ColumnDataSource(data=rank_to_sample_to_html)

    select_rank_sample_callback = CustomJS(args=dict(source=source), code="""
        mytable1.text = source.data[select_rank.value][select_sample.value];
    """)
    select_rank.js_on_change('value', select_rank_sample_callback)
    select_sample.js_on_change('value', select_rank_sample_callback)
    select_rank_sample_callback.args["select_rank"] = select_rank
    select_rank_sample_callback.args["select_sample"] = select_sample
    select_rank_sample_callback.args["mytable1"] = mytable1

    heatmap_legend = '<img src="heatmap_bar.png" /><div style="text-align:left;font-size: 11px;">Worst<span style="float:right;">Best</span><span style="margin-right: 36px;float:right;">Median</span></div>'
    heatmap_legend_div = Div(text=heatmap_legend, style={"width": "155px", "margin-bottom": "-10px"})
    # <input type="checkbox" id="expand-toggle" /><label for="expand-toggle" id="expand-btn">Toggle</label>
    return select_sample, select_rank, heatmap_legend_div, mytable1


def create_alpha_diversity_tab():
    imgs_shannon = '<img src="plot_shannon.png"/><img src="plot_shannon_diff.png"/>'
    div_plots_shannon = Div(text=imgs_shannon, css_classes=['bk-width-auto'])
    shannon_column = column(div_plots_shannon, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-width-auto-main'])
    shannon_panel = Panel(child=shannon_column, title='Shannon')
    tabs_plots = Tabs(tabs=[shannon_panel], css_classes=['bk-tabs-margin', 'bk-tabs-margin-lr'])
    return tabs_plots


def create_plots_html(plots_list):
    message_no_spdplot = 'Spider plots of performance require at least 3 profiles.'

    text = '<img src="spider_plot.png" />' if 'spider_plot' in plots_list else message_no_spdplot
    plot1 = Panel(child=Div(text=text), title='Relative performance', width=780)

    text = '<img src="spider_plot_recall_precision.png" />' if 'spider_plot_recall_precision' in plots_list else message_no_spdplot
    plot2 = Panel(child=Div(text=text), title='Absolute performance')

    tabs_plots = Tabs(tabs=[plot1, plot2], width=780, css_classes=['bk-tabs-margin'])
    return tabs_plots


def create_beta_diversity_tab(labels, plots_list):
    rank_to_img = {rank: [''] for rank in c.ALL_RANKS}
    for rank in c.ALL_RANKS:
        for label in labels:
            file = os.path.join("by_tool", label, 'beta_diversity_bc_' + rank)
            if file in plots_list:
                rank_to_img[rank][0] = rank_to_img[rank][0] + '<img src=' + '"' + file + '.png' + '"' + '/>'
    div_plots = Div(text=rank_to_img[c.SPECIES][0], css_classes=['bk-width-auto'])

    source = ColumnDataSource(data=rank_to_img)

    select2_rank_sample_callback = CustomJS(args=dict(source=source), code="""
        div_plots.text = source.data[select2_rank.value][0];
    """)

    select2_rank = Select(title="Taxonomic rank:", value=c.SPECIES, options=c.ALL_RANKS, css_classes=['bk-fit-content'])
    select2_rank.js_on_change('value', select2_rank_sample_callback)
    select2_rank_sample_callback.args["select2_rank"] = select2_rank
    select2_rank_sample_callback.args["div_plots"] = div_plots

    beta_div_column = column(select2_rank, div_plots, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-width-auto-main'])
    return beta_div_column


def create_gs_tab(plots_list, tabs_list):
    # Rarefaction curves panel
    imgs = '<img src="gold_standard/rarefaction_curves.png"/><img src="gold_standard/rarefaction_curves_log_scale.png"/>'
    div_plots_rarefaction = Div(text=imgs, css_classes=['bk-width-auto'])
    div_plots_text = Div(text="<div style='margin-top:18px; margin-bottom:0px;'><font color='navy'><ul style='list-style-type:square;margin-bottom:0;margin-top:0;'><li>OPAL always assumes that the samples are from the same environment.</li><li>Dotted lines are accumulation curves.</li></ul></font></div>", css_classes=['bk-width-auto'])
    gs_column_rarefaction = column(div_plots_text, div_plots_rarefaction, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-width-auto-main'])
    rarefaction_panel = Panel(child=gs_column_rarefaction, title="Rarefaction curves")

    # Proportions panel
    imgs_proportions = ''
    for rank in c.ALL_RANKS:
        if os.path.join('gold_standard', rank) in plots_list:
            fig_name = 'gold_standard/' + rank
            imgs_proportions = imgs_proportions + '<img src="' + fig_name + '.png" onclick="showlegend(this, \'' + rank + '_legend\')" class="proportions"/>'
            imgs_proportions = imgs_proportions + '<img src="' + fig_name + '_legend.png" style="visibility:hidden;display:none;" id="' + rank + '_legend" class="legend">'
    if len(imgs_proportions) > 0:
        imgs_proportions = "<div style='margin-top:18px; margin-bottom:16px;'><font color='navy'><u>Hint:</u> click on a plot for its legend and drag it around the screen as necessary. Click on the same plot again to hide the legend.</font></div>" + imgs_proportions
        div_plots = Div(text=imgs_proportions, css_classes=['bk-width-auto'])
        gs_column_prop = column(div_plots, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-width-auto-main'])
        proportions_panel = Panel(child=gs_column_prop, title="Proportions")

        tabs_plots = Tabs(tabs=[proportions_panel, rarefaction_panel], css_classes=['bk-tabs-margin', 'bk-tabs-margin-lr'])
        tabs_list.append(Panel(child=tabs_plots, title="Gold standard"))
    else:
        tabs_plots = Tabs(tabs=[rarefaction_panel], css_classes=['bk-tabs-margin', 'bk-tabs-margin-lr'])
        tabs_list.append(Panel(child=tabs_plots, title="Gold standard"))


def create_computing_efficiency_tab(pd_metrics, plots_list, tabs_list):
    if 'time_memory' not in plots_list:
        return

    styles = [{'selector': 'td', 'props': [('width', '77pt'), ('padding-left', '10px')]},
              {'selector': 'th', 'props': [('width', 'auto'), ('text-align', 'left'), ('padding-left', '10px')]},
              {'selector': '', 'props': [('width', 'max-content'), ('width', '-moz-max-content'), ('border-spacing', '0px')]}]
    df = pd_metrics[pd_metrics['metric'].isin(['time', 'memory'])]
    df = df.pivot_table(index=['tool'], columns='metric', values='value').reset_index().sort_values(by=['time'])
    df = df[['tool', 'time', 'memory']]
    df.rename(index=str, columns={'tool': 'Tool', 'time': 'Time (hours)', 'memory': 'Memory (GB)'}, inplace=True)

    html = df.style.set_table_styles(styles).hide_index().render()
    div_html = Div(text=html, css_classes=['bk-width-auto', 'bk-inline-block'])

    div_time_memory = Div(text='<img src="time_memory.png"/>', css_classes=['bk-width-auto', 'bk-inline-block'])
    column_time_memory = row(div_html, div_time_memory, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-width-auto-main', 'bk-inline-block'])
    tabs_list.append(Panel(child=column_time_memory, title="Computing efficiency"))


def create_html(pd_rankings, pd_metrics, labels, sample_ids_list, plots_list, output_dir, desc_text):
    col_rankings = create_rankings_html(pd_rankings)

    create_heatmap_bar(output_dir)

    select_sample, select_rank, heatmap_legend_div, mytable1 = create_metrics_table(pd_metrics, labels, sample_ids_list)

    tabs_plots = create_plots_html(plots_list)

    metrics_row = row(column(select_sample, select_rank, heatmap_legend_div, mytable1, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-height-auto', 'bk-inline-block']), column(tabs_plots, sizing_mode='scale_width', css_classes=['bk-width-auto', 'bk-inline-block']), css_classes=['bk-width-auto', 'bk-inline-block'], sizing_mode='scale_width')

    beta_div_column = create_beta_diversity_tab(labels, plots_list)

    tabs_list = [Panel(child=metrics_row, title="Metrics"),
                 Panel(child=col_rankings, title="Rankings"),
                 Panel(child=create_alpha_diversity_tab(), title="Alpha diversity"),
                 Panel(child=beta_div_column, title="Beta diversity")]

    create_computing_efficiency_tab(pd_metrics, plots_list, tabs_list)

    create_gs_tab(plots_list, tabs_list)

    tabs = Tabs(tabs=tabs_list, css_classes=['bk-tabs-margin'])

    title = create_title_div("main", "OPAL: Open-community Profiling Assessment tooL", " produced on {0} with OPAL version {1} ".format(
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), __version__))

    template = Template('''<!DOCTYPE html>
        <html lang="en">
            <head>
                <meta charset="utf-8">
                <title>OPAL: Open-community Profiling Assessment tooL</title>
                {{ js_resources }}
                {{ css_resources }}
                <style>.bk-fit-content {width: fit-content; width: -moz-fit-content;}
                .bk-width-auto {width: auto !important;}
                .bk-height-auto {height: auto !important;}
                .bk-inline-block {display: inline-block !important; float: left;}
                .bk-width-auto-main>div {width: -webkit-fill-available !important;}
                div.bk-width-auto-main {width: -webkit-fill-available !important;}
                .bk-tabs-margin{margin-top: 20px !important;}
                .bk-tabs-margin-lr{margin-left: 10px; margin-right: 10px}
                .bk-root {display: flex; justify-content: center;}
                .bk-padding-top {padding-top: 10px;}
                html {overflow: -moz-scrollbars-vertical; overflow-y: scroll;}
                .tooltip {
                    position: relative;
                    display: inline-block;
                    border-bottom: 1px dashed lightgray;
                    cursor: help;
                }
                .tooltip-sem {
                    border-bottom: 1px dotted black;
                }
                .tooltip .tooltiptext {
                    visibility: hidden;
                    width: 260px;
                    background-color: #555;
                    color: #fff;
                    text-align: center;
                    border-radius: 6px;
                    padding: 5px 0;
                    position: absolute;
                    z-index: 1;
                    bottom: 125%;
                    left: 50%;
                    margin-left: -60px;
                    opacity: 0;
                }
                .tooltip .tooltiptext::after {
                    content: "";
                    position: absolute;
                    top: 100%;
                    left: 50%;
                    margin-left: -5px;
                    border-width: 5px;
                    border-style: solid;
                    border-color: #555 transparent transparent transparent;
                }
                .tooltip:hover .tooltiptext {
                    visibility: visible;
                    opacity: 1;
                }
                .proportions {
                    cursor: pointer;
                }
                .legend {
                    position:absolute;
                    cursor: move;
                    z-index: 1;
                }
                </style>
            </head>
            <body>
                {{ div }}
                {{ script }}
                <script>
                    showlegend = function(img, elementid){
                        var x = document.getElementById(elementid);
                        if (x.style.visibility == 'visible') {
                            x.style.visibility = 'hidden';
                            x.style.display = 'none';
                        } else {
                            if (!x.style.top) {
                                x.style.top = img.offsetTop.toString().concat('px');
                            }
                            if (!x.style.left) {
                                x.style.left = img.offsetLeft.toString().concat('px');
                            }
                            x.style.visibility = 'visible';
                            x.style.display = 'initial';
                        }
                    }
                    function startDrag(e) {
                        if (!e) {
                            var e = window.event;
                        }
                        targ = e.target ? e.target : e.srcElement;
                        if (targ.className != 'legend') {return};
                        offsetX = e.clientX;
                        offsetY = e.clientY;
                        coordX = parseInt(targ.style.left);
                        coordY = parseInt(targ.style.top);
                        drag = true;
                        document.onmousemove=dragDiv;
                        return false;
                    }
                    function dragDiv(e) {
                        if (!drag) {return};
                        if (!e) { var e= window.event};
                        targ.style.left=coordX+e.clientX-offsetX+'px';
                        targ.style.top=coordY+e.clientY-offsetY+'px';
                        return false;
                    }
                    function stopDrag() {
                        drag=false;
                    }
                    window.onload = function() {
                        document.onmousedown = startDrag;
                        document.onmouseup = stopDrag;
                    }
                </script>
            </body>
        </html>
        ''')

    if desc_text:
        data_desc_div = Div(text="""<div style="text-align:left;font-size: 11pt;font-weight: bold;">{}""".format(desc_text), css_classes=['bk-width-auto'])
        html_columns = column(title, data_desc_div, tabs, sizing_mode='scale_width', css_classes=['bk-width-auto-main'])
    else:
        html_columns = column(title, tabs, sizing_mode='scale_width', css_classes=['bk-width-auto-main'])
    script, div = components(html_columns)
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    html = template.render(js_resources=js_resources,
            css_resources=css_resources,
            script=script,
            div=div)

    with open(os.path.join(output_dir, "results.html"), 'w') as f:
        f.write(html)
