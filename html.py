import os
from collections import defaultdict
from utils import constants as c
import pandas as pd
import datetime
import numpy as np
import math
import re
from jinja2 import Template

import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as pltx
import matplotlib.colors as pltc
from matplotlib.colors import rgb2hex

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
            sem_value = pd_sem_over_samples.get_value(index, df_columns[j])
            mean_value = '{:>.{precision}g}'.format(pd_mean_over_samples.get_value(index, df_columns[j]), precision=3)
            # if standard error is nan because there is only one sample, then just copy mean
            if np.isnan(sem_value):
                pd_mean_over_samples_str.set_value(index, df_columns[j], "{}".format(mean_value))
            else:
                sem_value = '{:>.{precision}g}'.format(sem_value, precision=3)
                sem_value = '<div class="tooltip tooltip-sem">{}<span class="tooltiptext">standard error: {}</span></div>'.format(sem_value, sem_value)
                pd_mean_over_samples_str.set_value(index, df_columns[j], "{} ({})".format(mean_value, sem_value))
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

    col_rankings = column([data_table,
                           row(weight_recall, weight_precision),
                           row(weight_l1norm, weight_unifrac),
                           p], css_classes=['bk-padding-top'])
    return col_rankings


def get_colors_and_ranges(name, all_values):
    color1 = 'dodgerblue'
    color2 = 'red'
    hue1 = 12
    hue2 = 240
    if name == c.PRECISION or name == c.RECALL or name == c.F1_SCORE or name == c.JACCARD:
        return color1, color2, hue1, hue2, 0, 1
    if name == c.FP or name == c.FN or name == c.UNIFRAC or name == c.UNW_UNIFRAC:
        return color2, color1, hue2, hue1, 0, max(all_values)
    if name == c.TP:
        return color1, color2, hue1, hue2, 0, max(all_values)
    if name == c.L1NORM:
        return color2, color1, hue2, hue1, 0, 2
    if name == c.BRAY_CURTIS:
        return color2, color1, hue2, hue1, 0, 1
    return color1, color2, hue1, hue2, max(all_values), min(all_values)


def get_heatmap_colors(pd_series):
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

    if math.isnan(min(values)):
        red = 'background-color: red'
        return [red for x in values] if not dropped_gs else [''] + [red for x in values]

    color1, color2, hue1, hue2, min_value, max_value = get_colors_and_ranges(pd_series.name, all_values)

    cm = sns.diverging_palette(hue1, hue2, sep=50, l=80, n=15, s=90, as_cmap=True)
    norm = pltc.Normalize(min_value, max_value)

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


def create_metrics_table(pd_metrics, labels):
    rank_to_sample_pd = get_rank_to_sample_pd(pd_metrics)

    set_sample_ids = set()
    for rank in rank_to_sample_pd:
        for sample_id in rank_to_sample_pd[rank]:
            set_sample_ids.add(sample_id)
    all_sample_ids = list(set_sample_ids)
    all_sample_ids.remove('(average over samples)')
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

    styles = [{'selector': 'td', 'props': [('width', '70pt')]},
              {'selector': 'th', 'props': [('width', '70pt'), ('text-align', 'left')]},
              {'selector': 'th:nth-child(1)', 'props': [('width', '120pt'), ('font-weight', 'normal')]},
              {'selector': '', 'props': [('width', 'max-content'), ('width', '-moz-max-content'), ('border-top', '1px solid lightgray'), ('border-spacing', '0px')]},
              {'selector': 'expand-toggle:checked ~ * .data', 'props': [('background-color', 'white !important')]}]
    styles_hidden_thead = styles + [{'selector': 'thead', 'props': [('display', 'none')]}]

    d = {c.RECALL: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.RECALL, c.RECALL),
         c.PRECISION: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.PRECISION, c.PRECISION),
         c.F1_SCORE: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.F1_SCORE, c.F1_SCORE),
         c.TP: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.TP, c.TP),
         c.FP: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.FP, c.FP),
         c.FN: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.FN, c.FN),
         c.JACCARD: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.JACCARD, c.JACCARD),
         c.UNIFRAC: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.UNIFRAC, c.UNIFRAC),
         c.UNW_UNIFRAC: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.UNW_UNIFRAC, c.UNW_UNIFRAC),
         c.L1NORM: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.L1NORM, c.L1NORM),
         c.BRAY_CURTIS: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.BRAY_CURTIS, c.BRAY_CURTIS),
         c.OTUS: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.OTUS, c.OTUS),
         c.SHANNON_DIVERSITY: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.SHANNON_DIVERSITY, c.SHANNON_DIVERSITY),
         c.SHANNON_EQUIT: '<div class="tooltip">{}<span class="tooltiptext">{}</span></div>'.format(c.SHANNON_EQUIT, c.SHANNON_EQUIT)}
    pattern = re.compile('|'.join(map(re.escape, d)))

    def translate(match):
        return d[match.group(0)]

    rank_to_sample_to_html = defaultdict(list)
    for rank in rank_to_sample_pd:
        for sample_id in all_sample_ids:
            if sample_id in rank_to_sample_pd[rank]:
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
                        html += mydf_metrics.style.apply(get_heatmap_colors, axis=1).set_precision(3).set_table_styles(this_style).render()
                    else:
                        html += mydf_metrics.style.set_precision(3).set_table_styles(this_style).render()
                    if rank == 'rank independent':
                        break
                    first_metrics = False
                html = pattern.sub(translate, html)
                rank_to_sample_to_html[rank].append(html)
            else:
                rank_to_sample_to_html[rank].append("")

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

    heatmap_legend = '<img src="heatmap_bar.png" /><div style="text-align:left;font-size: 11px;">Worst<span style="float:right;">Best</span></div>'
    heatmap_legend_div = Div(text=heatmap_legend, style={"width": "155px", "margin-bottom": "-10px"})
    # <input type="checkbox" id="expand-toggle" /><label for="expand-toggle" id="expand-btn">Toggle</label>
    return select_sample, select_rank, heatmap_legend_div, mytable1


def create_plots_html(output_dir):
    message_no_spdplot = 'Spider plots require at least 3 profiles.'

    text = '<img src="spider_plot.png" />' if os.path.exists(os.path.join(output_dir, 'spider_plot.png')) else message_no_spdplot
    plot1 = Panel(child=Div(text=text), title='Relative performance', width=780)

    text = '<img src="spider_plot_recall_precision.png" />' if os.path.exists(os.path.join(output_dir, 'spider_plot_recall_precision.png')) else message_no_spdplot
    plot2 = Panel(child=Div(text=text), title='Absolute performance')

    img_shannon = '<img src="plot_shannon.png" />'
    img_shannon_diff = '<img src="plot_shannon_diff.png" />'
    div_plot_shannon = Div(text=img_shannon)
    source = ColumnDataSource(data=dict(active=[img_shannon, img_shannon_diff]))
    checkbox_callback = CustomJS(args=dict(source=source), code="""
        if(this.active.length > 0) {
            div_plot_shannon.text = source.data['active'][1];
        } else {
            div_plot_shannon.text = source.data['active'][0];
        }
    """)
    checkbox_group = CheckboxGroup(labels=["Absolute differences to gold standard"], active=[], width=300, css_classes=['bk-checkbox-shannon'])
    checkbox_group.js_on_click(checkbox_callback)
    checkbox_callback.args["div_plot_shannon"] = div_plot_shannon

    shannon_column = column(checkbox_group, div_plot_shannon, responsive=True, css_classes=['bk-shannon-plots'])

    plot3 = Panel(child=shannon_column, title='Shannon')
    tabs_plots = Tabs(tabs=[plot1, plot2, plot3], width=780, css_classes=['bk-tabs-margin'])
    return tabs_plots


def create_html(pd_rankings, pd_metrics, labels, output_dir):
    col_rankings = create_rankings_html(pd_rankings)

    create_heatmap_bar(output_dir)

    select_sample, select_rank, heatmap_legend_div, mytable1 = create_metrics_table(pd_metrics, labels)

    tabs_plots = create_plots_html(output_dir)

    metrics_row = row(column(select_sample, select_rank, heatmap_legend_div, mytable1, responsive=True, css_classes=['bk-width-auto']), column(tabs_plots, responsive=True, css_classes=['bk-width-auto']), css_classes=['bk-width-auto'], responsive=True)

    tab1 = Panel(child=metrics_row, title="Metrics")
    tab2 = Panel(child=col_rankings, title="Rankings")

    tabs = Tabs(tabs=[tab1, tab2], css_classes=['bk-tabs-margin'])

    title = create_title_div("main", "OPAL: Profiling Assessment", " produced on {0} with OPAL version {1} ".format(
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), __version__))

    template = Template('''<!DOCTYPE html>
        <html lang="en">
            <head>
                <meta charset="utf-8">
                <title>OPAL: Profiling Assessment</title>
                {{ js_resources }}
                {{ css_resources }}
                <style>.bk-fit-content {width: fit-content; width: -moz-fit-content;}
                .bk-width-auto {width: auto !important;}
                .bk-width-auto-main>div {width: -webkit-fill-available !important;}
                div.bk-width-auto-main {width: -webkit-fill-available !important;}
                .bk-tabs-margin{margin-top: 20px !important;}
                .bk-root {display: flex; justify-content: center;}
                .bk-padding-top {padding-top: 10px;}
                .bk-shannon-plots > div:first-child {
                    padding-bottom: 0px;
                    padding-left: 20px;
                    padding-top: 14px;
                }
                .bk-shannon-plots > div:last-child {
                    padding-top: 0px;
                }
                .bk-checkbox-shannon input[type="checkbox"] {
                    margin: 2px 0 0;
                } 
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
                    width: 120px;
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
                </style>
            </head>
            <body>
                {{ div }}
                {{ script }}
            </body>
        </html>
        ''')

    html_columns = column(title, tabs, responsive=True, css_classes=['bk-width-auto-main'])
    script, div = components(html_columns)
    js_resources = INLINE.render_js()
    css_resources = INLINE.render_css()

    html = template.render(js_resources=js_resources,
            css_resources=css_resources,
            script=script,
            div=div)

    with open(os.path.join(output_dir, "results.html"), 'w') as f:
        f.write(html)
