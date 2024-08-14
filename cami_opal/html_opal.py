import os
import io
from collections import defaultdict
from cami_opal.utils import constants as c
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
from bokeh.models.widgets import TableColumn, Slider, Div, Select
from bokeh.models import DataTable, CustomJS, TabPanel, Tabs
from bokeh.resources import INLINE
from bokeh.plotting import ColumnDataSource
from bokeh.embed import file_html


TITLE = 'OPAL: Open-community Profiling Assessment tooL'
SUM_OF_SCORES = 'Sum of scores'

TEMPLATE = Template('''
{% from macros import embed %}
<!DOCTYPE html>
<html lang="en">
  {% block head %}
  <head>
  {% block inner_head %}
    <meta charset="utf-8">
    <title>{% block title %}{{ title | e if title else "Bokeh Plot" }}{% endblock %}</title>
  {%  block preamble -%}{%- endblock %}
  {%  block resources %}
    <style>
      html, body {
        box-sizing: border-box;
        height: 100%;
        margin: 8px;
        padding: 0;
      }
      .bk-fit-content {width: fit-content; width: -moz-fit-content;}
      .bk-display-block {display: block !important;}
      .bk-float-left {float: left;}
      .bk-tabs-margin{margin-top: 20px !important;}
      .bk-tabs-margin-lr{margin-left: 10px; margin-right: 10px}
      .bk-root {display: flex; justify-content: center;}
      .bk-padding-top {padding-top: 10px;}
      .bk-combo-box > div:first-child {
        width: auto !important;
        padding-right: 40px;
      }
    </style>
  {%   block css_resources -%}
    {{- bokeh_css if bokeh_css }}
  {%-  endblock css_resources %}
  {%   block js_resources -%}
    {{  bokeh_js if bokeh_js }}
  {%-  endblock js_resources %}
  {%  endblock resources %}
  {%  block postamble %}{% endblock %}
  {% endblock inner_head %}
  </head>
  {% endblock head%}
  {% block body %}
  <body>
  {%  block inner_body %}
  {%    block contents %}
  {%      for doc in docs %}
  {{        embed(doc) if doc.elementid }}
  {%-       for root in doc.roots %}
  {%          block root scoped %}
  {{            embed(root) }}
  {%          endblock %}
  {%        endfor %}
  {%      endfor %}
  {%    endblock contents %}
  {{ plot_script | indent(4) }}
  {%  endblock inner_body %}
  </body>
  {% endblock body%}
</html>
''')


TOOLTIPS = """
    <style>
        tr:hover {outline: 1px solid black;}
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
    </style>
    """


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

    fig.savefig(os.path.join(output_dir, 'heatmap_bar.png'), dpi=100, format='png', bbox_inches='tight', pad_inches=-.001, transparent=True)
    pltx.close(fig)


def create_title_div(id, name, info):
    div = Div(text="""<div style="text-align:left;font-size: 20pt;font-weight: bold;">{1}<span style="float: right;font-size: 10pt;font-weight:normal;">{2}</span>""".format(id, name, info), styles={'width': 'auto'})
    div.stylesheets = ["""
    .bk-clearfix {
      width: -webkit-fill-available;
      width: -moz-available;
    }
    """]
    return div


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

    pd_grouped = pd.concat([pd_grouped, pd_mean_over_samples_str]).sort_index()

    # copy unifrac values to every taxonomic rank
    for index, row in pd_grouped.iterrows():
        pd_grouped.loc[index, c.UNIFRAC] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNIFRAC]
        pd_grouped.loc[index, c.UNW_UNIFRAC] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNW_UNIFRAC]
        pd_grouped.loc[index, c.UNIFRAC_CAMI] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNIFRAC_CAMI]
        pd_grouped.loc[index, c.UNW_UNIFRAC_CAMI] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNW_UNIFRAC_CAMI]
    if c.UNIFRAC + c.UNFILTERED_SUF in pd_grouped.columns:
        for index, row in pd_grouped.iterrows():
            pd_grouped.loc[index, c.UNIFRAC + c.UNFILTERED_SUF] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNIFRAC + c.UNFILTERED_SUF]
            pd_grouped.loc[index, c.UNW_UNIFRAC + c.UNFILTERED_SUF] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNW_UNIFRAC + c.UNFILTERED_SUF]
            pd_grouped.loc[index, c.UNIFRAC_CAMI + c.UNFILTERED_SUF] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNIFRAC_CAMI + c.UNFILTERED_SUF]
            pd_grouped.loc[index, c.UNW_UNIFRAC_CAMI + c.UNFILTERED_SUF] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNW_UNIFRAC_CAMI + c.UNFILTERED_SUF]

    for (rank, sample), g in pd_grouped.groupby(['rank', 'sample']):
        rank_to_sample_pd[rank][sample] = g.reset_index().rename(columns={'tool': 'Tool'}).drop(['rank', 'sample'], axis=1).set_index('Tool').T
        # drop all metrics except unifrac for 'rank independent'
        if rank == 'rank independent':
            rank_to_sample_pd[rank][sample] = rank_to_sample_pd[rank][sample].drop(c.ALL_METRICS[4:])
    return rank_to_sample_pd


def get_formatted_pd_rankings(pd_rankings, labels):
    df_list = []
    for metric, g in pd_rankings.groupby(level=0):
        df = g.reset_index().sort_values('position')
        df_list.append(pd.DataFrame({metric: df['tool'].tolist(), 'score' + metric: df['position'].tolist()}))
    df_sum = pd_rankings.groupby(['tool'])['position'].sum().reset_index().sort_values('position')
    df_list.append(
        pd.DataFrame({SUM_OF_SCORES: df_sum['tool'].tolist(), 'score' + SUM_OF_SCORES: df_sum['position'].tolist()}))
    pd_show = pd.concat(df_list, axis=1)

    pd_plot = pd_rankings.reset_index().pivot(index='tool', columns='metric')
    pd_plot[SUM_OF_SCORES] = pd_plot.sum(axis=1)
    pd_plot = pd_plot.reindex(labels)
    return pd_show, pd_plot


def create_rankings_html(pd_rankings, ranks_scored, labels):
    pd_show, pd_plot = get_formatted_pd_rankings(pd_rankings, labels)

    source = ColumnDataSource(data=dict(x=labels,
                                        top=pd_plot[SUM_OF_SCORES],
                                        recall=pd_plot['position', c.RECALL],
                                        precision=pd_plot['position', c.PRECISION],
                                        l1norm=pd_plot['position', c.L1NORM],
                                        unifrac=pd_plot['position', c.UNIFRAC]))

    callback = CustomJS(args=dict(source=source), code="""
        const data = source.data;
        const wrecall = weight_recall.value;
        const wprecision = weight_precision.value;
        const wl1norm = weight_l1norm.value;
        const wunifrac = weight_unifrac.value;
        const topx = data['top'];
        const recall = data['recall'];
        const precision = data['precision'];
        const l1norm = data['l1norm'];
        const unifrac = data['unifrac'];
        
        for (let i = 0; i < topx.length; i++) {
            topx[i] = recall[i] * wrecall + precision[i] * wprecision + l1norm[i] * wl1norm + unifrac[i] * wunifrac;
        }
        
        source.change.emit();
    """)

    weight_recall = Slider(start=0, end=10, value=1, step=.1, title=c.RECALL + " weight")
    callback.args["weight_recall"] = weight_recall
    weight_recall.js_on_change('value', callback)

    weight_precision = Slider(start=0, end=10, value=1, step=.1, title=c.PRECISION + " weight")
    callback.args["weight_precision"] = weight_precision
    weight_precision.js_on_change('value', callback)

    weight_l1norm = Slider(start=0, end=10, value=1, step=.1, title=c.L1NORM + " weight")
    callback.args["weight_l1norm"] = weight_l1norm
    weight_l1norm.js_on_change('value', callback)

    weight_unifrac = Slider(start=0, end=10, value=1, step=.1, title=c.UNIFRAC + " weight")
    callback.args["weight_unifrac"] = weight_unifrac
    weight_unifrac.js_on_change('value', callback)

    p = figure(x_range=pd_plot.index.to_list(), width=1200, height=400, title=SUM_OF_SCORES + " - lower is better")
    p.vbar(x='x', top='top', source=source, width=0.5, bottom=0, color="firebrick")

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
    col_rankings = column([Div(text="<u>Hint 1:</u> click on the columns of scores for sorting.", styles={"width": "600px", "margin-bottom": "0px"}),
                           Div(text="Taxonomic ranks scored: " + ", ".join(
                               ranks_scored) + ". Scoring is only valid if all assessed tools have results for all the same samples and taxonomic ranks. Lower scores are better.",
                               styles={"width": "600px", "margin-bottom": "0px"}),
                           data_table,
                           Div(text="<u>Hint 2:</u> slide the bars to change the weight of the metrics.", styles={"width": "500px", "margin-top": "18px"}),
                           row(weight_recall, weight_precision),
                           row(weight_l1norm, weight_unifrac),
                           p])
    return col_rankings


def get_colors_and_ranges(name, all_values, df_metrics):
    color1 = 'dodgerblue'
    color2 = 'red'
    hue1 = 12
    hue2 = 240

    metrics = [c.PRECISION, c.RECALL, c.F1_SCORE, c.JACCARD]
    metrics = metrics + [metric + c.UNFILTERED_SUF for metric in metrics]
    if name in metrics:
        return color1, color2, hue1, hue2, 0, 1

    metrics = [c.FP, c.UNIFRAC, c.UNW_UNIFRAC, c.UNIFRAC_CAMI, c.UNW_UNIFRAC_CAMI]
    metrics = metrics + [metric + c.UNFILTERED_SUF for metric in metrics]
    if name in metrics:
        return color2, color1, hue2, hue1, 0, max(all_values)

    if name == c.TP or name == c.TP + c.UNFILTERED_SUF:
        return color1, color2, hue1, hue2, 0, max(all_values)

    if name == c.FN or name == c.FN + c.UNFILTERED_SUF:
        fn_values = df_metrics.loc[name, ].values
        # convert "<mean> (<standard error>)" to float of <mean>
        if len(fn_values) > 0 and isinstance(fn_values[0], str):
            fn_values = [float(x.split(' ')[0]) for x in fn_values]
        return color2, color1, hue2, hue1, 0, max(fn_values)
    if name == c.L1NORM or name == c.L1NORM + c.UNFILTERED_SUF:
        return color2, color1, hue2, hue1, 0, 2
    if name == c.BRAY_CURTIS or name == c.BRAY_CURTIS + c.UNFILTERED_SUF:
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

    if pd_series.name == c.SUM_ABUNDANCES or pd_series.name == c.SUM_ABUNDANCES + c.UNFILTERED_SUF:
        return ['background-color: white' for x in values]

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
    estimates_metrics = [c.SUM_ABUNDANCES, c.UNIFRAC_CAMI, c.UNW_UNIFRAC_CAMI, c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.BRAY_CURTIS]
    alpha_diversity_metrics = [c.OTUS, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT]
    rank_independent_metrics = [c.UNIFRAC_CAMI, c.UNW_UNIFRAC_CAMI, c.UNIFRAC, c.UNW_UNIFRAC]

    if c.FP + c.UNFILTERED_SUF in pd_metrics['metric'].values:
        presence_metrics = [[metric, metric + c.UNFILTERED_SUF] for metric in presence_metrics]
        presence_metrics = [metric for elem in presence_metrics for metric in elem]
        estimates_metrics = [[metric, metric + c.UNFILTERED_SUF] for metric in estimates_metrics]
        estimates_metrics = [metric for elem in estimates_metrics for metric in elem]
        alpha_diversity_metrics = [[metric, metric + c.UNFILTERED_SUF] for metric in alpha_diversity_metrics]
        alpha_diversity_metrics = [metric for elem in alpha_diversity_metrics for metric in elem]
        rank_independent_metrics = [[metric, metric + c.UNFILTERED_SUF] for metric in rank_independent_metrics]
        rank_independent_metrics = [metric for elem in rank_independent_metrics for metric in elem]

    all_metrics = [presence_metrics, estimates_metrics, alpha_diversity_metrics]

    presence_metrics_label = 'Presence/absence of taxa'
    estimates_metrics_label = 'Abundance estimates'
    alpha_diversity_metrics = 'Alpha diversity'
    all_metrics_labels = [presence_metrics_label, estimates_metrics_label, alpha_diversity_metrics]

    styles = [{'selector': 'td', 'props': [('width', '115pt')]},
              {'selector': 'th', 'props': [('width', '115pt'), ('text-align', 'left')]},
              {'selector': 'th:nth-child(1)', 'props': [('width', '130pt'), ('font-weight', 'normal')]},
              {'selector': '', 'props': [('width', 'max-content'), ('width', '-moz-max-content'), ('border-top', '1px solid lightgray'), ('border-spacing', '0px')]},
              {'selector': 'expand-toggle:checked ~ * .data', 'props': [('background-color', 'white !important')]}]
    styles_hidden_thead = styles + [{'selector': 'thead', 'props': [('display', 'none')]}]

    def get_html_dict(metrics):
        d_dict = {}
        for tuple in metrics:
            d_dict[tuple[0]] = '<div class="tooltip">{}<span class="tooltiptext">{}: {}</span></div>'.format(tuple[0], tuple[0], tuple[1])
        return d_dict

    metrics_tuples = [(c.RECALL, c.TOOLTIP_RECALL),
                      (c.PRECISION, c.TOOLTIP_PRECISION),
                      (c.F1_SCORE, c.TOOLTIP_F1_SCORE),
                      (c.TP, c.TOOLTIP_TP),
                      (c.FP, c.TOOLTIP_FP),
                      (c.FN, c.TOOLTIP_FN),
                      (c.JACCARD, c.TOOLTIP_JACCARD),
                      (c.UNIFRAC, c.TOOLTIP_UNIFRAC),
                      (c.UNW_UNIFRAC, c.TOOLTIP_UNW_UNIFRAC),
                      (c.UNIFRAC_CAMI, c.TOOLTIP_UNIFRAC_CAMI),
                      (c.UNW_UNIFRAC_CAMI, c.TOOLTIP_UNW_UNIFRAC_CAMI),
                      (c.L1NORM, c.TOOLTIP_L1NORM),
                      (c.BRAY_CURTIS, c.TOOLTIP_BRAY_CURTIS),
                      (c.OTUS, c.TOOLTIP_OTUS),
                      (c.SHANNON_DIVERSITY, c.TOOLTIP_SHANNON_DIVERSITY),
                      (c.SHANNON_EQUIT, c.TOOLTIP_SHANNON_EQUIT),
                      (c.SUM_ABUNDANCES, c.SUM_ABUNDANCES)]

    d = get_html_dict(metrics_tuples)

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
            html = io.StringIO()
            first_metrics = True
            for metrics, metrics_label in zip(all_metrics, all_metrics_labels):
                if rank == 'rank independent':
                    metrics = rank_independent_metrics
                    metrics_label = estimates_metrics_label
                html.write('<p style="margin-bottom: auto"><b>{}</b></p>'.format(metrics_label))
                mydf_metrics = mydf.loc[metrics]

                sorted_columns = [x for x in labels if x in mydf_metrics.columns]
                if c.GS in mydf_metrics.columns:
                    sorted_columns.insert(0, c.GS)

                mydf_metrics = mydf_metrics.loc[:, sorted_columns]

                if first_metrics:
                    this_style = styles
                else:
                    this_style = styles_hidden_thead
                if metrics_label == presence_metrics_label or metrics_label == estimates_metrics_label:
                    html.write(mydf_metrics.style.apply(get_heatmap_colors, df_metrics=mydf_metrics, axis=1).format(precision=3).set_table_styles(this_style).to_html())
                else:
                    html.write(mydf_metrics.style.format(precision=3).set_table_styles(this_style).to_html())
                if rank == 'rank independent':
                    break
                first_metrics = False
            html = pattern.sub(translate, html.getvalue())
            rank_to_sample_to_html[rank].append('{}<div style="margin-bottom:10pt;">{}</div>'.format(TOOLTIPS, html))

    mytable1 = Div(text="""<div>{}</div>""".format(rank_to_sample_to_html[c.ALL_RANKS[0]][0]))

    select_rank = Select(title="Taxonomic rank:", value=c.ALL_RANKS[0], options=c.ALL_RANKS + ['rank independent'])

    select_sample = Select(title="Sample:", value='0', options=list(zip(map(str, range(len(all_sample_ids))), all_sample_ids)))

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
    heatmap_legend_div = Div(text=heatmap_legend, styles={"width": "155px"})
    # <input type="checkbox" id="expand-toggle" /><label for="expand-toggle" id="expand-btn">Toggle</label>
    return select_sample, select_rank, heatmap_legend_div, mytable1


def create_alpha_diversity_tab():
    imgs_shannon = '<img src="plot_shannon.png"/><img src="plot_shannon_diff.png"/>'
    div_plots_shannon = Div(text=imgs_shannon)
    shannon_column = column(div_plots_shannon, sizing_mode='scale_width')
    shannon_panel = TabPanel(child=shannon_column, title='Shannon')
    tabs_plots = Tabs(tabs=[shannon_panel])
    return tabs_plots


def create_plots_html(plots_list):
    message_no_spdplot = 'Spider plots of performance require at least 3 profiles.'

    text = '<img src="spider_plot_relative.png" />' if 'spider_plot_relative' in plots_list else message_no_spdplot
    plot1 = TabPanel(child=Div(text=text), title='Relative performance')

    text = '<img src="spider_plot_absolute.png" />' if 'spider_plot_absolute' in plots_list else message_no_spdplot
    plot2 = TabPanel(child=Div(text=text), title='Absolute performance')

    tabs_plots = Tabs(tabs=[plot1, plot2], width=780)
    return tabs_plots


def create_beta_diversity_tab(labels, plots_list):
    rank_to_img = {rank: [''] for rank in c.ALL_RANKS}
    for rank in c.ALL_RANKS:
        for label in labels:
            file = os.path.join("by_tool", label.replace(' ', '_'), 'beta_diversity_bc_' + rank)
            if file in plots_list:
                rank_to_img[rank][0] = rank_to_img[rank][0] + '<img src=' + '"' + file + '.png' + '"' + '/>'
    div_plots = Div(text=rank_to_img[c.SPECIES][0])

    source = ColumnDataSource(data=rank_to_img)

    select2_rank_sample_callback = CustomJS(args=dict(source=source), code="""
        div_plots.text = source.data[select2_rank.value][0];
    """)

    select2_rank = Select(title="Taxonomic rank:", value=c.SPECIES, options=c.ALL_RANKS)
    select2_rank.js_on_change('value', select2_rank_sample_callback)
    select2_rank_sample_callback.args["select2_rank"] = select2_rank
    select2_rank_sample_callback.args["div_plots"] = div_plots

    beta_div_column = column(select2_rank, div_plots, sizing_mode='scale_width')
    return beta_div_column


def create_gs_tab(plots_list, tabs_list):
    # Rarefaction curves panel
    imgs = '<img src="gold_standard/rarefaction_curves.png"/><img src="gold_standard/rarefaction_curves_log_scale.png"/>'
    div_plots_rarefaction = Div(text=imgs)
    div_plots_text = Div(text="<div style='margin-top:18px; margin-bottom:0px;'><ul style='list-style-type:square;margin-bottom:0;margin-top:0;'><li>OPAL always assumes that the samples are from the same environment.</li><li>Dotted lines are accumulation curves.</li></ul></div>")
    gs_column_rarefaction = column(div_plots_text, div_plots_rarefaction, sizing_mode='scale_width')
    rarefaction_panel = TabPanel(child=gs_column_rarefaction, title="Rarefaction curves")

    # Proportions panel
    imgs_proportions = ''
    for rank in c.ALL_RANKS:
        if os.path.join('gold_standard', rank) in plots_list:
            fig_name = 'gold_standard/' + rank
            imgs_proportions = imgs_proportions + '<img src="' + fig_name + '.png"/>'
            imgs_proportions = imgs_proportions + '<img src="' + fig_name + '_legend.png">'
    if len(imgs_proportions) > 0:
        div_plots = Div(text=imgs_proportions)

        gs_column_prop = column(div_plots, sizing_mode='scale_width')
        proportions_panel = TabPanel(child=gs_column_prop, title="Proportions")

        tabs_plots = Tabs(tabs=[proportions_panel, rarefaction_panel])
        tabs_list.append(TabPanel(child=tabs_plots, title="Gold standard"))
    else:
        tabs_plots = Tabs(tabs=[rarefaction_panel])
        tabs_list.append(TabPanel(child=tabs_plots, title="Gold standard"))


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
    div_html = Div(text=html)

    div_time_memory = Div(text='<img src="time_memory.png"/>')
    column_time_memory = row(div_html, div_time_memory, sizing_mode='scale_width')
    tabs_list.append(TabPanel(child=column_time_memory, title="Computing efficiency"))


def create_html(pd_rankings, ranks_scored, pd_metrics, labels, sample_ids_list, plots_list, output_dir, desc_text):
    col_rankings = create_rankings_html(pd_rankings, ranks_scored, labels)

    create_heatmap_bar(output_dir)

    select_sample, select_rank, heatmap_legend_div, mytable1 = create_metrics_table(pd_metrics, labels, sample_ids_list)

    tabs_plots = create_plots_html(plots_list)

    metrics_row = column(column(row(select_sample, select_rank), heatmap_legend_div, mytable1, sizing_mode='scale_width'),
                         column(tabs_plots, sizing_mode='scale_width'), sizing_mode='scale_width')

    beta_div_column = create_beta_diversity_tab(labels, plots_list)

    tabs_list = [TabPanel(child=metrics_row, title="Metrics"),
                 TabPanel(child=col_rankings, title="Rankings"),
                 TabPanel(child=create_alpha_diversity_tab(), title="Alpha diversity"),
                 TabPanel(child=beta_div_column, title="Beta diversity")]

    create_computing_efficiency_tab(pd_metrics, plots_list, tabs_list)

    create_gs_tab(plots_list, tabs_list)

    tabs = Tabs(tabs=tabs_list)

    title = create_title_div("main", TITLE, " produced on {0} with OPAL version {1} ".format(
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), __version__))

    if desc_text:
        data_desc_div = Div(text="""<div style="text-align:left;font-size: 11pt;font-weight: bold;">{}""".format(desc_text))
        html_columns = column(title, data_desc_div, tabs, sizing_mode='scale_width')
    else:
        html_columns = column(title, tabs, sizing_mode='scale_width')

    html = file_html(models=html_columns, resources=INLINE, title=TITLE, template=TEMPLATE)

    with open(os.path.join(output_dir, "results.html"), 'w') as f:
        f.write(html)
