import os
from collections import defaultdict
from utils import constants as c
import pandas as pd
import datetime
import numpy as np
import math

import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as pltx
import matplotlib.colors as pltc
from matplotlib.colors import rgb2hex

from version import __version__

from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource
from bokeh.layouts import widgetbox, layout, column, row
from bokeh.models.widgets import TableColumn, Slider, Div, Select, Panel, Tabs
from bokeh.models import (ColorBar,
                          Text,
                          BasicTicker,
                          HoverTool,
                          FuncTickFormatter,
                          DataTable,
                          widgets,
                          CustomJS)


DIV_WIDTH = 1200
DIV_HEIGHT = None


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
    div = Div(text="""<div style="text-align:left;font-size: 20pt;font-weight: bold;">{1}&nbsp;&nbsp;&nbsp;&nbsp;<span style="font-size: 10pt;font-weight:normal;">{2}</span>""".format(id, name, info),
              width=DIV_WIDTH, height=DIV_HEIGHT)
    return [div]


def get_columns(labels, current_columns):
        columns = [c.GS]
        for label in labels:
            if label in current_columns:
                columns.append(label)
        return columns


def get_rank_to_sample_pd(pd_metrics):
    rank_to_sample_pd = defaultdict(dict)

    # copy pd without gold standard
    pd_copy = pd_metrics.copy()
    #pd_copy = pd_copy[pd_copy.tool != c.GS]

    # label ranks, so that the rows don't get lost
    pd_copy.loc[pd_copy[pd.isnull(pd_copy['rank'])].index, 'rank'] = 'rank independent'

    # transform table
    pd_grouped = pd_copy.pivot_table(index=['rank', 'tool', 'sample'], columns='metric', values='value')

    # add sample "(average over samples)"
    pd_mean_over_samples = pd_grouped.groupby(['rank', 'tool'], sort=False).mean().reset_index()
    pd_mean_over_samples['sample'] = '(average over samples)'
    pd_mean_over_samples.set_index(['rank', 'tool', 'sample'], inplace=True)
    pd_grouped = pd.concat([pd_grouped, pd_mean_over_samples])

    # copy unifrac values to every taxonomic rank
    pd_grouped_copy = pd_grouped.copy()
    for index, row in pd_grouped_copy.iterrows():
        pd_grouped.loc[index][c.UNIFRAC] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNIFRAC]
        pd_grouped.loc[index][c.UNW_UNIFRAC] = pd_grouped.loc[('rank independent', index[1], index[2])][c.UNW_UNIFRAC]

    for (rank, sample), g in pd_grouped.groupby(['rank', 'sample']):
        rank_to_sample_pd[rank][sample] = g.reset_index().rename(columns={'tool': 'Tool'}).drop(['rank', 'sample'], axis=1).set_index('Tool').T
        # drop all metrics except unifrac
        if rank == 'rank independent':
            rank_to_sample_pd[rank][sample] = rank_to_sample_pd[rank][sample].drop(c.ALL_METRICS[2:])
    return rank_to_sample_pd


def create_html(pd_rankings, pd_metrics, labels, output_dir):
    sum_of_scores = 'Sum of scores'
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
        pd.DataFrame({sum_of_scores: df_sum['tool'].tolist(), 'score' + sum_of_scores: df_sum['position'].tolist()}))
    df_list_unsorted_pos.append(
        pd.DataFrame({sum_of_scores: df_sum_unsorted_pos['tool'].tolist(), 'score' + sum_of_scores: df_sum_unsorted_pos['position'].tolist()}))

    pd_show = pd.concat(df_list, axis=1)
    pd_show_unsorted_pos = pd.concat(df_list_unsorted_pos, axis=1)

    table_source = ColumnDataSource(pd_show)

    columns = [
        TableColumn(field=sum_of_scores, title=sum_of_scores, sortable=False),
        TableColumn(field='score' + sum_of_scores, title='', width=50),
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

    top = [float(x) for x in pd_show_unsorted_pos['score' + sum_of_scores]]
    source = ColumnDataSource(data=dict(x=pd_show_unsorted_pos[sum_of_scores].tolist(),
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


    #weight_recall = TextInput(title=c.RECALL + " weight", value='1.0', callback=callback)
    weight_recall = Slider(start=0, end=10, value=1, step=.1, title=c.RECALL + " weight", callback=callback)
    callback.args["weight_recall"] = weight_recall

    weight_precision = Slider(start=0, end=10, value=1, step=.1, title=c.PRECISION + " weight", callback=callback)
    callback.args["weight_precision"] = weight_precision

    weight_l1norm = Slider(start=0, end=10, value=1, step=.1, title=c.L1NORM + " weight", callback=callback)
    callback.args["weight_l1norm"] = weight_l1norm

    weight_unifrac = Slider(start=0, end=10, value=1, step=.1, title=c.UNIFRAC + " weight", callback=callback)
    callback.args["weight_unifrac"] = weight_unifrac


    # p.vbar(x=metrics_list, width=0.5, bottom=0, top=pd_show['score'+all_metrics].tolist(), color="firebrick")
    p = figure(x_range=pd_show_unsorted_pos[sum_of_scores].tolist(), plot_width=800, plot_height=400, title=sum_of_scores + " - lower is better")
    # topx[:] = [0 for x in pd_show['score' + all_metrics]]
    # p.vbar(x=pd_show[all_metrics].tolist(), width=0.5, bottom=0, top="top", color="firebrick", source=source)
    # source = ColumnDataSource(data=dict(top=[]))
    p.vbar(x='x', top='top', source=source, width=0.5, bottom=0, color="firebrick")
    #p.vbar(x=pd_show[all_metrics].tolist(), width=0.5, bottom=0, color="firebrick", top=top)


    create_heatmap_bar(output_dir)

    rank_to_sample_pd = get_rank_to_sample_pd(pd_metrics)

    set_sample_ids = set()
    for rank in rank_to_sample_pd:
        for sample_id in rank_to_sample_pd[rank]:
            set_sample_ids.add(sample_id)
    all_sample_ids = list(set_sample_ids)

    presence_metrics = [c.RECALL, c.PRECISION, c.F1_SCORE, c.TP, c.FP, c.FN, c.JACCARD]
    estimates_metrics = [c.UNIFRAC, c.UNW_UNIFRAC, c.L1NORM, c.BRAY_CURTIS]
    alpha_diversity_metics = [c.OTUS, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT]
    rank_independent_metrics = [c.UNIFRAC, c.UNW_UNIFRAC]
    all_metrics = [presence_metrics, estimates_metrics, alpha_diversity_metics]

    presence__metrics_label = 'Presence/absence of taxa'
    estimates_metrics_label = 'Abundance estimates'
    alpha_diversity_metics = 'Alpha diversity'
    all_metrics_labels = [presence__metrics_label, estimates_metrics_label, alpha_diversity_metics]

    styles = [{'selector': 'td', 'props': [('width', '70pt')]},
              {'selector': 'th', 'props': [('width', '70pt'), ('text-align', 'left')]},
              {'selector': 'th:nth-child(1)', 'props': [('width', '120pt'), ('font-weight', 'normal')]},
              {'selector': '', 'props': [('width', 'max-content'), ('width', '-moz-max-content'), ('border-top', '1px solid lightgray'), ('border-spacing', '0px')]},
              {'selector': 'expand-toggle:checked ~ * .data', 'props': [('background-color', 'white !important')]}]
    styles_hidden_thead = styles + [{'selector': 'thead', 'props': [('display', 'none')]}]

    def get_colors_and_ranges(pd_series, max_value_w_gs):
        values = pd_series.tolist()
        color1 = 'dodgerblue'
        color2 = 'red'
        hue1 = 12
        hue2 = 240
        if pd_series.name == c.PRECISION or pd_series.name == c.RECALL or pd_series.name == c.F1_SCORE or pd_series.name == c.JACCARD:
            return color1, color2, hue1, hue2, 0, 1
        if pd_series.name == c.FP or pd_series.name == c.FN:
            return color2, color1, hue2, hue1, 0, max(values)
        if pd_series.name == c.TP:
            return color1, color2, hue1, hue2, 0, max_value_w_gs
        return color1, color2, hue1, hue2, max(values), min(values)

    def get_heatmap_colors(pd_series):
        values = pd_series.tolist()
        max_value_w_gs = max(values)

        dropped_gs = False
        if pd_series.index[0] == c.GS:
            pd_series.drop(c.GS)
            values = values[1:]
            dropped_gs = True

        if math.isnan(min(values)):
            red = 'background-color: red'
            return [red for x in values] if not dropped_gs else [''] + [red for x in values]

        color1, color2, hue1, hue2, min_value, max_value = get_colors_and_ranges(pd_series, max_value_w_gs)

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
                    if metrics_label == presence__metrics_label:
                        html += mydf_metrics.style.apply(get_heatmap_colors, axis=1).set_precision(3).set_table_styles(this_style).render()
                    else:
                        html += mydf_metrics.style.set_precision(3).set_table_styles(this_style).render()
                    if rank == 'rank independent':
                        break
                    first_metrics = False
                rank_to_sample_to_html[rank].append(html)
            else:
                rank_to_sample_to_html[rank].append("")

    mytable1 = Div(text="""<div>{}</div>""".format(rank_to_sample_to_html[c.ALL_RANKS[0]][0]))

    #radio_group = RadioGroup(labels=["Option 1", "Option 2", "Option 3"], active=0)
    select_rank = Select(title="Taxonomic rank:", value=c.ALL_RANKS[0], options=c.ALL_RANKS + ['rank independent'])

    # zip(range(len(min_completeness)), min_completeness)
    select_sample = Select(title="Sample:", value='0', options=list(zip(map(str, range(len(all_sample_ids))), all_sample_ids)))

    source = ColumnDataSource(data=rank_to_sample_to_html)

    #select_rank_sample_callback = CustomJS(args=dict(select_rank=select_rank, select_sample=select_sample, mytable1=mytable1), code="""
    select_rank_sample_callback = CustomJS(args=dict(source=source), code="""
        mytable1.text = source.data[select_rank.value][select_sample.value];
    """)
    #select.js_on_change('active', radio_callback)
    select_rank.js_on_change('value', select_rank_sample_callback)
    select_sample.js_on_change('value', select_rank_sample_callback)
    select_rank_sample_callback.args["select_rank"] = select_rank
    select_rank_sample_callback.args["select_sample"] = select_sample
    select_rank_sample_callback.args["mytable1"] = mytable1

    #mytable2 = Div(children=mytable1, text="my text2")

    col_rankings = column([data_table,
                           weight_recall,
                           weight_precision,
                           weight_l1norm,
                           weight_unifrac,
                           p])

    heatmap_legend = '<img src="heatmap_bar.png" /><div style="text-align:left;font-size: 11px;">Worst<span style="float:right;">Best</span></div>'
    checkbox = Div(text=heatmap_legend, style={"width": "155px", "margin-bottom": "-10px"})
    # <input type="checkbox" id="expand-toggle" /><label for="expand-toggle" id="expand-btn">Toggle</label>

    plot1 = Panel(child=Div(text='<img src="spider_plot.png" />'), title='Spider plot (1)')
    plot2 = Panel(child=Div(text='<img src="spider_plot_recall_precision.png" />'), title='Spider plot (2)')
    plot3 = Panel(child=Div(text='<img src="plot_shannon.png" />'), title='Shannon')
    tabs_plots = Tabs(tabs=[plot1, plot2, plot3], width=400)

    tab1 = Panel(child=column(select_rank, checkbox, mytable1, tabs_plots), title="Metrics")
    tab2 = Panel(child=col_rankings, title="Rankings")

    tabs = Tabs(tabs=[tab1, tab2])

    title = create_title_div("main", "OPAL: Profiling Assessment", " produced on {0} with OPAL version {1} ".format(
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), __version__))

    l = layout(title, select_sample, tabs) #, sizing_mode='scale_width'
    # l = layout([row(select_rank, select_sample),
    #             mytable1,
    #             data_table,
    #             weight_recall,
    #             weight_precision,
    #             weight_l1norm,
    #             weight_unifrac,
    #             p])

    # l = row(
    #     widgetbox(weight_recall),
    #     column(p)
    # )

    html = file_html(l, CDN, "OPAL")  # sizing_mode='scale_width'
    file = open(os.path.join(output_dir, "results.html"), 'w+')
    file.write(html)
    file.close()


