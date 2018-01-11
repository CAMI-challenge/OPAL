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
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn, TextInput, Slider, Div, Select, Panel, Tabs
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


def create_title_div(id, name, info):
    div = Div(text="""<h1 id="{0}">{1}</h1><p>{2}</p>""".format(id, name, info),
              width=DIV_WIDTH, height=DIV_HEIGHT)
    return [div]


def get_rank_to_sample_pd(pd_metrics):
    rank_to_sample_pd = defaultdict(dict)

    # copy pd without gold standard
    pd_copy = pd_metrics.copy()
    pd_copy = pd_copy[pd_copy.tool != c.GS]

    # label ranks, so that the rows don't get lost
    pd_copy.loc[pd_copy[pd.isnull(pd_copy['rank'])].index, 'rank'] = 'rank independent'

    # transform table
    pd_grouped = pd_copy.pivot_table(index=['rank', 'tool', 'sample'], columns='metric', values='value')

    # copy unifrac values to every taxonomic rank
    pd_grouped_copy = pd_grouped.copy()
    for index, row in pd_grouped_copy.iterrows():
        if index[1] == c.GS:
            continue
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


    rank_to_sample_pd = get_rank_to_sample_pd(pd_metrics)

    set_sample_ids = set()
    for rank in rank_to_sample_pd:
        for sample_id in rank_to_sample_pd[rank]:
            set_sample_ids.add(sample_id)
    all_sample_ids = list(set_sample_ids)

    alpha_diversity = [c.L1NORM, c.SHANNON_DIVERSITY, c.SHANNON_EQUIT]
    beta_diversity = [c.UNIFRAC, c.UNW_UNIFRAC, c.JACCARD, c.BRAY_CURTIS]
    binary_metrics = [c.RECALL, c.PRECISION, c.F1_SCORE, c.TP, c.FP, c.FN]
    rank_independent_metrics = [c.UNIFRAC, c.UNW_UNIFRAC]
    all_metrics = [alpha_diversity, beta_diversity, binary_metrics]

    alpha_diversity_label = 'Alpha diversity'
    beta_diversity_label = 'Beta diversity'
    binary_metrics_label = 'Binary metrics'
    all_metrics_labels = [alpha_diversity_label, beta_diversity_label, binary_metrics_label]

    styles = [{'selector': 'td', 'props': [('width', '70pt')]},
              {'selector': 'th', 'props': [('width', '70pt'), ('text-align', 'left')]},
              {'selector': 'th:nth-child(1)', 'props': [('width', '120pt'), ('font-weight', 'normal')]},
              {'selector': '', 'props': [('width', 'max-content'), ('width', '-moz-max-content'), ('border-top', '1px solid lightgray'), ('border-spacing', '0px')]}]
    styles_hidden_thead = styles + [{'selector': 'thead', 'props': [('display', 'none')]}]

    def color_negative_red(pd_series):
        values = pd_series.tolist()

        min_value = min(values)
        if math.isnan(min_value):
            return ['' for x in values]

        range1 = 0
        range2 = 240

        if pd_series.name == c.PRECISION or pd_series.name == c.RECALL or pd_series.name == c.F1_SCORE:
            min_value = 0
            max_value = 1
        else:
            min_value = round(min_value, 3) if not math.isnan(min_value) else 0
            max_value = max(values)
            max_value = round(max_value, 3) if not math.isnan(max_value) else 0
            if pd_series.name == c.FP or pd_series.name == c.FN:
                range1 = 240
                range2 = 0
                min_value = 0

        norm = pltc.Normalize(min_value, max_value)
        #norm = pltc.PowerNorm(gamma=1, vmin=min_value, vmax=max_value)

        normed = norm(values)
        cm = sns.diverging_palette(range1, range2, sep=50, l=60, n=5, as_cmap=True)
        heatmap_colors = [rgb2hex(x) for x in pltx.cm.get_cmap(cm)(normed)]
        #print(heatmap_colors)
        return ['background-color: %s' % color for color in heatmap_colors]

        # color = 'red' if val > 0.5 else 'white'
        # return 'background-color: %s' % color

    # cm = sns.light_palette("green", as_cmap=True)


    rank_to_sample_to_html = defaultdict(list)
    for rank in rank_to_sample_pd:
        for sample_id in all_sample_ids:
            if sample_id in rank_to_sample_pd[rank]:
                mydf = rank_to_sample_pd[rank][sample_id]
                mydf.index.name = None
                # mydf = mydf.applymap(lambda x: '{:.3f}'.format(x) if isinstance(x, float) and not x.is_integer() else x) # float_format='%.3f'
                html = ''
                first_metrics = True
                for metrics, metrics_label in zip(all_metrics, all_metrics_labels):
                    if rank == 'rank independent':
                        metrics = rank_independent_metrics
                        metrics_label = beta_diversity_label
                    html += '<p style="margin-bottom: auto"><b>{}</b></p>'.format(metrics_label)
                    mydf_metrics = mydf.loc[metrics]
                    if first_metrics:
                        html += mydf_metrics.style.apply(color_negative_red, axis=1).set_precision(3).set_table_styles(styles).render()
                    else:
                        html += mydf_metrics.style.apply(color_negative_red, axis=1).set_precision(3).set_table_styles(styles_hidden_thead).render()
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

    # TODO: get rankings out of taxonomic rankings tab

    col_rankings = column([data_table,
                           weight_recall,
                           weight_precision,
                           weight_l1norm,
                           weight_unifrac,
                           p])

    tab1 = Panel(child=column(select_rank, mytable1), title="Metrics")
    tab2 = Panel(child=col_rankings, title="Rankings")

    tabs = Tabs(tabs=[tab1, tab2])

    title = create_title_div("main", "OPAL: Profiling Assessment", " produced on {0} with OPAL version {1} ".format(
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), __version__))

    l = layout(title, select_sample, tabs)
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
    file = open(os.path.join(output_dir, "rankings.html"), 'w+')
    file.write(html)
    file.close()


