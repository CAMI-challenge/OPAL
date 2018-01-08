import os
from utils import constants as c
import pandas as pd

from bokeh.embed import file_html
from bokeh.resources import CDN
from bokeh.plotting import figure, output_file, show
from bokeh.models import ColumnDataSource
from bokeh.layouts import widgetbox, layout, column, row
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn, TextInput, Slider, Div, RadioGroup
from bokeh.models import (ColorBar,
                          Text,
                          BasicTicker,
                          HoverTool,
                          FuncTickFormatter,
                          DataTable,
                          widgets,
                          CustomJS)


def create_html(pd_rankings, pd_metrics, labels, output_dir):
    all_metrics = 'Sum of scores'
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
        pd.DataFrame({all_metrics: df_sum['tool'].tolist(), 'score' + all_metrics: df_sum['position'].tolist()}))
    df_list_unsorted_pos.append(
        pd.DataFrame({all_metrics: df_sum_unsorted_pos['tool'].tolist(), 'score' + all_metrics: df_sum_unsorted_pos['position'].tolist()}))

    pd_show = pd.concat(df_list, axis=1)
    pd_show_unsorted_pos = pd.concat(df_list_unsorted_pos, axis=1)

    table_source = ColumnDataSource(pd_show)

    columns = [
        TableColumn(field=all_metrics, title=all_metrics, sortable=False),
        TableColumn(field='score' + all_metrics, title='', width=50),
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

    top = [float(x) for x in pd_show_unsorted_pos['score' + all_metrics]]
    source = ColumnDataSource(data=dict(x=pd_show_unsorted_pos[all_metrics].tolist(),
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
    p = figure(x_range=pd_show_unsorted_pos[all_metrics].tolist(), plot_width=800, plot_height=400, title=all_metrics + " - lower value is better")
    # topx[:] = [0 for x in pd_show['score' + all_metrics]]
    # p.vbar(x=pd_show[all_metrics].tolist(), width=0.5, bottom=0, top="top", color="firebrick", source=source)
    # source = ColumnDataSource(data=dict(top=[]))
    p.vbar(x='x', top='top', source=source, width=0.5, bottom=0, color="firebrick")
    #p.vbar(x=pd_show[all_metrics].tolist(), width=0.5, bottom=0, color="firebrick", top=top)


    l = layout([data_table,
                weight_recall,
                weight_precision,
                weight_l1norm,
                weight_unifrac,
                p])

    # l = row(
    #     widgetbox(weight_recall),
    #     column(p)
    # )

    html = file_html(l, CDN, "OPAL")  # sizing_mode='scale_width'
    file = open(os.path.join(output_dir, "rankings.html"), 'w+')
    file.write(html)
    file.close()




