# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 17:00:44 2018

@author: Caiyd
"""

import click
import pandas as pd
import seaborn as sns
from bokeh.plotting import figure, ColumnDataSource, output_file, save
from bokeh.models import HoverTool, BoxZoomTool, ResetTool, WheelZoomTool, PanTool, SaveTool, ZoomInTool, ZoomOutTool



def load_dataset(datafile, xcol: str, ycol: str, tags: tuple, groupby):
    """
    tabular形式有header
    tab分割
    """
    cols = list(tags) + [xcol, ycol]
    if groupby:
        cols.append(groupby)
    df = pd.read_csv(datafile, sep='\t', usecols=cols)
    return df



def plot(df, xcol, ycol, tags, groupby, colors, outprefix):
    output_file(f'{outprefix}.html')

    tooltips = [(f"({xcol},{ycol})", f"(@{xcol}, @{ycol})")]
    for tag in tags:
        tooltips.append((f"{tag}", f"@{tag}"))
    hover = HoverTool(tooltips=tooltips)

    p = figure(title="", tools=[hover, BoxZoomTool(), ResetTool(), WheelZoomTool(), PanTool(), SaveTool(), ZoomInTool(), ZoomOutTool()],
               toolbar_location="below", toolbar_sticky=False,
               plot_width=800, plot_height=600)

    if groupby:
        for group, color in zip(df[groupby].unique(), colors):
            source = ColumnDataSource(df.loc[df[groupby] == group, :])
            p.circle(x=xcol, y=ycol, size=10, alpha=0.8,
                     color=color, source=source, legend=group)
        p.legend.location = "top_left"
        p.legend.click_policy="hide"
    else:
        source = ColumnDataSource(df)
        p.circle(x=xcol, y=ycol, size=10, alpha=0.8,
                 color=colors, source=source)


    save(p)




@click.command()
@click.option('--datafile', help='用于绘图的数据文件, tsv格式, 第一行为header')
@click.option('--xcol', help='x轴上展示的变量(header中的一个字段名)')
@click.option('--ycol', help='y轴上展示的变量')
@click.option('--tags', '-t', help='每个数据点需要展示的标签字段名(可多选, e.g -t IID -t GID)', multiple=True)
@click.option('--groupby', help='根据某一字段来分类, 标注不同颜色, default is False', default=False, type=str)
@click.option('--outprefix', help='输出html文件的前缀(输出{outprefix}.html)')
def main(datafile, xcol, ycol, tags, groupby, outprefix):
    df = load_dataset(datafile, xcol, ycol, tags, groupby)
    if groupby:
        ngroup = len(df[groupby].unique())
        colors = sns.color_palette("Set2", ngroup).as_hex()
    else:
        colors = 'blue'
    plot(df, xcol, ycol, tags, groupby, colors, outprefix)


if __name__ == '__main__':
    main()



