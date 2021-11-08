import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
import math
import scanpy as sc
import os 
from plotnine import *



def gen_mpl_labels(adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None):
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}

    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)

    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()
        ]
    else:
        texts = [ax.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()]

    adjust_text(texts, **adjust_kwargs)
    

def arrowed_spines(fig, ax, extra_x=0):

    xmin, xmax = ax.get_xlim() 
    ymin, ymax = ax.get_ylim()

    # removing the default axis on all sides:
    for side in ['bottom','right','top','left']:
        ax.spines[side].set_visible(True)

    # removing the axis ticks
    plt.xticks([]) # labels 
    plt.yticks([])
    ax.xaxis.set_ticks_position('none') # tick markers
    ax.yaxis.set_ticks_position('none')

    # get width and height of axes object to compute 
    # matching arrowhead length and width
    dps = fig.dpi_scale_trans.inverted()
    bbox = ax.get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height

    # manual arrowhead width and length
    hw = 1./50.*(ymax-ymin) 
    hl = 1./50.*(xmax-xmin)
    lw = 0.6 # axis line width
    ohg = 0.3 # arrow overhang

    # compute matching arrowhead length and width
    yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width 
    yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height

    # draw x and y axis
    ax.arrow(xmin, ymin, (xmax-xmin)/2, 0., fc='k', ec='k', lw = lw, 
             head_width=hw, head_length=hl, overhang = ohg, 
             length_includes_head= True, clip_on = False) 

    ax.arrow(xmin, ymin, 0., (ymax-ymin)/2, fc='k', ec='k', lw = lw, 
             head_width=yhw, head_length=yhl, overhang = ohg, 
             length_includes_head= True, clip_on = False)
    
    ax.text(xmin+0.1, ymin-0.1, 'UMAP1', va='top', fontsize=8)
    ax.text(xmin-0.8-extra_x, ymin+0.1, 'UMAP2', va='bottom', rotation='vertical', fontsize=8)

    
def plot_umap(adata, groupby, alpha=1, size=10, fontsize=12, legend_loc=None):
    '''Plot umap group by some obs columns'''
    if legend_loc:
        ax = sc.pl.umap(adata, color=groupby, show=False, legend_loc=legend_loc, frameon=False, alpha=alpha, size=size, 
                        legend_fontsize=10, legend_fontweight='normal', edgecolor='none')
        extra_x= 0# 0.2 # right margin 
    else:
        extra_x=0
        ax = sc.pl.umap(adata, color=groupby, show=False, legend_loc=legend_loc, frameon=False, alpha=alpha, size=size, edgecolor='none')
        gen_mpl_labels(
            adata,
            groupby,
            ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
            text_kwargs=dict(fontsize=fontsize, fontweight='normal'))
    fig = ax.get_figure()
    fig.tight_layout()
    arrowed_spines(fig, ax, extra_x)
    fig.savefig("plots/%s.pdf" %groupby, bbox_inches='tight')


    
def plot_volcano(data, genes_volcano, adjust=False, threshold=0.05, fontsize=12, pointsize=8, textsize=10, ticksize=10, figsize=(3.5, 4,5), **kwargs):
    ''' Volcano plot from sig_gene table'''
    #plt.figure(figsize=(7, 10))
    fig, ax = plt.subplots(figsize=figsize)
    xs_up, ys_up, xs_down, ys_down, xs_not, ys_not = [], [], [], [], [], []
    xs_highlight, ys_highlight, texts = [], [], []
    for x, y, l in zip(data['logfoldchanges'], data['pvals_adj'], data['names']):
        if y >= threshold:
            xs_not.append(x)
            ys_not.append(y)
        elif x>0:
            xs_up.append(x)
            ys_up.append(y)
        else:
            xs_down.append(x)
            ys_down.append(y)
        
        if l in genes_volcano:
            texts.append(ax.text(x, -np.log10(y), l, size=textsize))
            xs_highlight.append(x)
            ys_highlight.append(y)
    ax.scatter(xs_not, -np.log10(ys_not), c='grey', edgecolor=(1,1,1,0), label='Not sig', s=pointsize)
    ax.scatter(xs_up, -np.log10(ys_up), c='g', edgecolor=(1,1,1,0), label='FDR<5%, Up', s=pointsize)
    ax.scatter(xs_down, -np.log10(ys_down), c='r', edgecolor=(1,1,1,0), label='FDR<5%, Down', s=pointsize)
    #ax.scatter(xs_highlight, -np.log10(ys_highlight), c='none', edgecolors='k',s=pointsize)
    ax.scatter(xs_highlight, -np.log10(ys_highlight), c='k', s=pointsize)
    
    ax.legend(fontsize=fontsize)
    ax.set_xlabel('$log_2(Fold\:change)$', fontsize=fontsize)
    ax.set_ylabel('$-log_{10}(pval\_adj)$', fontsize=fontsize)
    ylim = int(math.ceil(max(-np.log10(data['pvals_adj'])) * 1.07))
    ax.set_ylim(0, ylim)
    ax.set_axisbelow(True)
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)
    
    if adjust:
        # expand_points=(1, 5) for the original smaller sets of genes
        adjust_text(texts, autoalign='xy', expand_points=(1, 2), arrowprops=dict(arrowstyle="-", color='k', lw=0.6), **kwargs)
        


def dotplot(df, title=''):
    '''Dot plot from given df'''
    description_list = df.sort_values(['gene_ratio', 'description'], ascending=True)['description']
    x = (ggplot(df, aes(x = 'gene_ratio', y = 'description')) 
        + geom_point(aes(size = 'count', color = 'pvals_adj')) 
        + theme_bw(base_size = 12) 
        + scale_colour_gradient(low="red") 
        + ylab('') 
        + ggtitle(title)
        + scale_y_discrete(limits=description_list))

    return(x)


def parse_string_filename(infile):
    directory = os.path.dirname(infile)
    filename = os.path.basename(infile)
    comparison, tag, up_down, ext = filename.split('.')
    up_down = up_down.split('-')[1]
    return(directory, comparison, tag, up_down)

def dotplot_file_to_file(infile, top_n=30):
    ''' parse string output tsv file and make a dotplot to file'''
    directory, comparison, tag, up_down = parse_string_filename(infile)
    outfile = '%s/%s.%s.%s.dotplot.pdf' %(directory, comparison, tag, up_down)
    df = pd.read_csv(infile, delimiter='\t')
    df['gene_ratio'] = df['observed gene count'] / df['background gene count']
    df.rename(columns={'false discovery rate': 'pvals_adj', 
                       'term description': 'description',
                       'observed gene count': 'count'}, inplace=True)
    df = df.sort_values('pvals_adj', ascending=True)
    df_top = df.head(top_n)
    height = min(top_n/30*10, 25) # max height 25
    dotplot(df_top, title='%s-%s' %(comparison, up_down)).save(outfile, format='pdf', width=3.5, height=height, dpi=300)
    

adjust_text_dict = {
    'autoalign': 'xy',
    'expand_points': (2, 2),
    'arrowprops': {
        'arrowstyle': '-',
        'color': 'black',
        'lw': 0.6
    }
}


def revigo_plot(df, df_topn, size=8, title=''):
    ''' RevigoR result to plot'''
    p1 = ggplot(data=df)
    p1 = p1 + geom_point(aes('PlotX', 'PlotY', colour='Value', size='LogSize'), alpha=1, shape='o') + scale_size_area()
    p1 = p1 + scale_color_gradient(low='red', high='yellow', name='$log_{10}(pvalue)$') + scale_size(range=[1,10], name='$log_{10}(GO~size)$') + theme_classic(base_size=size)
    p1 = p1 + geom_text(aes('PlotX', 'PlotY', label='Name'), data=df_topn, size=size, adjust_text=adjust_text_dict)
    p1 = p1 + labs(x = "Semantic space X", y = "Semantic space Y") 
    p1 = p1 + ggtitle(title)
    return p1


def parse_revigo_filename(infile):
    directory = os.path.dirname(infile)
    filename = os.path.basename(infile)
    comparison, tag, up_down, _, _ = filename.split('.')
    return(directory, comparison, tag, up_down)


def revigo_file_to_file(infile, n=10, size=8):
    ''' parse revigo output tabular file and make a revigo plot to file'''
    directory, comparison, tag, up_down = parse_revigo_filename(infile)
    outfile = '%s/%s.%s.%s.revigo.pdf' %(directory, comparison, tag, up_down)
    df = pd.read_csv(infile, skipinitialspace = True, quotechar = '"')
    df = df.sort_values('Dispensability', ascending=False)
    df_topn = df.tail(n) # df sorted descending
    df['Name'] = df['Name'].str.capitalize()
    fig_size = min(n/10*3.5, 20) 
    revigo_plot(df, df_topn, size=size, title='%s-%s' %(comparison, up_down)).save(outfile, format='pdf', width=fig_size, height=fig_size, dpi=300)