#from __future__ import print_function

import csv
import os
import scipy.io
import numpy as np

import PyQt5

import pandas as pd
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

import magic

import palantir

import slalom
import pdb
import pandas as pd
import scipy as SP
from slalom import plotFactors, plotRelevance, plotLoadings, saveFA, dumpFA

from scipy.spatial import ConvexHull

import sys

import warnings

with warnings.catch_warnings():
    # catch experimental ipython widget warning
    warnings.simplefilter('ignore')
    import seaborn as sns
    sns.set(context="paper", style='ticks',
            font_scale=1.5, font='Bitstream Vera Sans')

def load_matrix(matrix_dir):
    #genes_path = 'Donor 1/Sample 1/genes_new.csv'
    genes_path = os.path.join(matrix_dir, "genes.tsv")
    gene_ids = np.array([row[0] for row in csv.reader(open(genes_path), delimiter="\t")])
    gene_names = np.array([row[1] for row in csv.reader(open(genes_path), delimiter="\t")])

    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    barcodes = [row[0]+'-'+matrix_dir for row in csv.reader(open(barcodes_path), delimiter="\t")]

    return magic.io.load_mtx(os.path.join(matrix_dir,'matrix.mtx'), cell_axis='col',
                         gene_names=gene_ids, cell_names=barcodes, sparse=False)

def clean_data(data, cell_min_molecules=1000, genes_min_cells=10):
    ms = data.sum(axis=1)
    cs = data.sum()
    return data.loc[ms.index[ms > cell_min_molecules], cs.index[cs > genes_min_cells]]


def load_txt_2(data,annoFiles, niceTerms=True,annoDBs='MSigDB',dataFile_delimiter=',', verbose=True):
    """Load input file for slalom from txt files.
    Loads an txt files and extracts all the inputs required by slalom
    Args:
        dataFile (str): Strong containing the file name of the text file with the expression levels
        dataFile_delimiter (str): delimiter for reading the data_file. Defaults to ','.
        annoFiles (str, list): Either string containing the file name of the txt file with the gene set annotations or a list containing several anotation files. Each line in
                               in an annotattion file corresponds one gene set; a line starts with the name of the gene set and is followed by the annotated genes.
        annoDBs (str, list)      : database file (MsigDB/REACTOME). If several annotation files are provided this hast to be a list of the same length.
        niceTerms    (bool): Indicates whether to nice terms (omit prefix, capitalize, shorten). Defaults to true.
        dataFile_delimiter (str): Delimiter used in dataFile; defaults to ','.
        verbose (bool): Show progress on loading terms (defaults to True).
    Returns:
        An dictionary containing all the inputs required by slalom.
    """
    if type(annoFiles)==np.dtype('str'):
        annoFiles = [annoFiles]

    if type(annoDBs)==np.dtype('str'):
        annoDBs = [annoDBs]

    if type(niceTerms)==np.dtype('bool'):
        niceTerms = [niceTerms]

    if len(annoFiles)>1:
        if len(niceTerms)==1:
            niceTerms = rep(niceTerms,len(annoFiles))

    if not len(annoDBs)==len(annoFiles):
        raise Exception('annoFiles and annoDBs should have the same length')




    #read data file
    df = data
    if verbose==True:
        print('Data file loaded')
    Ilist = list()
    termsList = list()
    i_file = 0
    for annoFile in annoFiles:
        if not os.path.exists(annoFile):
            raise Exception('annotation file (%s) not found' % annoFile)

        annoDB = annoDBs[i_file].lower()
        if not annoDB in ['msigdb','reactome', 'custom']:
            raise Exception('database (db) needs to be either msigdb, reactome or custom')


        with open(annoFile) as f:
            content = [x.strip('\n') for x in f.readlines()]

        content = [anno.split() for anno in content]

        terms = []
        annotated_genes = []
        for anno in content:
            terms.append(anno[0])
            if annoDB=='msigdb':
                anno_lower = [gene.title() for gene in anno[2:]]
            else:
                anno_lower = [gene.title() for gene in anno[1:]]

            annotated_genes.append(anno_lower)
        I = pd.DataFrame(SP.zeros((df.shape[0], len(terms))), index=[ind.title() for ind in df.index], columns=terms)

        for i_anno in range(len(terms)):
            anno_expressed = list()
            for g in annotated_genes[i_anno]:
                if g in I.index:
                    anno_expressed.append(g)
            I.loc[anno_expressed,terms[i_anno]]=1.
            if verbose==True  and SP.mod(i_anno,50)==0:
                print('%i terms out of %i terms loaded for current annotation file' % (i_anno, len(terms)))

        if niceTerms[i_file]==True:
            if annoDB=='msigdb':
                substring='HALLMARK_'
            elif annoDB=='reactome':
                substring='REACTOME_'
            else:
                substring=' '

            terms = [term[term.find(substring)+len(substring):30] for term in terms]
            terms = [term.capitalize().replace('_',' ') for term in terms]
        Ilist.append(I.values)
        termsList.append(terms)
        i_file+=1
        if verbose==True:
            print('Processed annotation file',annoFile)

    data_out = {}
    data_out['terms'] = SP.hstack(termsList)
    data_out['Y'] = df.values.T
    data_out['I'] = SP.hstack(Ilist)
    data_out['genes'] = list(df.index)
    data_out['lab'] = df.columns
    return data_out

def load_txt(df,annoFiles, niceTerms=True,annoDBs='MSigDB',dataFile_delimiter=',', verbose=True):
    """Load input file for slalom from txt files.
    Loads an txt files and extracts all the inputs required by slalom
    Args:
        dataFile (str): Strong containing the file name of the text file with the expression levels
        dataFile_delimiter (str): delimiter for reading the data_file. Defaults to ','.
        annoFiles (str, list): Either string containing the file name of the txt file with the gene set annotations or a list containing several anotation files. Each line in
                               in an annotattion file corresponds one gene set; a line starts with the name of the gene set and is followed by the annotated genes.
        annoDBs (str, list)      : database file (MsigDB/REACTOME). If several annotation files are provided this hast to be a list of the same length.
        niceTerms    (bool): Indicates whether to nice terms (omit prefix, capitalize, shorten). Defaults to true.
        dataFile_delimiter (str): Delimiter used in dataFile; defaults to ','.
        verbose (bool): Show progress on loading terms (defaults to True).
    Returns:
        An dictionary containing all the inputs required by slalom.
    """
    annoFiles = [annoFiles]

    annoDBs = [annoDBs]

    niceTerms = [niceTerms]

    if len(annoFiles)>1:
        if len(niceTerms)==1:
            niceTerms = rep(niceTerms,len(annoFiles))

    if not len(annoDBs)==len(annoFiles):
        raise Exception('annoFiles and annoDBs should have the same length')

    if verbose==True:
        print('Data file loaded')
    Ilist = list()
    termsList = list()
    i_file = 0
    for annoFile in annoFiles:
        if not os.path.exists(annoFile):
            raise Exception('annotation file (%s) not found' % annoFile)

        annoDB = annoDBs[i_file].lower()
        if not annoDB in ['msigdb','reactome', 'custom']:
            raise Exception('database (db) needs to be either msigdb, reactome or custom')


        with open(annoFile) as f:
            content = [x.strip('\n') for x in f.readlines()]

        content = [anno.split() for anno in content]

        terms = []
        annotated_genes = []
        for anno in content:
            terms.append(anno[0])
            if annoDB=='msigdb':
                anno_lower = [gene.title() for gene in anno[2:]]
            else:
                anno_lower = [gene.title() for gene in anno[1:]]

            annotated_genes.append(anno_lower)
        I = pd.DataFrame(SP.zeros((df.shape[0], len(terms))), index=[ind.title() for ind in df.index], columns=terms)

        for i_anno in range(len(terms)):
            anno_expressed = list()
            for g in annotated_genes[i_anno]:
                if g in I.index:
                    anno_expressed.append(g)
            I.loc[anno_expressed,terms[i_anno]]=1.
            if verbose==True  and SP.mod(i_anno,50)==0:
                print('%i terms out of %i terms loaded for current annotation file' % (i_anno, len(terms)))

        if niceTerms[i_file]==True:
            if annoDB=='msigdb':
                substring='HALLMARK_'
            elif annoDB=='reactome':
                substring='REACTOME_'
            else:
                substring=' '

            terms = [term[term.find(substring)+len(substring):30] for term in terms]
            terms = [term.capitalize().replace('_',' ') for term in terms]
        Ilist.append(I.values)
        termsList.append(terms)
        i_file+=1
        if verbose==True:
            print('Processed annotation file',annoFile)

    data_out = {}
    data_out['terms'] = SP.hstack(termsList)
    data_out['Y'] = df.values.T
    data_out['I'] = SP.hstack(Ilist)
    data_out['genes'] = list(df.index)
    data_out['lab'] = df.columns
    return data_out

def cell_cycle_correction(data, annoFile = 'cell_cycle/geneset.gmt', annoDB   = 'custom'):
    if not os.path.exists(annoFile):
        raise Exception("Annotation file not found")
    data_slalom = load_txt(df=data.T,annoFiles=annoFile,annoDBs=annoDB)
    print ("Loaded {:d} cells, {:d} genes".format(data_slalom['Y'].shape[0],data_slalom['Y'].shape[1]))
    print ("Annotation: {:d} terms".format(len(data_slalom['terms'])))
    #I: indicator matrix that assigns genes to pathways
    I = data_slalom['I'] #if loaded from the hdf file change to I = data['IMSigDB']
    #Y: log expresison values
    Y = data_slalom['Y']
    #terms: ther names of the terms
    terms = data_slalom['terms']

    #gene_ids: the ids of the genes in Y
    gene_ids = data_slalom['genes']

    #initialize FA instance, here using a Gaussian noise model and fitting 3 dense hidden factors
    FA = slalom.initFA(Y, terms,I, gene_ids=gene_ids, noise='gauss', nHidden=3, minGenes=10)

    #model training
    FA.train()

    #print diagnostics
    FA.printDiagnostics()

    corrected_data = FA.regressOut(terms=['M phase', 'Dna replication', 'Chromosome segregation','M phase of mitotic cell cycle'])

    full_matrix = data.copy()
    annotated_genes = np.array(data_slalom['genes'])[np.sum(data_slalom['I'], axis=1) != 0]
    full_matrix[annotated_genes] = corrected_data

    return full_matrix

def get_fig(fig=None, ax=None, figsize=[4, 4]):
    """fills in any missing axis or figure with the currently active one
    :param ax: matplotlib Axis object
    :param fig: matplotlib Figure object
    """
    if not fig:
        fig = plt.figure(figsize=figsize)
    if not ax:
        ax = plt.gca()
    return fig, ax

def plot_cell_clusters(tsne, clusters, ax, outline_clusters=None):
    """Plot cell clusters on the tSNE map
    :param tsne: tSNE map
    :param clusters: Results of the determine_cell_clusters function
    """

    print(clusters)

    # Cluster colors
    n_clusters = len(set(clusters))
#    cluster_colors = pd.Series(sns.color_palette(
#        'hls', n_clusters), index=set(clusters))

    cluster_colors = pd.Series(['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'][:n_clusters], index = set(clusters))
    print(cluster_colors)
    # Set up figure
    n_cols = 6
    n_rows = int(np.ceil(n_clusters / n_cols))
    fig = plt.figure(figsize=[2*n_cols, 2 * (n_rows + 2)])
    gs = plt.GridSpec(n_rows + 2, n_cols,
                      height_ratios=np.append([0.75, 0.75], np.repeat(1, n_rows)))

    # Clusters
#    ax = plt.subplot(gs[0:2, 2:4])

    print('plot clusters')
    ax.scatter(tsne['x'], tsne['y'], s=3,
               color=cluster_colors[clusters[tsne.index]])
    ax.set_axis_off()

    # Branch probabilities
    if outline_clusters is not None:
        for i, cluster in enumerate(set(outline_clusters)):
            #row = int(np.floor(i/n_cols))
            #ax = plt.subplot(gs[row+2, i % n_cols])
            #ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3, color='lightgrey')
            cells = outline_clusters.index[outline_clusters == cluster]
    #        ax.scatter(tsne.loc[cells, 'x'], tsne.loc[cells, 'y'],
    #                   s=3, color=cluster_colors[cluster])

    #        x1, y1 = np.random.normal(loc=5, scale=2, size=(2,15))
    #        x2, y2 = np.random.normal(loc=8, scale=2.5, size=(2,13))

    #        plt.scatter(x1, y1)
    #        plt.scatter(x2, y2)

            encircle(tsne.loc[cells, 'x'], tsne.loc[cells, 'y'], ec='k', fc=cluster_colors[i], alpha=0.2, ax=ax)
    #        encircle(x1, y1, ec="k", fc="gold", alpha=0.2)
    #        encircle(x2, y2, ec="orange", fc="none")

#        ax.set_axis_off()
#        ax.set_title(cluster, fontsize=10)

def encircle(x,y, ax=None, **kw):
    if not ax: ax=plt.gca()
    p = np.c_[x,y]
    hull = ConvexHull(p)
    poly = plt.Polygon(p[hull.vertices,:], **kw)
    ax.add_patch(poly)

def plot_tsne(tsne, cells, fig=None, ax=None):
    """Plot tSNE projections of the data
    :param fig: matplotlib Figure object
    :param ax: matplotlib Axis object
    :param title: Title for the plot
    """
    #fig, ax = get_fig(fig=fig, ax=ax)
    #points = ax.scatter(tsne['x'], tsne['y'], s=5, c=cells, cmap=matplotlib.cm.Spectral_r)
    #ax.set_axis_off()
    #plt.colorbar(points, label='cell', ax=ax, ticks=np.linspace(0,4,5))
    #plt.colorbar(cells)

    cmap=matplotlib.cm.Spectral_r
    norm = matplotlib.colors.BoundaryNorm(np.arange(-0.5,5.5,1), cmap.N)
    plt.scatter(tsne['x'], tsne['y'],c=cells,cmap=cmap,norm=norm,s=100,edgecolor='none')
    plt.colorbar(ticks=np.linspace(0,4,5))
    plt.show()

    return fig, ax

def plot_palantir_results(pr_res, tsne):
    """ Plot Palantir results on tSNE
    """

    # Set up figure
    n_branches = pr_res.branch_probs.shape[1]
    n_cols = 6
    n_rows = int(np.ceil(n_branches / n_cols))
    fig = plt.figure(figsize=[2*n_cols, 2 * (n_rows + 2)])
    gs = plt.GridSpec(n_rows + 2, n_cols,
                      height_ratios=np.append([0.75, 0.75], np.repeat(1, n_rows)))

    # Pseudotime
    ax = plt.subplot(gs[0:2, 1:3])
    ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3,
               cmap=matplotlib.cm.plasma, c=pr_res.pseudotime[tsne.index])
    ax.set_axis_off()
    ax.set_title('Pseudotime')

    # Entropy
    ax = plt.subplot(gs[0:2, 3:5])
    ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3,
               cmap=matplotlib.cm.plasma, c=pr_res.entropy[tsne.index])
    ax.set_axis_off()
    ax.set_title('Differentiation potential')

    # Branch probabilities
    order = [2, 3, 1, 4, 0, 5]
    row = 2
    for i, branch in enumerate(pr_res.branch_probs.columns):
        row = int(np.floor(i/n_cols))
        ax = plt.subplot(gs[row+2, order[i]])
        ax.scatter(tsne.loc[:, 'x'], tsne.loc[:, 'y'], s=3,
                   cmap=matplotlib.cm.plasma, c=pr_res.branch_probs.loc[tsne.index, branch])
        ax.set_axis_off()
        ax.set_title(branch[20:], fontsize=10)


def plot_gene_expression(data, tsne, genes):
    """ Plot gene expression on tSNE maps
    :param genes: Iterable of strings to plot on tSNE
    """

    not_in_dataframe = set(genes).difference(data.columns)
    if not_in_dataframe:
        if len(not_in_dataframe) < len(genes):
            print('The following genes were either not observed in the experiment, '
                  'or the wrong gene symbol was used: {!r}'.format(not_in_dataframe))
        else:
            print('None of the listed genes were observed in the experiment, or the '
                  'wrong symbols were used.')
            return 0

    # remove genes missing from experiment
    genes = set(genes).difference(not_in_dataframe)

    # Plot
    fig = FigureGrid(len(genes), 5)

    for g, ax in zip(genes, fig):
        points = ax.scatter(tsne['x'], tsne['y'], c=data.loc[tsne.index,  g], s=3,
                   cmap=matplotlib.cm.Spectral_r)
#        plt.colorbar(points, label='expression', ax=ax)
        #plt.colorbar(z1_plot,cax=ax1)
        ax.set_axis_off()
        ax.set_title(g)

class FigureGrid:
    """
    Generates a grid of axes for plotting
    axes can be iterated over or selected by number. e.g.:
    >>> # iterate over axes and plot some nonsense
    >>> fig = FigureGrid(4, max_cols=2)
    >>> for i, ax in enumerate(fig):
    >>>     plt.plot(np.arange(10) * i)
    >>> # select axis using indexing
    >>> ax3 = fig[3]
    >>> ax3.set_title("I'm axis 3")
    """

    # Figure Grid is favorable for displaying multiple graphs side by side.

    def __init__(self, n: int, max_cols=3, scale=3):
        """
        :param n: number of axes to generate
        :param max_cols: maximum number of axes in a given row
        """

        self.n = n
        self.nrows = int(np.ceil(n / max_cols))
        self.ncols = int(min((max_cols, n)))
        figsize = self.ncols * scale, self.nrows * scale

        # create figure
        self.gs = plt.GridSpec(nrows=self.nrows, ncols=self.ncols)
        self.figure = plt.figure(figsize=figsize)

        # create axes
        self.axes = {}
        for i in range(n):
            row = int(i // self.ncols)
            col = int(i % self.ncols)
            self.axes[i] = plt.subplot(self.gs[row, col])

    def __getitem__(self, item):
        return self.axes[item]

    def __iter__(self):
        for i in range(self.n):
            yield self[i]

    def save(self):
        self.figure.savefig('output/'+id+'/gene_expression.png')
        plt.close(self.figure)

"""
import multiprocessing
from multiprocessing import Process
import time

def deg_multi(df, n, m):
    manager = multiprocessing.Manager()
    return_dict = manager.dict()

    r_process = Process(target = deg, args=(df, n, m, return_dict,))
    r_process.start()

    while r_process.is_alive():
        print("R running...")    # expect to be able to see this print from main process,
        time.sleep(10)            # while R does work in second process

    r_process.join()

    print('return dict values:')
    print(return_dict.values())

    return return_dict['return']
    """

import pickle

def deg(df, n, m):#, return_dict):
    print('deg utils')
    print(n)
    print(m)
    from rpy2.robjects import pandas2ri
    from rpy2.robjects import r
    import rpy2.robjects as robjects
    robjects.r('''
    f <- function(df, n, m) {
                library(scde)
                sg <- factor(c(rep(c('b'),n),rep(c('g'),m)), levels = c("b", "g"))
                # the group factor should be named accordingly
                names(sg) <- colnames(df)
                o.ifm <- scde.error.models(counts = df, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
                o.prior <- scde.expression.prior(models = o.ifm, counts = df, length.out = 400, show.plot = FALSE)
                groups <- factor(c(rep(c('b'),n),rep(c('g'),m)), levels = c("b", "g"))
                names(groups) <- row.names(o.ifm)
                # run differential expression tests on all genes.
                ediff <- scde.expression.difference(o.ifm, df, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)
                # ediff[order(ediff$Z, decreasing  =  TRUE), ]
                # ediff[abs(ediff$Z) > 1.96, ]
                ediff
        }
    ''')

    r_deg = robjects.globalenv['f']

    pandas2ri.activate()

    """
    import pickle
    with open('degres.pickle', 'rb') as input_file:
        res = pickle.load(input_file)
#        print(res)
        pd_df = pandas2ri.ri2py(res)
        pd_df.columns = res.colnames
        pd_df.index = res.rownames
        return pd_df
        """

    # get p values for the differentailly expressed genes. Check the following link:
    #https://hms-dbmi.github.io/scw/differential-expression.html
    print(n)
    print(m)

    res = r_deg(df, n, m)

    pd_df = pandas2ri.ri2py(res)
    pd_df.columns = res.colnames
    pd_df.index = res.rownames

    print('TEST')

#    return_dict['return'] = pd_df

    return pd_df

def goea(gene_ids, gene_symbols, trajectory, cluster, out_dir): ## list of genes represented by their ensembl id and gene symbol
    ## load ontologies

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    from goatools.obo_parser import GODag
    obodag = GODag("goea/go-basic.obo")

    ## load associations
    from goatools.associations import read_ncbi_gene2go
    geneid2gos_human = read_ncbi_gene2go("goea/gene2go", taxids=[9606])

    ## background gene set
    from goea.genes_NCBI_9606_ProteinCoding import GENEID2NT as GeneID2nt_human

    ## GOEA object
    from goatools.go_enrichment import GOEnrichmentStudy
    goeaobj = GOEnrichmentStudy(
        GeneID2nt_human.keys(), # List of mouse protein-coding genes
        geneid2gos_human, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method

    geneid2symbol = {}

    for gene_symbol in gene_symbols:
        for id in GeneID2nt_human.keys():
            if GeneID2nt_human[id][5] == gene_symbol:
                geneid2symbol[int(id)] = gene_symbol

    #from PyEntrezId import Conversion
    #for (gene_id, gene_symbol) in zip(gene_ids, gene_symbols):
#    id = Conversion('herman.netskar@gmail.com')
#        gene_id = id.convert_ensembl_to_entrez(gene_id) ## get entrez
#        geneid2symbol[int(gene_id)] = gene_symbol

    ## identify which id correspond to the genes in the cluster

    ## Run GOEA
    # 'p_' means "pvalue". 'fdr_bh' is the multipletest method we are currently using.

    """
    import rpy2
    from rpy2.robjects import r, pandas2ri

    from rpy2.robjects import pandas2ri
    import rpy2.robjects as robjects
    robjects.r('''
    f <- function(geneNames) {
                library(clusterProfiler)
                kk <- enrichKEGG(geneNames)
                as.data.frame(kk)
        }
    ''')

    r_enrich = robjects.globalenv['f']
    """

#    print(r_enrich.r_repr())

    gene_names = np.array(list(geneid2symbol.keys()))

    print(gene_names)

    """
    pandas2ri.activate()

    res = r_enrich(gene_names)

    res = r_enrich(gene_names, organism="hsa", pvalueCutoff=0.5, pAdjustMethod="BH", qvalueCutoff=0.1)

    print(res)

    print(pandas2ri.ri2py(res))

    return
    """

    geneids_study = geneid2symbol.keys()

    with open(out_dir + '/' + trajectory[-8:] + 'cluster ' + str(cluster) + 'genes.txt', 'w') as f:
        for gene in geneids_study:
            f.write("%s\n" % gene)

    goea_results_all = goeaobj.run_study(geneids_study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

    ## Write the result to file
    goeaobj.wr_xlsx(out_dir + '/' + trajectory[-8:] + 'cluster ' + str(cluster) + 'goea_symbols.xlsx', goea_results_sig, itemid2name=geneid2symbol)
    goeaobj.wr_xlsx(out_dir + '/' + trajectory[-8:] + 'cluster ' + str(cluster) + 'goea_geneids.xlsx', goea_results_sig)

def factor_analysis(norm_df, id_to_name, tsne, annoFile, annoDB):
    df_factors = None
    factors = []

    if not os.path.exists(annoFile):
        raise Exception("Annotation file not found")

    norm_df_copy = norm_df.copy()
    norm_df_copy.columns = np.array([id_to_name[gene_id] for gene_id in norm_df_copy.columns.values])

    data = load_txt_2(data=norm_df_copy.T,annoFiles=annoFile,annoDBs=annoDB)

    print ("Loaded {:d} cells, {:d} genes".format(data['Y'].shape[0],data['Y'].shape[1]))
    print ("Annotation: {:d} terms".format(len(data['terms'])))

    #I: indicator matrix that assigns genes to pathways
    I = data['I'] #if loaded from the hdf file change to I = data['IMSigDB']
    #Y: log expresison values
    Y = data['Y']
    #terms: ther names of the terms
    terms = data['terms']

    #gene_ids: the ids of the genes in Y
    gene_ids = data['genes']

    #initialize FA instance, here using a Gaussian noise model and fitting 3 dense hidden factors
    FA = slalom.initFA(Y, terms,I, gene_ids=gene_ids, noise='gauss', nHidden=3, minGenes=15)

    #model training
    FA.train()

    #print diagnostics
    FA.printDiagnostics()

    #plot results
    #fig = plotRelevance(FA)

    print(terms[0])

    X = FA.getX(terms=[terms[0]])

    if True: ## is the first factor
        df_factors = pd.DataFrame(X, columns=[terms[0]])
        df_factors.index = tsne.index
    else:
        df_factors[terms[0]] = X
    factors.append(terms[0])

#print(df_factors)

    plot_gene_expression(df_factors, tsne, factors)
    plt.show()
#    plt.savefig(os.path.join(out_dir,'factors.pdf'))


    #get factors; analogous getters are implemented for relevance and weights (see docs)
    #X = FA.getX(terms=['G2m checkpoint','P53 pathway'])
    print(FA.getX())
    print(FA.getX().shape)

    #scatter plot of the top two factors
    #fig = plotFactors(X=X, lab=data['lab'], terms=['G2m checkpoint','P53 pathway'], isCont=False)

    #visualize changes for the G2m checkpoint
    #fig = plotLoadings(FA, 'G2m checkpoint', n_genes=20)
