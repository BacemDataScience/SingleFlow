import pandas as pd
import numpy as np
import magic
import os
import csv

import scprep

import slalom

import scipy as SP

import matplotlib
import matplotlib.pyplot as plt

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

def cell_cycle_correction(data, annoFile, annoDB   = 'custom'):
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

def clean_data(data, cell_min_molecules=1000, genes_min_cells=10):
    ms = data.sum(axis=1)
    cs = data.sum()
    return data.loc[ms.index[ms > cell_min_molecules], cs.index[cs > genes_min_cells]]

def load_matrix(matrix_dir):
    #genes_path = 'Donor 1/Sample 1/genes_new.csv'
    genes_path = os.path.join(matrix_dir, "genes.tsv")
    gene_ids = np.array([row[0] for row in csv.reader(open(genes_path), delimiter="\t")])
    gene_names = np.array([row[1] for row in csv.reader(open(genes_path), delimiter="\t")])

    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    barcodes = [row[0]+'-'+matrix_dir for row in csv.reader(open(barcodes_path), delimiter="\t")]

    data = magic.io.load_mtx(os.path.join(matrix_dir,'matrix.mtx'), cell_axis='col',
                         gene_names=gene_ids, cell_names=barcodes, sparse=False)
#    scprep.plot.plot_library_size(data, cutoff=1500)
#    plt.show()

    return data

def load_data(inputs):
    genes_path = os.path.join(inputs[0], "genes.tsv")
    gene_ids = np.array([row[0] for row in csv.reader(open(genes_path), delimiter="\t")])
    gene_names = np.array([row[1] for row in csv.reader(open(genes_path), delimiter="\t")])
    id_to_name = dict(zip(gene_ids, gene_names))

    data = [load_matrix(folder) for folder in inputs]
    data = pd.concat(data)

    print('Shape of data')
    print(data.shape)

    print('Using default values for cleaning the data')
    data = clean_data(data)

    print('Shape of cleaned data')
    print(data.shape)

    return data, id_to_name
