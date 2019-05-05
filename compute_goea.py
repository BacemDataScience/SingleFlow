import os
import pandas as pd
import numpy as np

def goea(gene_ids, gene_symbols, trajectory, cluster, goea_dir, out_dir): ## list of genes represented by their ensembl id and gene symbol
    ## load ontologies

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    from goatools.obo_parser import GODag
    obodag = GODag(os.path.join(goea_dir, 'go-basic.obo'))

    ## load associations
    from goatools.associations import read_ncbi_gene2go
    geneid2gos_human = read_ncbi_gene2go(os.path.join(goea_dir, 'gene2go'), taxids=[9606])

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
