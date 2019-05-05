import numpy as np
import pandas as pd

def deg(df, cell_clusters, a, b):
    c1 = cell_clusters[a]
    c2 = cell_clusters[b]
    c1_df = df.loc[c1].T
    c1_df.columns = a + '_' + c1_df.columns
    c2_df = df.loc[c2].T
    c2_df.columns = b + '_' + c2_df.columns


    print('DEG analysis of ' + a + ' and ' + b)

    n = c1_df.shape[1]
    m = c2_df.shape[1]

    print('Cluster ', a, ' contains ', n, ' cells out of a total of ', df.shape[0], ' cells.')
    print('Cluster ', b, ' contains ', m, ' cells out of a total of ', df.shape[0], ' cells.')

    cluster_df = c1_df.join(c2_df)

    print('deg utils')
    print(n)
    print(m)
    from rpy2.robjects import pandas2ri
    from rpy2.robjects import r
    import rpy2.robjects as robjects
    robjects.r('''
    f <- function(df, a, b, n, m) {
                library(scde)
                sg <- factor(c(rep(c(a),n),rep(c(b),m)), levels = c(a, b))
                # the group factor should be named accordingly
                names(sg) <- colnames(df)
                o.ifm <- scde.error.models(counts = df, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
                o.prior <- scde.expression.prior(models = o.ifm, counts = df, length.out = 400, show.plot = FALSE)
                groups <- factor(c(rep(c(a),n),rep(c(b),m)), levels = c(a, b))
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

    res = r_deg(cluster_df, a,b, n, m)

    pd_df = pandas2ri.ri2py(res)
    pd_df.columns = res.colnames
    pd_df.index = res.rownames

    return pd_df
