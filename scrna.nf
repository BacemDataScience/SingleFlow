#!/usr/bin/env nextflow

params.inputFiles = "Donor 2/Sample 14,Donor 2/Sample 15,Donor 2/Sample 16,Donor 2/Sample 17,Donor 2/Sample 18"
params.colors = "b,g,r,c,m,y,k,w"
params.colorsCustom = "black,yellow,orange,b,c,r,g,w,m"
params.cellCycleCorrect = false
params.customClusters = false
params.deg = false
params.outdir = 'results'
params.genesFile = false
params.rnaVelocity = false
params.phenographClusters = false
params.geneTrends = false
params.geneExpression = false

log.info "scRNA pipeline running"

process loadData {
    output:
    file 'data.pkl' into loadedData
    file 'gene_names.npy' into geneNames
    file 'id_to_name.npy' into idToNameCh

    echo true
    cache 'deep'

    script:
    """
    #!/usr/bin/env pythonw

    import sys
    import numpy as np
    import os

    sys.path.append('$baseDir')

    from load_data import load_data

    input_files = "${params.inputFiles}".split(',')
    input_files = ['/'.join(('${baseDir}'.split('/')[:-1])) + '/' + file for file in input_files]

    df, id_to_name = load_data(input_files)

    if not '${params.exclude}' == 'null':
        df = df.drop(np.load('${baseDir}/${params.exclude}'))

    print('shape after exclude')
    print(df.shape)

    gene_names = np.array([id_to_name[gene_id] for gene_id in df.columns.values])

    df.to_pickle('data.pkl')
    np.save('gene_names.npy', gene_names)
    np.save('id_to_name.npy', id_to_name)

    out_dir = os.path.join('${baseDir}', '${params.outdir}')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    """
}

process glmpca {
    input:
    file df from loadedData

    output:
    file 'glmpca.pkl' into glmPCACh

    when:
    params.glmpca

    script:
    """
    #!/usr/bin/env Rscript
    require(reticulate)
    source_python('${baseDir}/pickle_reader.py')

    df <- t(read_pickle_file('$df'))

    source('${baseDir}/glmpca.R')

    print(df)

    pca <- glmpca(df, 10)

    print(pca)

    write_pickle_df(pca, 'glmpca.pkl')

    """

}

process normalizeData {
    input:
    file df from loadedData

    output:
    file 'norm_df.pkl' into normalizedData, normalizedData2

    cache 'deep'

    when:
    !params.glmpca

    script:
    """
    #!/usr/bin/env pythonw
    import pandas as pd
    import palantir

    import matplotlib
    import matplotlib.pyplot as plt

    df = pd.read_pickle('${df}')

    norm_df = palantir.preprocess.normalize_counts(df)
    norm_df = palantir.preprocess.log_transform(norm_df)

    norm_df.to_pickle('norm_df.pkl')

    """
}

process cellCycleCorrect {
    input:
    file normDf from normalizedData
    file genes from geneNames

    output:
    file 'corrected_df.pkl' into cccData

    when:
    params.cellCycleCorrect

    exec:
    println 'Cell cycle correction'

    script:
    """
    #!/usr/bin/env pythonw
    import pandas as pd
    import numpy as np
    import sys

    sys.path.append('$baseDir')

    from load_data import cell_cycle_correction

    norm_df = pd.read_pickle('${normDf}')

    gene_names = np.load('${genes}')
    gene_ids = norm_df.columns

    norm_df.columns = gene_names
    norm_df = cell_cycle_correction(norm_df, '${baseDir}/cell_cycle/geneset.gmt')
    norm_df.columns = gene_ids

    norm_df.to_pickle('corrected_df.pkl')
    """
}

process preprocessing {
    input:
    file df from cccData.mix(normalizedData)

    output:
    file 'pca_df.pkl' into pcaData, pcaData2, pcaData3

    echo true
    cache 'deep'

    script:
    """
    #!/usr/bin/env pythonw
    import palantir
    import pandas as pd
    import numpy as np
    import os

    import matplotlib
    import matplotlib.pyplot as plt

    norm_df = pd.read_pickle('${df}')

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    pca_projections, _ = palantir.utils.run_pca(norm_df)
    pca_projections.to_pickle('pca_df.pkl')

    """
}

process tsneEmbedding {
    input:
    file pca from pcaData3.mix(glmPCACh)

    output:
    file 'dm_res.npy' into dmData
    file 'ms_df.pkl' into msData
    file 'tsne_df.pkl' into tsneData1, tsneData2, tsneData4, tsneData5, tsneData6, tsneData7, tsneData8, tsneData9, tsneData10, tsneData11, tsneData12

    cache 'deep'
    echo true

    script:
    """
    #!/usr/bin/env python
    import palantir
    import pandas as pd
    import numpy as np
    import os

    import matplotlib
    import matplotlib.pyplot as plt

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    pca_projections = pd.read_pickle('$pca')

    dm_res = palantir.utils.run_diffusion_maps(pca_projections)
    np.save('dm_res.npy', dm_res)

    ms_data = palantir.utils.determine_multiscale_space(dm_res)
    ms_data.to_pickle('ms_df.pkl')

    tsne = palantir.utils.run_tsne(ms_data)
    tsne.to_pickle('tsne_df.pkl')

    palantir.plot.plot_tsne(tsne)
    plt.savefig(os.path.join(out_dir, 'tsne.svg'))

    """

}

process plotSubsets {
    input:
    file tsneDf from tsneData5

    script:
    """
    #!/usr/bin/env python
    import matplotlib
    import matplotlib.pyplot as plt
    import palantir
    import pandas as pd
    import numpy as np
    import os
    import sys
    sys.path.append('$baseDir')
    from plot import highlight_cells_on_tsne

    tsne = pd.read_pickle('${tsneDf}')

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    samples = "${params.inputFiles}".split(',')
    for sample in samples:
        cluster = tsne.index.values[np.flatnonzero(np.core.defchararray.find(tsne.index.values.astype(str),sample)!=-1)]
        highlight_cells_on_tsne(tsne, cluster)
        plt.savefig(os.path.join(out_dir, 'subset' + sample.split('/')[-1] + '.svg'))
    """
}

process magicImputation {
    input:
    file normDf from normalizedData
    file dmDict from dmData
    file genes from geneNames

    output:
    file 'magic_df.pkl' into magicImputed

    cache 'deep'

    when:
    params.imputation == 'MAGIC'

    script:
    """
    #!/usr/bin/env pythonw
    import palantir
    import pandas as pd
    import numpy as np
    import os

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    gene_names = np.load('${genes}')

    norm_df = pd.read_pickle('${normDf}')
    dm_res = np.load('${dmDict}').item()

    imp_df = palantir.utils.run_magic_imputation(norm_df, dm_res)
    imp_df.to_pickle('magic_df.pkl')

    """
}

process saverImputation {
    output:
    file 'saver_df.pkl' into saverImputed

    cache 'deep'

    when:
    params.imputation == 'SAVER'

    script:
    """
    #!/usr/bin/env pythonw
    import palantir
    """

}

process dcaImputation {
    output:
    file 'dca_df.pkl' into dcaImputed

    cache 'deep'

    when:
    params.imputation == 'DCA'

    script:
    """
    #!/usr/bin/env pythonw
    import palantir
    """

}

process imputedData {
    input:
    file impDf from magicImputed.mix(saverImputed, dcaImputed)

    output:
    file 'imp_df.pkl' into imputedData, imputedData2, imputedData3, imputedData4

    cache 'deep'

    script:
    """
    cp $impDf imp_df.pkl
    """

}

process palantir {
    input:
    file msDf from msData
    file impDf from imputedData
    file tsneDf from tsneData1

    output:
    file 'pr_res.pkl' into pr_res_ch, pr_res_ch2, pr_res_ch3, pr_res_ch4, pr_res_ch5, pr_res_ch6

    echo true
    cache 'deep'

    script:
    """
    #!/usr/bin/env pythonw
    import pandas as pd
    import palantir
    import matplotlib
    import matplotlib.pyplot as plt
    import pickle
    import os

    tsne = pd.read_pickle('${tsneDf}')
    ms_data = pd.read_pickle('${msDf}')
    imp_df = pd.read_pickle('${impDf}')

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    #ENSG00000136997 - MYC
    start_cell = imp_df['ENSG00000136997'].idxmax() ## lowest expression of FGFBP2

    pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500, n_jobs=1)

    with open('pr_res.pkl', 'wb') as output_file:
        pickle.dump(pr_res, output_file)

    palantir.plot.plot_palantir_results(pr_res, tsne)
    plt.savefig(os.path.join(out_dir, 'pr_res.svg'))

    palantir.plot.highlight_cells_on_tsne(tsne, [start_cell])
    plt.savefig(os.path.join(out_dir, 'pr_start_cell.svg'))

    palantir.plot.highlight_cells_on_tsne(tsne, pr_res.branch_probs.columns)
    plt.savefig(os.path.join(out_dir, 'pr_terminal_cell.svg'))

    """
}

process custom_clusters {
    input:
    file tsneDf from tsneData2

    output:
    file 'cell_clusters.npy' into customClusterCh, customClusterCh2

    echo true

    cache false

    when:
    params.customClusters

    script:
    """
    #!/usr/bin/env python
    from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5

    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    from PyQt5.QtCore import QThread, pyqtSignal, Qt
    import pandas as pd
    import numpy as np
    import os
    import sys
    sys.path.append('$baseDir')

    from TsneWindow import TsneWindow

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    tsne = pd.read_pickle('${tsneDf}')

    cell_clusters = []

    colors = "${params.colorsCustom}".split(',')
    qapp = QtWidgets.QApplication(sys.argv)
    app = TsneWindow(tsne, cell_clusters, colors, out_dir)
    app.show()
    qapp.exec_()

    exclude = []
    cell_clusters_dict = dict()
    for i in range(len(cell_clusters)):
        cell_clusters_dict[colors[i]] = cell_clusters[i].values
        exclude = np.append(exclude, np.array(cell_clusters[i].values))

    np.save('${baseDir}/exclude_custom.npy', exclude)

    np.save('cell_clusters.npy', cell_clusters_dict)

    """
}

process phenograph {
    input:
    file tsneDf from tsneData4
    file pcaDf from pcaData

    output:
    file 'cell_clusters.npy' into phenoClustersCh, phenoClustersCh2
    file 'phenograph.pkl' into phenoClustersDfCh

    echo true
    cache 'deep'

    when:
    params.phenographClusters

    script:
    """
    #!/usr/bin/env python
    import matplotlib
    import matplotlib.pyplot as plt
    import palantir
    import pandas as pd
    import numpy as np
    import os
    import sys
    sys.path.append('$baseDir')
    from plot import plot_cell_clusters

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    colors = "${params.colors}".split(',')

    pca_projections = pd.read_pickle('${pcaDf}')
    tsne = pd.read_pickle('${tsneDf}')

    clusters = palantir.utils.determine_cell_clusters(pca_projections)
    clusters.to_pickle('phenograph.pkl')

    plot_cell_clusters(tsne, clusters, colors)

    plt.savefig(os.path.join(out_dir, 'cell_clusters.svg'))

    cell_clusters = []
    for cluster in set(clusters):
        cells = clusters.index[clusters == cluster]
        cell_clusters.append(cells)

    cell_clusters_dict = dict()
    for i in range(len(cell_clusters)):
        cell_clusters_dict[colors[i]] = cell_clusters[i].values

    np.save('cell_clusters.npy', cell_clusters_dict)

    """
}

process pseudotimeBoxplot {
    input:
    file phenograph from phenoClustersDfCh
    file prRes from pr_res_ch

    echo true
    cache 'deep'

    when:
    params.phenographClusters

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np
    import os
    import pickle
    import matplotlib
    import matplotlib.pyplot as plt

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    with open('${prRes}', 'rb') as input_file:
        pr_res = pickle.load(input_file)
    clusters = pd.read_pickle('${phenograph}')
    samples = "${params.inputFiles}".split(',')
    colors = "${params.colors}".split(',')


    data_dict = dict()
    colors_dict = dict()
    for i, sample in enumerate(samples):
        data_dict[sample] = pr_res.pseudotime[np.flatnonzero(
                                        np.core.defchararray.find(pr_res.pseudotime.index.values.astype(str),sample)!=-1)].values
        pheno = clusters[np.flatnonzero(
                                        np.core.defchararray.find(pr_res.pseudotime.index.values.astype(str),sample)!=-1)].values
        color_data = []
        for j in range(len(pheno)):
            color_data.append(colors[pheno[j]])
        colors_dict[sample] = color_data

    vals = [data_dict[k] for k in samples]
    colors = [colors_dict[k] for k in samples]
    xs = []
    for i in range(len(samples)):
        xs.append(np.random.normal(i+1, 0.04, len(vals[i])))

    plt.boxplot(vals, labels=samples, showfliers=False)

    ngroup = len(vals)
    for x, val, color in zip(xs, vals, colors):
        plt.scatter(x, val, c=color, alpha=0.4)

    plt.savefig(os.path.join(out_dir, 'pseudotime_boxplot.svg'))

    if not '${params.clustersComposition}' == 'null':
        cell_clusters = np.load('${baseDir}/${params.clustersComposition}').item()
        for cluster in cell_clusters.keys():
            print(cluster)
            composition = dict()
            composition_subsets = dict()
            for cell in cell_clusters[cluster]:
                sample = cell.split('/')[-1]
                if not sample in composition_subsets:
                    composition_subsets[sample] = 1
                else:
                    composition_subsets[sample] = composition_subsets[sample] + 1
                if not clusters[cell] in composition:
                    composition[clusters[cell]] = 1
                else:
                    composition[clusters[cell]] = composition[clusters[cell]] + 1
            print(composition_subsets)
            print(composition)

    """
}

process kMeans {
    input:
    file pcaDf from pcaData2
    file tsneDf from tsneData9

    echo true

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import sklearn
    import sklearn.cluster
    import palantir
    import matplotlib
    import matplotlib.pyplot as plt
    import os

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    pca = pd.read_pickle('${pcaDf}')
    tsne = pd.read_pickle('${tsneDf}')
    mat = pca.values
    #initial = mat[:len(mat):int(len(mat)/4),]
    #initial = mat[:4,]
    km = sklearn.cluster.KMeans(n_clusters=4)#, init=initial)
    km.fit(mat)
    labels = km.labels_
    km_clusters = pd.Series(labels, index=pca.index)
    palantir.plot.plot_cell_clusters(tsne, km_clusters)
    plt.savefig(os.path.join(out_dir, 'k_means_clusters.svg'))

    """
}

process counts {
    input:
    file df from loadedData
    file tsneDf from tsneData8

    echo true

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import matplotlib
    import matplotlib.pyplot as plt
    import os

    df = pd.read_pickle('${df}')
    tsne = pd.read_pickle('${tsneDf}')
    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    counts = df.sum(axis=1)

    ax = plt.subplot(111)

    ax.scatter(tsne['x'], tsne['y'], c=counts.loc[tsne.index], s=3,
                   cmap=matplotlib.cm.Spectral_r)

    plt.savefig(os.path.join(out_dir, 'counts.svg'))

    """
}

process geneExpression {
    input:
    file impDf from imputedData2
    file genes from geneNames
    file tsneDf from tsneData6

    echo true

    cache false

    when:
    params.geneExpression

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np
    import palantir
    import os
    import matplotlib
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('$baseDir')
    from plot import plot_gene_expression

    out_dir = os.path.join('${baseDir}', '${params.outdir}')
    file_name = os.path.join('${baseDir}', '${params.geneExpression}')

    imp_df = pd.read_pickle('${impDf}')
    imp_df.columns = np.load('${genes}')
    tsne = pd.read_pickle('${tsneDf}')

    with open(file_name) as input_file:
        genes = input_file.read().splitlines()
        genes = [gene for gene in genes if gene]

    palantir.plot.plot_gene_expression(imp_df, tsne, genes)
    plt.savefig(os.path.join(out_dir, 'gene_expression.svg'))

    gene_expression_dir = os.path.join(out_dir, 'gene_expression')
    if not os.path.exists(gene_expression_dir):
        os.mkdir(gene_expression_dir)

    for gene in genes:
        if not plot_gene_expression(imp_df, tsne, [gene]) == 0:
            plt.savefig(os.path.join(gene_expression_dir, gene + '.svg'))
    """
}

process factorAnalysis {
    input:
    file impDf from imputedData4
    file genes from geneNames
    file tsneDf from tsneData11
    file idToName from idToNameCh
    file normDf from normalizedData2

    echo true

    cache false

    when:
    params.factorAnalysis

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np
    import palantir
    import os
    import matplotlib
    import matplotlib.pyplot as plt
    import sys
    sys.path.append('$baseDir')
    from plot import plot_gene_expression
    from factor_analysis import factor_analysis

    out_dir = os.path.join('${baseDir}', '${params.outdir}')
    file_name = os.path.join('${baseDir}', '${params.factorAnalysis}')

    id_to_name = np.load('${idToName}').item()
    norm_df = pd.read_pickle('${normDf}')
    imp_df = pd.read_pickle('${impDf}')
    imp_df.columns = np.load('${genes}')
    tsne = pd.read_pickle('${tsneDf}')

    with open(file_name) as input_file:
        genes = input_file.read().splitlines()
        genes = [gene for gene in genes if gene]
        factor_name = file_name.split('.')[0]
    palantir.plot.plot_gene_expression(imp_df, tsne, genes)
    plt.savefig(os.path.join(out_dir, factor_name.split('/')[-1] + '.svg'))

    factor = [factor_name.split('/')[-1]] + genes
    new_file_name = factor_name + '.csv'
    print('Factor: ')
    print(factor)
    print('New file: ', new_file_name)
    with open(new_file_name, 'w') as output_file:
        for item in factor:
            output_file.write('%s\t' % item)
    factor_analysis(norm_df, id_to_name, tsne, new_file_name, 'custom', out_dir)

    """
}

process rnaVelocity {
    input:
    file tsneDf from tsneData7

    echo true

    when:
    params.rnaVelocity

    script:
    """
    #!/usr/bin/env python
    import sys
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import loompy
    import velocyto as vcy
    import os
    import pandas as pd
    import numpy as np

    def slicer_vectorized(a,start,end):
        b = a.view((str,1)).reshape(len(a),-1)[:,start:end]
        return np.fromstring(b.tostring(),dtype=(str,end-start))

    out_dir = os.path.join('${baseDir}', '${params.outdir}')
    vlm = vcy.VelocytoLoom(os.path.join('$baseDir', '${params.rnaVelocity}'))
    tsne = pd.read_pickle('${tsneDf}')
    include = np.array([x[-9:-3].lower() + x[-2:] + ':' + x[:16] for x in tsne.index])
    has = slicer_vectorized(vlm.ca['CellID'].astype('U'),0,25)
    bool_in = [(x in include) for x in has]
    vlm.filter_cells(bool_in)

    print(include)
    print(has)
    vlm.ca['Clusters'] = np.array([int(x[6:8])-14 for x in vlm.ca['CellID']])
    vlm.ca['ClusterName'] = np.array([x[:8] for x in vlm.ca['CellID']])

    print(vlm.ca['Clusters'])
    print(vlm.ca['ClusterName'])

    colors = "${params.colors}".split(',')
    colors = np.array([[ 0.95,  0.6,  0.1], [0.85,  0.3,  0.1], [ 0.8,  0.02,  0.1], [0.81,  0.43,  0.72352941], [0.2,  0.53,  0.71]])
    colors_dict = dict()
    samples = "${params.inputFiles}".split(',')
    for i, sample in enumerate(samples):
        sample = sample.split('/')[-1].split(' ')[0].lower() + str(sample.split('/')[-1].split(' ')[1])
        colors_dict[sample] = colors[i]

    print(colors_dict)

    vlm.set_clusters(vlm.ca["ClusterName"], cluster_colors_dict=colors_dict)

    vlm.ts = np.column_stack([tsne['x'], tsne['y']])

    vlm.plot_fractions()
    plt.savefig(os.path.join(out_dir, 'spliced_unspliced_fractions.svg'))
    plt.close()

    vlm.score_detection_levels(min_cells_express=10)
    vlm.filter_genes(by_detection_levels=True)

    vlm.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)
    vlm.filter_genes(by_cv_vs_mean=True)
    plt.savefig(os.path.join(out_dir, 'cv_s.svg'))
    plt.close()

    vlm._normalize_S(relative_size=vlm.S.sum(0), target_size=vlm.S.sum(0).mean())
    vlm._normalize_U(relative_size=vlm.U.sum(0), target_size=vlm.U.sum(0).mean())

    vlm.perform_PCA()
    plt.plot(np.cumsum(vlm.pca.explained_variance_ratio_)[:100])
    n_comps = np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_))>0.002))[0][0]
    plt.axvline(n_comps, c="k")
    print(n_comps)
    plt.savefig(os.path.join(out_dir, 'pca.svg'))
    plt.close()

    k = 250
    vlm.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=16)

    vlm.fit_gammas(limit_gamma=False, fit_offset=False)

    vlm.predict_U()
    vlm.calculate_velocity()
    vlm.calculate_shift(assumption="constant_velocity") ## or constant unspliced
    vlm.extrapolate_cell_at_t(delta_t=1.)

    vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
                             n_neighbors=200, knn_random=True, sampled_fraction=0.5)
    vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)
    vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)
    plt.figure(None,(17,7))
    plt.subplot(121)
    vlm.plot_grid_arrows(quiver_scale=10,
        scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True},
        min_mass=5.5, angles='xy', scale_units='xy',
        headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1.5,
        plot_random=False, scale_type="relative")
    plt.subplot(122)
    for zz, (k, v) in enumerate(list(colors_dict.items())[::-1]):
        plt.scatter(1, zz/2., c=v, marker="s", lw=0., edgecolor="k", s=220)
        plt.text(1.1, zz/2., k, fontdict={"va":"center", "size":15})
    plt.xlim(0.8,2)
    plt.ylim(-1,8)
    plt.axis("off")
    plt.savefig(os.path.join(out_dir, 'rna_velocity.svg'))
    plt.close()


    """
}

process computeGeneTrends {
    input:
    file prRes from pr_res_ch2
    file impDf from imputedData3

    output:
    file 'gene_trends.pkl' into computedGT, computedGT2

    echo true
    cache 'deep'

    when:
    params.geneTrends

    script:
    """
    #!/usr/bin/env python
    import palantir
    import pickle
    import pandas as pd

    with open('${prRes}', 'rb') as input_file:
        pr_res = pickle.load(input_file)

    imp_df = pd.read_pickle('${impDf}')

    gene_trends = palantir.presults.compute_gene_trends(pr_res, imp_df)

    with open('gene_trends.pkl', 'wb') as output_file:
        pickle.dump(gene_trends, output_file)

    """
}

process chooseTrajectoryTrendsLate {
    input:
    file prRes from pr_res_ch5
    file gtRes from computedGT2

    output:
    file 'late_trends.pkl' into trendsChLate
    file 'late_clusters.pkl' into geneClustersChLate

    echo true
    cache 'deep'

    when:
    params.geneTrendsLate

    script:
    """
    #!/usr/bin/env python
    import pickle
    import palantir
    import pandas as pd
    import numpy as np
    import os

    from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    from PyQt5.QtCore import QThread, pyqtSignal, Qt

    import matplotlib
    import matplotlib.pyplot as plt

    import sys
    sys.path.append('$baseDir')
    from ChooseTrajectoryWin import ChooseTrajectoryWin

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    with open('${prRes}', 'rb') as input_file:
        pr_res = pickle.load(input_file)

    with open('${gtRes}', 'rb') as input_file:
        gene_trends = pickle.load(input_file)

    start_pseudotime = float('${params.geneTrendsLate}')

    trajectories = []

    qapp = QtWidgets.QApplication(sys.argv)
    app = ChooseTrajectoryWin(pr_res.branch_probs.columns, trajectories)
    app.show()
    qapp.exec_()

    trajectory = trajectories[0]
    print(trajectory, 'was chosen as the terminal cell for gene trend analysis')

    trends = gene_trends[trajectory]['trends']
    trends = trends.loc[:,trends.columns > start_pseudotime]

    gene_clusters = palantir.presults.cluster_gene_trends(trends, k=1000)
    print('Clusters computed')

    with open('late_trends.pkl', 'wb') as output_file:
        pickle.dump(trends, output_file)

    with open('late_clusters.pkl', 'wb') as output_file:
        pickle.dump(gene_clusters, output_file)

    palantir.plot.plot_gene_trend_clusters(trends, gene_clusters)

    plt.savefig(os.path.join(out_dir, 'gene_clusters_late.svg'))

    """
}

process geneTrendsClustersLate {
    input:
    file trends from trendsChLate
    file clusters from geneClustersChLate
    file idToName from idToNameCh
    file tsneDf from tsneData12
    file prRes from pr_res_ch6

    echo true

    cache false

    script:
    """
    #!/usr/bin/env python
    import pickle
    import palantir
    import pandas as pd
    import numpy as np
    import os

    from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    from PyQt5.QtCore import QThread, pyqtSignal, Qt

    import matplotlib
    import matplotlib.pyplot as plt

    import sys
    sys.path.append('$baseDir')
    from GeneTrendsWindow import GeneTrendsWindow

    out_dir = os.path.join('${baseDir}', '${params.outdir}')
    goea_dir = os.path.join('${baseDir}', 'goea')

    tsne = pd.read_pickle('${tsneDf}')
    with open('${prRes}', 'rb') as input_file:
        pr_res = pickle.load(input_file)
    cell_clusters = []
    colors = "${params.colors}".split(',')

    with open('${trends}', 'rb') as input_file:
        trends = pickle.load(input_file)

    with open('${clusters}', 'rb') as input_file:
        clusters = pickle.load(input_file)

    clusters[clusters == 3] = 0
    clusters[clusters == 6] = 0
    clusters[clusters == 4] = 3
    clusters[clusters == 5] = 4
    clusters[clusters == 7] = 4

    palantir.plot.plot_gene_trend_clusters(trends, clusters)
    plt.savefig(os.path.join(out_dir, 'gene_clusters_late.svg'))

    id_to_name = np.load('${idToName}').item()

    for i in np.unique(clusters.values):
        gene_ids = np.array(clusters[clusters == i].index.values)
        genes = [id_to_name[id] for id in gene_ids]
        path = os.path.join(out_dir, 'genes_cluster_late' + str(i) + '.txt')
        path_ids = os.path.join(out_dir, 'gene_ids_cluster_late' + str(i) + '.txt')
        np.savetxt(path, genes, fmt='%s')
        np.savetxt(path_ids, gene_ids, fmt='%s')

    qapp = QtWidgets.QApplication(sys.argv)
    app = GeneTrendsWindow(trends, pr_res, clusters, id_to_name, out_dir, goea_dir, tsne, cell_clusters, colors)
    app.show()
    qapp.exec_()

    """
}

process chooseTrajectoryTrends {
    input:
    file prRes from pr_res_ch3
    file gtRes from computedGT

    output:
    file 'trajectory_trends.pkl' into trendsCh
    file 'gene_clusters.pkl' into geneClustersCh

    echo true
    cache 'deep'

    when:
    params.geneTrends

    script:
    """
    #!/usr/bin/env python
    import pickle
    import palantir
    import pandas as pd
    import numpy as np
    import os

    from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    from PyQt5.QtCore import QThread, pyqtSignal, Qt

    import matplotlib
    import matplotlib.pyplot as plt

    import sys
    sys.path.append('$baseDir')
    from ChooseTrajectoryWin import ChooseTrajectoryWin

    out_dir = os.path.join('${baseDir}', '${params.outdir}')

    with open('${prRes}', 'rb') as input_file:
        pr_res = pickle.load(input_file)

    with open('${gtRes}', 'rb') as input_file:
        gene_trends = pickle.load(input_file)

    trajectories = []

    qapp = QtWidgets.QApplication(sys.argv)
    app = ChooseTrajectoryWin(pr_res.branch_probs.columns, trajectories)
    app.show()
    qapp.exec_()

    trajectory = trajectories[0]

    print(trajectory, 'was chosen as the terminal cell for gene trend analysis')

    trends = gene_trends[trajectory]['trends']
    gene_clusters = palantir.presults.cluster_gene_trends(trends, k=1000)

    with open('trajectory_trends.pkl', 'wb') as output_file:
        pickle.dump(trends, output_file)

    with open('gene_clusters.pkl', 'wb') as output_file:
        pickle.dump(gene_clusters, output_file)

    palantir.plot.plot_gene_trend_clusters(trends, gene_clusters)

    plt.savefig(os.path.join(out_dir, 'gene_clusters.svg'))

    """

}

process geneTrendsClusters {
    input:
    file trends from trendsCh
    file clusters from geneClustersCh
    file idToName from idToNameCh
    file tsneDf from tsneData10
    file prRes from pr_res_ch4

    output:
    file 'cell_clusters_pseudo.npy' into pseudotimeClusterCh, pseudotimeClusterCh2
    echo true

    cache false

    script:
    """
    #!/usr/bin/env python
    import pickle
    import palantir
    import pandas as pd
    import numpy as np
    import os

    from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    from PyQt5.QtCore import QThread, pyqtSignal, Qt

    import matplotlib
    import matplotlib.pyplot as plt

    import sys
    sys.path.append('$baseDir')
    from GeneTrendsWindow import GeneTrendsWindow

    out_dir = os.path.join('${baseDir}', '${params.outdir}')
    goea_dir = os.path.join('${baseDir}', 'goea')

    tsne = pd.read_pickle('${tsneDf}')
    with open('${prRes}', 'rb') as input_file:
        pr_res = pickle.load(input_file)
    cell_clusters = []
    colors = "${params.colors}".split(',')

    with open('${trends}', 'rb') as input_file:
        trends = pickle.load(input_file)

    with open('${clusters}', 'rb') as input_file:
        clusters = pickle.load(input_file)

    clusters[clusters == 4] = 2
    clusters[clusters == 6] = 3
    clusters[clusters == 5] = 4

    palantir.plot.plot_gene_trend_clusters(trends, clusters)
    plt.savefig(os.path.join(out_dir, 'gene_clusters.svg'))

    id_to_name = np.load('${idToName}').item()

    for i in np.unique(clusters.values):
        gene_ids = np.array(clusters[clusters == i].index.values)
        genes = [id_to_name[id] for id in gene_ids]
        path = os.path.join(out_dir, 'genes_cluster' + str(i) + '.txt')
        path_ids = os.path.join(out_dir, 'gene_ids_cluster' + str(i) + '.txt')
        np.savetxt(path, genes, fmt='%s')
        np.savetxt(path_ids, gene_ids, fmt='%s')

    qapp = QtWidgets.QApplication(sys.argv)
    app = GeneTrendsWindow(trends, pr_res, clusters, id_to_name, out_dir, goea_dir, tsne, cell_clusters, colors)
    app.show()
    qapp.exec_()

    cell_clusters_dict = dict()
    for i in range(len(cell_clusters)):
        cell_clusters_dict[colors[i]] = cell_clusters[i].values
    np.save('cell_clusters_pseudo.npy', cell_clusters_dict)

    """
}

process deg {
    input:
    file clusters from phenoClustersCh.mix(customClusterCh, pseudotimeClusterCh)
    file df from loadedData
    file idToName from idToNameCh

    output:
    file '*_*_df.pkl' into degDf mode flatten

    echo true

    cache false

    when:
    params.deg

    script:
    """
    #!/usr/bin/env python

    import itertools
    import numpy as np
    import pandas as pd

    df = pd.read_pickle('${df}')
    cell_clusters = np.load('${clusters}').item()
    id_to_name = np.load('${idToName}').item()

    for combination in itertools.combinations(cell_clusters.keys(), 2):
        a = combination[0]
        b = combination[1]
        c1 = cell_clusters[a]
        c2 = cell_clusters[b]
        c1_df = df.loc[c1].T
        c1_df.columns = [a]*c1_df.shape[1]
        c2_df = df.loc[c2].T
        c2_df.columns = [b]*c2_df.shape[1]
        cluster_df = c1_df.join(c2_df)
        cluster_df.to_pickle(a + '_' + b + '_df.pkl')

    """
}

/*
process pklToRdataSCDE {
    input:
    file x from degDf

    output:
    file '*.RData' into rDataSCDE mode flatten

    when:
    params.scde

    script:
    """
    #!/usr/bin/env python
    from rpy2 import robjects
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    import pandas as pd

    file = '$x'

    df = pd.read_pickle(file)
    colnames = df.columns
    df.columns = [colnames[i]+str(i) for i in range(len(colnames))]
    r_data = pandas2ri.py2ri(df)
    robjects.r.assign(file.split('.')[0], r_data)
    robjects.r("save("+file.split('.')[0]+", file='{}')".format(file.split('.')[0]+'.RData'))

    """
}
*/

process scde {
    input:
    //file x from rDataSCDE
    file x from degDf

    echo true
    cache false

    when:
    params.scde

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import os
    from string import Template

    out_dir = os.path.join('${baseDir}', '${params.outdir}')
    out_dir = os.path.join(out_dir, 'scde')
    df_file = '${x}'
    df = pd.read_pickle(df_file)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    df.to_pickle(os.path.join(out_dir, df_file))

    df.to_csv(os.path.join(out_dir, df_file.split('.')[0]+'.csv'))

    d = {'C1' : df_file.split('_')[0], 'C2' : df_file.split('_')[1]}
    filein = open(os.path.join('${baseDir}', 'scde_job'))
    src = Template(filein.read())
    result = src.substitute(d)

    with open(os.path.join(out_dir, 'scde_job_'+ d['C1'] + '_' + d['C2']), 'w') as output_file:
        output_file.write(result)


    """
}

process createDfExp {
    input:
    file clusters from phenoClustersCh2.mix(customClusterCh2, pseudotimeClusterCh2)
    file df from loadedData
    file idToName from idToNameCh

    output:
    file 'dfExpSet.rds' into dfExpSetCh

    echo true
    cache 'deep'

    when:
    params.deconvolution

    script:
    """
    #!/usr/bin/env Rscript

    require(reticulate)
    require(bseqsc)
    require(Biobase)
    source_python('${baseDir}/pickle_reader.py')

    df <- read_pickle_file('$df')
    clusters <- read_npy('$clusters')

    pData <- do.call(rbind,lapply(clusters,data.frame))
    colnames(pData) <- c('CellID')
    pData\$CellType <- rownames(pData)
    rownames(pData) <- pData\$CellID
    pData\$CellType <- unlist(strsplit(pData\$CellType, '.', fixed = TRUE))[c(TRUE, FALSE)]
    pData\$CellSubset <- as.integer(substring(rownames(pData), first=73, last=74))
    pData <- pData[match(rownames(df),rownames(pData)),]

    metadata <- data.frame(labelDescription=c("Patient gender","Case control status",
                                              "Tumor progress on XYZ scale"),
                                              row.names=colnames(pData))

    phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=metadata)

    dfExpSet <- ExpressionSet(as.matrix(t(df)), phenoData=phenoData)

    saveRDS(dfExpSet, file = "dfExpSet.rds")
    """

}

process createSignatureMatrix {
    input:
    file df from dfExpSetCh

    output:
    file 'signature_basis.rds' into signatureMatCh

    echo true

    when:
    params.deconvolution

    script:
    """
    #!/usr/bin/env Rscript

    library(bseqsc)
    library(Biobase)
    require(reticulate)
    library('xbioc')
    source_python('${baseDir}/pickle_reader.py')

    scRnaData <- readRDS(file = "${df}")

    print(scRnaData)

    markerGenes <- list()

    markerGenes[['b']] <- c('ENSG00000204428', 'ENSG00000112658',
    'ENSG00000148358', 'ENSG00000180549', 'ENSG00000172893',
       'ENSG00000156052')

    markerGenes[['g']] <- c('ENSG00000181450', 'ENSG00000225234',
    'ENSG00000171848','ENSG00000042445', 'ENSG00000172671',
       'ENSG00000138185', 'ENSG00000166167', 'ENSG00000176148',
       'ENSG00000100479',
       'ENSG00000168350',
       'ENSG00000136404', 'ENSG00000261177',
       'ENSG00000175106','ENSG00000168517', 'ENSG00000266473',
       'ENSG00000186665', 'ENSG00000149474',
       'ENSG00000101464', 'ENSG00000244005', 'ENSG00000185198',
       'ENSG00000205784', 'ENSG00000131944',
       'ENSG00000126246', 'ENSG00000225975', 'ENSG00000160285')

    markerGenes[['r']] <- c('ENSG00000115109', 'ENSG00000176171',
    'ENSG00000111057', 'ENSG00000110944', 'ENSG00000160285', 'ENSG00000187583')

    markerGenes[['c']] <- c('ENSG00000063127', 'ENSG00000167618', 'ENSG00000181894',
     'ENSG00000156239', 'ENSG00000134905', 'ENSG00000198176', 'ENSG00000135045',
     'ENSG00000185504', 'ENSG00000091157', 'ENSG00000134202')

    require("ggplot2")
    plot_total <- plotCellTotals(scRnaData, 'CellType', 'CellSubset')
    ggsave(file='${baseDir}/${params.outdir}/plot_cell_totals.svg', plot=plot_total, width=10, height=8)

    B <- bseqsc_basis(scRnaData, markerGenes, clusters = 'CellType', samples = 'CellSubset', ct.scale = TRUE)
    B_plot <- plotBasis(B, markerGenes, Colv = NA, Rowv = NA, layout = '_', col = 'Blues')

    saveRDS(B, file = "signature_basis.rds")


    """
}

process deconvolution {
    input:
    file B from signatureMatCh
    file df from loadedData

    echo true

    when:
    params.deconvolution

    script:
    """
    #!/usr/bin/env Rscript
    library(bseqsc)
    library(Biobase)
    require(reticulate)
    source_python('${baseDir}/pickle_reader.py')

    B <- readRDS(file = "${B}")

    # create bulk sample
    df <- read_pickle_file('$df')
    bulkRnaDf <- data.frame(colSums(df))
    bulkRnaDf\$test <- colSums(df)
    bulkRna <- as.matrix(bulkRnaDf)

    print(head(bulkRnaDf))

    print('fit')
    fit <- bseqsc_proportions(bulkRna, B, verbose = TRUE)

    print('fit done')
    print(t(coef(fit)))

    """
}

process degWindow {
    input:
    file idToName from idToNameCh

    echo true
    cache false

    when:
    params.degWin

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    import numpy as np
    import os
    import sys
    sys.path.append('$baseDir')

    from DEGWindow import DEGWindow

    from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5

    from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    from PyQt5.QtCore import QThread, pyqtSignal, Qt
    from matplotlib.figure import Figure

    df_file = os.path.join('$baseDir', '${params.degWin}')
    goea_dir = os.path.join('${baseDir}', 'goea')

    if '${params.adjP}' == 'null':
        adj_p_vals = True
    else:
        adj_p_vals = '${params.adjP}' == 'true'

    cluster1, cluster2 = df_file.split('/')[-1].split('_')[0:2]
    df = pd.read_csv(df_file, sep=' ')
    id_to_name = np.load('${idToName}').item()
    out_dir = os.path.join('${baseDir}', '${params.outdir}')
    out_dir = os.path.join(out_dir, 'deg')

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    qapp = QtWidgets.QApplication(sys.argv)
    app = DEGWindow(df, id_to_name, cluster1, cluster2, out_dir, goea_dir, adj_p_vals=adj_p_vals)
    app.show()
    qapp.exec_()
    """
}
