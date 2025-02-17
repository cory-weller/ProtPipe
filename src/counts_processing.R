#!/usr/bin/env Rscript
# R/4
#proteomics analysis for DIA-NN and Spectronaut quantity estimates

#### ARG PARSING ###################################################################################
library(optparse)
pwd = getwd()
optparse_indent = '\n                '
option_list = list( 
    make_option(
        "--pgfile",
        default=NULL,
        help=paste(
            'Input file of Protein Group Intensity (from DIA-NN or Spectronaut)',
            'Required.',
            sep=optparse_indent
        )
    ),
    make_option(
        "--out",
        dest="outdir",
        default='output', 
        help=paste(
            'Directory to direct all output. Directory will be created if does not exist).',
            'Defaults to the current working directory:',
            pwd,
            sep=optparse_indent
        )
    ),
    make_option(
        "--labelgene",
        dest="labelgene",
        default=NULL, 
        help='Gene to always label in output plots'
    ),
    make_option(
        "--base",
        dest="log_base",
        default=10, 
        help='Base for log transformation of intensity data. Default: 10'
    ),
    make_option(
        "--normalize",
        default='shift',
        type='character',
        help=paste(
            'shift: adjust sample intensities to match global median by adding a constant',
            'scale: adjust sample intensities to match global median by multiplicative scaling',
            'none: do not normalize',
            sep=optparse_indent
        )
    ),
    make_option(
        "--exclude",
        default=NULL,
        type='character',
        help=paste(
            'semicolon-separated string of files to exclude from analysis'
        )
    ),
    make_option(
        "--sds",
        dest = 'sds',
        default=3,
        type='numeric',
        help=paste(
            'Filter out samples with protein group counts > N standard deviations from the mean.',
            'Increase to higher values for greater tolerance of variance in protein group counts.',
            'Default: 3',
            sep=optparse_indent
        )
    ),
    make_option(
        "--minintensity",
        dest = 'minintensity',
        default=0,
        type='numeric',
        help='Minimum LINEAR (not log) intensity. Default: 0'
    ),
    make_option(
        "--fdr",
        dest = 'fdr_threshold',
        default=0.05,
        type='numeric',
        help=paste(
            'False Discovery Rate threshold for differential abundance analysis.',
            'Default: 0.05',
            sep=optparse_indent
        )
    ),
    make_option(
        "--foldchange",
        dest = 'foldchange',
        default=5,
        type='numeric',
        help=paste(
            'Minimum LINEAR fold change [NOT log, as log base can be modified] for labeling',
            'protein groups in differential abundance analysis. Default: 5 (equivalent to',
            'log10 fold-change threshold of 0.699)',
            sep=optparse_indent
        )
    ),
    make_option(
        "--imputation",
        action = 'store_true',
        default=FALSE, 
        type='logical',
        help=paste(
            'Applies data imputation. Not yet implimented.',
            sep=optparse_indent
        )
    ),
    make_option(
        "--design",
        default=NULL,  
        help=paste(
            'Comma- or tab-delimited, three-column text file specifying the experimental design.',
            'File should contain headers. Header names do not matter; column order DOES matter.',
            'Columns order: <sample_name> <condition> <control>',
            sep=optparse_indent
        )
    ),
    make_option(
        "--neighbors",
        default=15,
        type='numeric',
        help=paste(
            'N Neighbors to use for UMAP. Default: 15',
            sep=optparse_indent
        )
    ),
    make_option(
        "--dry",
        action = 'store_true',
        default=FALSE, 
        type='logical',
        help=paste(
            'Applies data imputation. Not yet implimented.',
            sep=optparse_indent
        )
    )
)



usage_string <- "Rscript %prog --pgfile [filename] --design [filename] [other options] "
opt <- parse_args(OptionParser(usage = usage_string, option_list))


source('src/functions.R')


# Set to TRUE when running interactively for debugging, to set test opts
if(FALSE) {
    opt <- list()
    opt$pgfile <-  'ANXA11_redux/report.pg_matrix.tsv'
    opt$design <-  'ANXA11_redux/design2.tsv'
    opt$outdir <-  'ANXA11_redux/'
    opt$sds <-  3
    opt$normalize <-  'shift'
    opt$log_base <-  2
    #opt$exclude <-  'ANXA11_EMV_1;ANXA11_EMV_2;ANXA11_EMV_3;ANXA11_EMV_4'
    opt$fdr_threshold <- 0.05
    opt$foldchange <- 10
    opt$minintensity <- 0
    opt$neighbors <- 15
    opt$sds <- 3
}

badargs <- FALSE

if(opt$dry) {
    cat("INFO: Quitting due to --dry run\n")
    quit(status=0)
}

if(! opt$normalize %in% c('shift','scale','none')) {
    cat("ERROR: --normalize must be 'shift' 'scale' or 'none'\n")
    badargs <- TRUE
}

if (is.null(opt$pgfile)) {
    cat("ERROR: --pgfile <file> must be provided\n")
    badargs <- TRUE
}

if (is.null(opt$design)) {
    cat("ERROR: --design <file> must be provided\n")
    badargs <- TRUE
}

if (badargs == TRUE) {
    quit(status=1)
}


#### PACKAGES ######################################################################################
package_list = c('ggplot2', 'data.table', 'corrplot', 'umap', 'magick', 'ggdendro', 'ecodist',
                 'ggbeeswarm', 'ggrepel', 'ggthemes', 'foreach')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
    cat("INFO: All packages successfully loaded\n")
} else {
    cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
options(warn = defaultW)    # Turn warnings back on



opt$lfc_threshold <- log(opt$foldchange, base=opt$log_base)
cat(paste0('INFO: LFC threshold of log[', opt$log_base, '](Intensity) > ', opt$lfc_threshold, '\n'))
cat(paste0('INFO: FDR threshold of ', opt$fdr_threshold, '\n'))


#### MAKE DIRS #####################################################################################

QC_dir <- paste0(opt$outdir, '/QC/')
if(! dir.exists(QC_dir)){
    dir.create(QC_dir, recursive = T)
}


cluster_dir <- paste0(opt$outdir, '/clustering/')
if(! dir.exists(cluster_dir)){
    dir.create(cluster_dir, recursive = T)
}


DI_dir <- paste0(opt$outdir, '/differential_intensity/')
if(! dir.exists(DI_dir)){
    dir.create(DI_dir, recursive = T)
}


#### IMPORT AND FORMAT DATA#########################################################################

tryTo(paste0('INFO: Reading input file ', opt$pgfile),{
    dat <- fread(opt$pgfile)
}, paste0('ERROR: problem trying to load ', opt$pgfile, ', does it exist?'))

tryTo(paste0('INFO: Massaging data from ', opt$pgfile, ' into a common style format for processing'), {
    dat <- standardize_format(dat)
}, 'ERROR: failed! Check for missing/corrupt headers?')

tryTo(paste0('INFO: Trimming extraneous column name info'), {
    setnames(dat, trim_colnames(dat))
}, 'ERROR: failed! Check for missing/corrupt headers?')


# exclude samples in opt$exclude
if (! is.null(opt$exclude)) {
    opt$exclude <- strsplit(opt$exclude, split=';')[[1]]

    tryTo(paste('INFO: excluding samples', opt$exclude),{
    for(i in opt$exclude) {
        dat[, (i) := NULL]
    }
}, paste0('ERROR: problem trying to load ', opt$pgfile, ', does it exist?'))
}

tryTo(paste0('INFO: Converting to long format'), {
    dat.long <- melt_intensity_table(dat)
}, 'ERROR: failed! Check for missing/corrupt headers?')


tryTo(paste0('INFO: Applying Log[base', opt$log_base, '](value+1) transformation to intensities'),{
    dat.long <- dat.long[, Intensity := log((Intensity + 1), base=opt$log_base)]
},'ERROR: failed! Was your log base numeric and > 1?')


tryTo('INFO: Excluding all unquantified or zero intensities', {
    dat.long <- dat.long[! is.na(Intensity)][Intensity != 0]
}, 'ERROR: failed!')



#### QC ############################################################################################

## Filtering

tryTo(paste0('INFO: Sorting samples by median intensity'),{
    increasing_levels <- as.character(dat.long[, list('median'=median(Intensity)), by=Sample][order(median)]$Sample)
    dat.long[, Sample := factor(Sample, levels=increasing_levels)]
}, 'ERROR: failed!')


tryTo(paste0('INFO: Applying Filter Log[', opt$log_base, '](Intensity) > ',opt$minintensity),{
    dat.long <- dat.long[Intensity > opt$minintensity]
}, 'ERROR: failed!')


tryTo('INFO: Plotting intensity distribution',{
    plot_flip_beeswarm(dat.long, QC_dir, 'intensities.png', plot_title='Un-normalized intensities')
    # plot_density(dat.long, QC_dir, 'intensity_density.png')
    # plot_density(dat.long.normalized, QC_dir, 'intensity_density_normalized.png')
}, 'ERROR: failed!')


## Normalization takes place by default, and can be modified with the --normalize flag. See opts.
if (opt$normalize == 'none') {
    cat('Skipping median-normalization due to --normalize none\n')
} else {
    if (opt$normalize == 'shift') {

        tryTo('INFO: Calculating median-normalized intensities by shifting sample intensities',{
            dat.long <- shift_normalize_intensity(dat.long)
        }, 'ERROR: failed!')

        tryTo('INFO: Plotting shift-normalized intensity distributions',{
            plot_flip_beeswarm(dat.long, QC_dir, 'intensities_shift_normalized.png', plot_title='shift-normalized intensities')
        }, 'ERROR: failed!')

    } else if (opt$normalize == 'scale') {

        tryTo('INFO: Calculating median-normalized intensities by scaling sample intensities',{
            dat.long <- scale_normalize_intensity(dat.long)
        }, 'ERROR: failed!')

        tryTo('INFO: Plotting scale-normalized intensity distributions',{
            plot_flip_beeswarm(dat.long, QC_dir, 'intensities_scale_normalized.png', plot_title='scale-normalized intensities')
        }, 'ERROR: failed!')
    }

    tryTo('INFO: Re-generating wide table with normalized intensities',{
        original_colorder <- colnames(dat)
        dat <- dcast(dat.long, Protein_Group+Genes+First_Protein_Description~Sample, value.var='Intensity')
        setcolorder(dat, original_colorder)
    }, 'ERROR: failed!')
}


tryTo('INFO: Replacing NA values with 0',{
    replace_NAs(dat, colnames(dat[,-c(1:3)]), 0)
    dat.long[is.na(Intensity), Intensity := 0]
}, 'ERROR: failed!')



# pgcounts represents the distribution of Protein Groups with Intensity > 0
# Visually, it is represented as a bar plot with x=sample, y=N, ordered by descending N
# Get counts of [N=unique gene groups with `Intensity` > 0]
tryTo('INFO: Tabulating protein group counts',{
    pgcounts <- dat.long[, .N, by=Sample]
    # Order samples by ascending counts
    pgcounts[, Sample := factor(Sample, levels=pgcounts[order(-N), Sample])]
    ezwrite(pgcounts, QC_dir, 'protein_group_nonzero_counts.tsv')
    plot_pg_counts(pgcounts, QC_dir, 'protein_group_nonzero_counts.png')
}, 'ERROR: failed!')


# pgthresholds represents the decay in number of unique protein groups per sample as
# the minimum Intensity threshold is incremented. Visually represented as a line plot.
tryTo('INFO: Calculating protein group counts by minimum intensity thresholds',{
    pgthresholds <- foreach(threshold=0:(1+max(dat.long$Intensity)), .combine='rbind') %do% {
        dat.tmp <- dat.long[, list('N'=sum(Intensity > threshold)), by=Sample]
        dat.tmp[, 'Threshold' := threshold]
        return(dat.tmp)
    }
    ezwrite(pgthresholds, QC_dir, 'protein_group_thresholds.tsv')
    plot_pg_thresholds(pgthresholds, QC_dir, 'protein_group_thresholds.png')
}, 'ERROR: failed!')


tryTo('INFO: Plotting sample intensity correlations',{
    dat.correlations <- get_spearman(dat)
    ezwrite(dat.correlations, QC_dir, 'sample_correlation.tsv')
    plot_correlation_heatmap(dat.correlations, QC_dir, 'sample_correlation.png')
}, 'ERROR: failed!')


tryTo('INFO: Importing experimental design',{
    design <- fread(opt$design, header=TRUE)
    setnames(design, c('sample_name', 'condition', 'control'))
}, 'ERROR: failed!')

if (!is.null(opt$exclude)) {
    tryTo(paste('INFO: excluding samples', opt$exclude),{
        design <- design[! sample_name %in% opt$exclude]
    }, 'ERROR: failed!')
}

tryTo('INFO: Validating experimental design',{
    print(design[])
    cat('\n')
    conditions <- unique(design$condition)
    for (condition.i in conditions) {
        samples <- design[condition == condition.i, sample_name]
        control <- unique(design[condition == condition.i, control])
        if (length(control) != 1) {
            cat(paste0('ERROR: condition ', condition.i, ' maps to multiple controls: ', control, '\n'))
            cat(paste0('       Check the design matrix and esure no more than one control label per condition\n'))
            quit(exit=1)
        } else {
            cat(paste0('INFO: condition ', condition.i, ' maps to control ', control, '\n'))
        }
    }
    cat(paste0('INFO: all conditions pass check (i.e. map to one control condition)\n'))
}, 'ERROR: failed!')


## Exclude samples with N protein groups < opt$sds away from mean
## Default value: 3 standard deviations, modifiable with --sds [N]
tryTo('INFO: Identifying samples with protein group count outliers',{
    cat(paste0('INFO: defining outliers as samples with [N protein groups] > ', opt$sds, ' standard deviations from the mean\n'))
    stdev <- sd(pgcounts[,N])
    mean_count <- mean(pgcounts[,N])
    min_protein_groups <- floor(mean_count - (opt$sds * stdev))
    max_protein_groups <- ceiling(mean_count + (opt$sds * stdev))
    cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']'))
    low_count_samples <- as.character(pgcounts[N < min_protein_groups, Sample])
    high_count_samples <- as.character(pgcounts[N > max_protein_groups, Sample])
    if(length(low_count_samples)==0) {
        cat('\nINFO: No low group count samples to remove\n')
    } else {
        cat(paste0('\nINFO: Pruning low-count outlier ', low_count_samples))
        cat('\n\n')
        print(pgcounts[Sample %in% low_count_samples])
        cat('\n')
        dat[, c(low_count_samples) := NULL]    # remove sample columns from wide table
        dat.long <- dat.long[! (Sample %in% low_count_samples)] # remove rows from long table
    }
    if(length(high_count_samples)==0) {
        cat('INFO: No high group count samples to remove\n')
    } else {
        cat(paste0('\nINFO: Pruning high-count outlier ', high_count_samples))
        cat('\n')
        print(pgcounts[Sample %in% high_count_samples])
        dat[, c(high_count_samples) := NULL]    # remove sample columns from wide table
        dat.long <- dat.long[! (Sample %in% high_count_samples)] # remove rows from long table
    }
}, 'ERROR: failed!')



#### CLUSTERING ####################################################################################


# PCA
tryTo('INFO: running PCA and plotting first two components',{
    pca <- get_PCs(dat, design)
    ezwrite(pca$components, cluster_dir, 'PCA.tsv')
    ezwrite(pca$summary, cluster_dir, 'PCA_summary.tsv')
    plot_PCs(pca, cluster_dir, 'PCA.png')
}, 'ERROR: failed!')


# Hierarchical Clustering
tryTo('INFO: running Hierarchical Clustering',{
    plot_hierarchical_cluster(dat, cluster_dir, 'hierarchical_cluster.png')
}, 'ERROR: failed!')


# UMAP
tryTo('INFO: running UMAP',{
    umap <- get_umap(dat, opt$neighbors)
    ezwrite(umap, cluster_dir, 'UMAP.tsv')
    plot_umap(umap, cluster_dir, 'UMAP.png')
}, 'ERROR: failed!')



#### DIFFERENTIAL INTENSITY ########################################################################

tryTo('INFO: Re-generating raw intensities from normalized log values for T-test',{
    # Transform from log to linear
    dat.raw <- convert_log_to_raw(dat, opt$log_base)
    # Convert 1 back to 0
    replace_values(dat.raw, colnames(dat.raw[,-c(1:3)]), 1, 0)
}, 'ERROR: failed!')

tryTo('INFO: Running differential intensity t-tests',{
    for (treatment in conditions) {
        treatment_sample_names <- intersect(colnames(dat), design[condition == treatment, sample_name])
        if (length(treatment_sample_names)==0) {next}
        control <- unique(design[condition == treatment, control])
        control_sample_names <- colnames(dat)[colnames(dat) %like% control]
        if (length(control_sample_names)==0) {next}
        if(treatment != control) {
            t_test <- do_t_test(dat.raw, treatment_sample_names, control_sample_names)
            ezwrite(t_test[order(p.adj)], DI_dir, paste0(treatment, '-vs-', control, '.tsv'))
            plot_volcano(t_test, treatment, control, opt$log_base, opt$lfc_threshold, opt$fdr_threshold, DI_dir, opt$labelgene)
        }
    }
}, 'ERROR: failed!')

quit()

## TODO:
# Evaluate normalization with samples from very different batches?
# Calculate differential intensities using DESeq2?
# Function to compare q-values from spectronaut VS DIA-NN?
