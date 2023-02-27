#!/usr/bin/env Rscript
# R/4
#proteomics analysis for DIA-NN and Spectronaut quantity estimates

#### PACKAGES ######################################################################################
package_list = c('ggplot2', 'data.table', 'corrplot', 'umap', 'magick',
                 'ggbeeswarm', 'ggrepel', 'optparse', 'ggthemes', 'foreach')
cat("INFO: Loading required packages\n      ")
cat(paste(package_list, collapse='\n      ')); cat('\n')

defaultW <- getOption("warn"); options(warn = -1)   # Temporarily disable warnings for quiet loading
if(all((lapply(package_list, require, character.only=TRUE)))) {
    cat("INFO: All packages successfully loaded\n")
} else {
    cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}
options(warn = defaultW)    # Turn warnings back on

#### DEFAULTS ######################################################################################
intensity_min_threshold <- 0


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
        "--base",
        dest="log_base",
        default=10, 
        help='Base for log transformation of intensity data. Default: 10'
    ),
    make_option(
        "--normalize",
        default='median',
        type='character',
        help=paste(
            'shift: adjust sample intensities to match global median by adding a constant',
            'scale: adjust sample intensities to match global median by multiplicative scaling',
            'none: do not normalize',
            sep=optparse_indent
        )
    ),
    make_option(
        "--sds",
        dest = 'sds',
        default=3,
        type='numeric',
        help=paste(
            'Filter out samples with protein group counts > N standard deviations from the mean.',
            'Default: 3',
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
print(opt)


if(TRUE) {
    opt <- list(); 
    opt$pgfile <- 'output/report.pg_matrix.tsv'
    opt$outdir <- 'output'
    opt$design <- 'example/design_matrix.csv'
    opt$log_base <- 2
    opt$dry <- FALSE
    opt$normalize <- 'shift'
}


if(opt$dry) {
    cat("INFO: Quitting due to --dry run\n")
    quit(status=0)
}

if(! opt$normalize %in% c('shift','scale','none')) {
    cat("ERROR: --normalize must be 'shift' 'scale' or 'none'\n")
    quit(status=1)
}


#### FUNCTIONS #####################################################################################
# Within data.table `DT`, 
# for `sdcols` specified columns, 
# replaces all NA with `newvalue`
replace_NAs <- function(DT, sdcols, newvalue) {
  DT[, (sdcols) := lapply(.SD, function(x) {ifelse(is.na(x),newvalue,x)}), .SDcols=sdcols]
}

tryTo <- function(infomsg='', cmd,  errormsg='ERROR: failed!') {
    # tryTo is a simple tryCatch wrapper, taking a command + message.
    # If the try block fails, an error is printed and R quits.
    cat(paste0(infomsg, '\n'))
    tryCatch(cmd, 
        error = function(e){
            cat(paste0(errormsg, '\n', e))
            quit(status=1)
        },
        finally = cat('')
    )
}

ezwrite <- function(x, output_dir, output_filename) {
    # Wrapper for fwrite that uses standard TSV output defaults.
    # Concatenates output directory and filename for final output location.
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    fwrite(x, file=paste0(output_dir, '/', output_filename),
        quote=F,
        row.names=F,
        col.names=T,
        sep='\t')
    
}

standardize_format <- function(DT.original) {
    DT <- copy(DT.original)
    if("Protein.Ids" %in% colnames(DT)) {
        DT[, 'Protein.Ids' := NULL]
        DT[, 'Protein.Names' := NULL]
        setnames(DT, 'Protein.Group', 'Protein_Group')
        setnames(DT, 'First.Protein.Description', 'First_Protein_Description')
    } else if('PG.ProteinGroups' %in% colnames(DT)) {
        setnames(DT, 'PG.ProteinGroups', 'Protein_Groups')
        setnames(DT, 'PG.Genes', 'Genes')
        setnames(DT, 'PG.ProteinDescriptions', 'First_Protein_Description')
        # Use only first protein description
        DT[, 'First_Protein_Description' := tstrsplit(First_Protein_Description, split=';')[1]]
    }

    # Remove leading directories for sample names
    # e.g. /path/to/sample1.mzML -> sample1.mzML
    setnames(DT, basename(colnames(DT)))

    # Remove trailing file extensions
    extensions <- '.mzML$|.mzml$|.RAW$|.raw$|.dia$|.DIA$'
    extension_samplenames <-  colnames(DT)[colnames(DT) %like% extensions]
    trimmed_samplenames <- gsub(extensions, '', extension_samplenames)
    setnames(DT, extension_samplenames, trimmed_samplenames)
    return(DT[])
}


melt_table <- function(DT) {
    info_cols <- c('Protein_Groups', 'Genes', 'First_Protein_Description')
    DT.long <- melt(DT, 
    measure.vars=colnames(DT[,-c(1:3)]),
    variable.name='Sample',
    value.name='Intensity')
    return(DT.long)
}

filter_intensity <- function(DT, threshold) {
    return(DT[Intensity > threshold])
}

plot_pg_counts <- function(DT, output_dir, output_filename) {
    n_samples <- nrow(DT)
    g <- ggplot(DT, aes(x=Sample, y=N)) +
                geom_bar(stat='identity', position='dodge', aes(y=N)) +
                geom_text(aes(label=N, y=N + (0.05*max(pgcounts$N)))) +
                coord_flip() +
                theme_few() +
                labs(x='Sample', y=paste0('N Protein Groups where Log[', opt$log_base, '](Intensity)> 0'))
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
   
    ggsave(g,filename=paste0(output_dir, output_filename), width=20, height=2*n_samples, units='cm')
}

plot_pg_thresholds <- function(DT, output_dir, output_filename) {
    # G
    g <- ggplot(DT, aes(x=Threshold, y=N, color=Sample)) +
            geom_line() +
            geom_point(shape=21, alpha=0.5) +
            theme_few() +
            labs(x=paste0('Minimum Log[', opt$log_base, '](Intensity) Threshold'),
            y=paste0('N Protein Groups where Log[', opt$log_base, '](Intensity)> 0')) 
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g,
        filename=paste0(output_dir, output_filename),
        width=20,
        height=20,
        units='cm'
    )
}

plot_density <- function(DT.original, output_dir, output_filename) {
    DT <- copy(DT.original)
    intensity_median <- median(DT[,Intensity])
    n_samples <- length(unique(DT[,Sample]))
    dat.quantiles <- DT[, list(
                    'q025'=quantile(Intensity, 0.025),
                    'q25'=quantile(Intensity, 0.25),
                    'q50'=quantile(Intensity, 0.50),
                    'q75'=quantile(Intensity, 0.75),
                    'q975'=quantile(Intensity, 0.975)
                    ), by=Sample]

    dat.legend <- melt(dat.quantiles[which.min(as.numeric(Sample))], measure.vars=colnames(dat.quantiles[,-1]))
    dat.legend[, qlabel := tstrsplit(variable, split='q')[2]]
    dat.legend[, qlabel := paste0('0.', qlabel)]
    dat.legend[, qlabel := as.numeric(qlabel)]

    g <- ggplot(DT, linetype='solid', aes(x=Intensity)) +
        geom_density(fill='gray80') +
        theme_few() +
        facet_grid(Sample~., switch='y') +
        geom_vline(xintercept=intensity_median, color='red') +
        geom_vline(data=dat.quantiles, linetype='solid',  alpha=0.7, aes(xintercept=q50))+
        geom_vline(data=dat.quantiles, linetype='dashed', alpha=0.7,  aes(xintercept=q25))+
        geom_vline(data=dat.quantiles, linetype='dashed', alpha=0.7,  aes(xintercept=q75))+
        geom_vline(data=dat.quantiles, linetype='dotted', alpha=0.7,  aes(xintercept=q025))+
        geom_vline(data=dat.quantiles, linetype='dotted', alpha=0.7,  aes(xintercept=q975))+
        theme(strip.text.y.left = element_text(angle = 0, hjust=0.5, vjust=0.5)) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) +
        labs(x=paste0('Log[', opt$log_base, '](Intensity)'), title='Intensity Distribution across Samples') +
        geom_label(data=dat.legend, aes(x=value, y=0.285, label=qlabel)) +
        theme(panel.border = element_blank()) +
        ylim(0,0.3)
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g, 
        filename=paste0(output_dir, output_filename),
        height=2.5*n_samples,
        width=30,
        units='cm'
    )
}

shift_normalize_intensity <- function(DT.original) {
    DT <- copy(DT.original)
    # Get global median of intensity values
    global_median <- median(DT[, Intensity])
    DT[, 'sample_median' := median(Intensity), by=Sample]
    DT[, 'global_median' := global_median]
    DT[, 'NormInt' := Intensity - (sample_median - global_median)]
    DT[, c('Intensity', 'sample_median', 'global_median') := NULL]
    setnames(DT, 'NormInt', 'Intensity')
    return(DT[])
}

scale_normalize_intensity <- function(DT.original) {
    DT <- copy(DT.original)
    # Get global median of intensity values
    global_median <- median(DT[, Intensity])
    DT[, 'sample_median' := median(Intensity), by=Sample]
    DT[, 'global_median' := global_median]
    DT[, 'NormInt' := Intensity * (global_median / sample_median)]
    DT[, c('Intensity', 'sample_median', 'global_median') := NULL]
    setnames(DT, 'NormInt', 'Intensity')
    return(DT[])
}

plot_flip_beeswarm <- function(DT, output_dir, output_filename, plot_title) {
    n_samples <- length(unique(DT$Sample))
    intensity_median <- median(DT[, Intensity])
    dat.quantiles <- DT[, list(
                'q025'=quantile(Intensity, 0.025),
                'q25'=quantile(Intensity, 0.25),
                'q50'=quantile(Intensity, 0.50),
                'q75'=quantile(Intensity, 0.75),
                'q975'=quantile(Intensity, 0.975)
                ), by=Sample]

    dat.legend <- melt(dat.quantiles[which.min(as.numeric(Sample))], measure.vars=colnames(dat.quantiles[,-1]))
    dat.legend[, qlabel := tstrsplit(variable, split='q')[2]]
    dat.legend[, qlabel := paste0('0.', qlabel)]
    dat.legend[, qlabel := as.numeric(qlabel)]

    g <- ggplot(DT, aes(x=0, y=Intensity)) +
        geom_beeswarm(alpha=0.5, shape=21) +
        facet_grid(.~Sample, switch='x') +
        theme_few() +
        geom_hline(data=dat.quantiles, color='gray50', linetype='solid',  alpha=0.7, aes(yintercept=q50))+
        geom_hline(data=dat.quantiles, color='gray50', linetype='dashed', alpha=0.7,  aes(yintercept=q25))+
        geom_hline(data=dat.quantiles, color='gray50', linetype='dashed', alpha=0.7,  aes(yintercept=q75))+
        geom_hline(data=dat.quantiles, color='gray50', linetype='dotted', alpha=0.7,  aes(yintercept=q025))+
        geom_hline(data=dat.quantiles, color='gray50', linetype='dotted', alpha=0.7,  aes(yintercept=q975))+
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        labs(title=plot_title, x='Sample', y=paste0('Log[', opt$log_base, '](Intensity)')) +
        theme(strip.text.x.bottom = element_text(angle = 90, hjust=0.5, vjust=0.5)) +
        scale_y_continuous(position='right') +
        theme(axis.text.y.right=element_text(angle=90, hjust=0.5)) +
        theme(axis.title.x.bottom = element_text(angle = 180, hjust=0.5, vjust=0.5)) +
        theme(axis.title.y.right = element_text(angle = 90, hjust=0.5, vjust=0.5)) +
        geom_hline(yintercept=intensity_median, color='red', linetype='dashed', alpha=1)
    ggsave(g, filename='.beeswarm.tmp.png', height=15, width=2*n_samples, units='cm')
    tmpimage <- image_read('.beeswarm.tmp.png')
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    image_rotate(tmpimage, 90) %>% image_write(paste0(output_dir, output_filename))
    invisible(file.remove('.beeswarm.tmp.png')) # invisible suppresses 'TRUE' being printed

}

plot_correlation_heatmap <- function(DT.corrs, output_dir, output_filename) {
    n_samples <- length(unique(DT.corrs[,SampleA]))
    max_limit <- max(DT.corrs$Spearman)
    min_limit <- min(DT.corrs$Spearman)
    mid_limit <- as.numeric(format(((max_limit + min_limit) / 2), digits=3))
    g <- ggplot(DT.corrs, aes(x=SampleA, y=SampleB, fill=Spearman, label=Spearman)) +
    geom_tile() +
    geom_text(color='gray10') + 
    theme_few() +
    scale_fill_gradient2(low = "skyblue", high = "tomato1", mid = "white", 
                        midpoint = mid_limit, limit = c(min_limit,max_limit),
                        space = "Lab", breaks=c(min_limit, mid_limit, max_limit),
                        name="Spearman\nCorrelation\n") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g,
    filename=paste0(output_dir, output_filename),
    height=1.2*n_samples,
    width=2*n_samples,
    units='cm')

}


get_correlations <- function(DT.original) {
    DT <- copy(DT.original)
    #### Pairwise correlations between sample columns
    dt.samples <- DT[,-c(1:3)]     # Ignore info columns (subset to only intensity values)
    dt.corrs <- cor(as.matrix(na.omit(dt.samples)), method='spearman')  

    # Convert to lower triangle only
    # (no need for full correlation matrix with diagonal or repeated comparisons)
    dt.corrs[lower.tri(dt.corrs, diag=T)] <- NA
    dt.corrs <- as.data.table(dt.corrs, keep.rownames=T)
    dt.corrs <- melt(dt.corrs, measure.vars=dt.corrs[,rn], value.name='Spearman')
    dt.corrs <- dt.corrs[! is.na(Spearman)]
    setnames(dt.corrs, c('rn', 'variable'), c('SampleA','SampleB'))

    # Format correlations as 3 digits
    dt.corrs[, Spearman := as.numeric(format(Spearman, digits=3))]

    # Adjust levels such that both axes plot samples in the same order
    dt.corrs.levels <- sort(as.character(unique(dat.long$Sample)))
    dt.corrs[, SampleA := factor(SampleA, levels=dt.corrs.levels)]
    dt.corrs[, SampleB := factor(SampleB, levels=dt.corrs.levels)]
    return(dt.corrs[])
}

get_PCs <- function(DT) {
    out <- list()
    cluster.dat <- DT[,-c(1:3)]   # Subset to only sample intensity columns
    replace_NAs(cluster.dat, colnames(cluster.dat), 0)    # Replace NA values with 0
    cluster.dat <- t(cluster.dat)
    pca <- prcomp(cluster.dat, center = TRUE, scale. = TRUE)
    out$summary <- as.data.table(t(summary(pca)$importance), keep.rownames=T)
    setnames(out$summary, c('component','stdv','percent','cumulative'))
    pca <- as.data.frame(pca$x)
    pca <- setDT(pca, keep.rownames=T)[]
    setnames(pca, 'rn', 'Sample')
    out$components <- merge(pca, design, by.x= 'Sample', by.y='sample_name')
    return(out)
}

plot_PCs <- function(PCA, output_dir, output_filename) {
    PCA$summary
    PCA$components

    pc1_label <- as.character(format(100*PCA$summary[component=='PC1', 'percent'], digits=3))
    pc1_label <- paste0('PC1: ', pc1_label, '%')

    pc2_label <- as.character(format(100*PCA$summary[component=='PC2', 'percent'], digits=3))
    pc2_label <- paste0('PC2: ', pc2_label, '%')

    g <- ggplot(PCA$components, aes(x=PC1, y=PC2, color=condition)) +
        geom_point() +
        stat_ellipse(level=0.95) +
        theme_few() +
        labs(x=pc1_label, y=pc2_label)

    cat(paste0('   -> ', output_dir, output_filename, '\n'))
   
    ggsave(g,filename=paste0(output_dir, output_filename), width=18, height=18, units='cm')
}

#### MAKE DIRS #####################################################################################
QC_dir <- paste0(opt$outdir, '/QC/')
if(! dir.exists(QC_dir)){
    dir.create(QC_dir, recursive = T)
}

cluster_dir <- paste0(opt$outdir, '/clustering/')
if(! dir.exists(cluster_dir)){
    dir.create(cluster_dir, recursive = T)
}

DA_dir <- paste0(opt$outdir, '/differential_intensity/')
if(! dir.exists(DA_dir)){
    dir.create(DA_dir, recursive = T)
}

#### IMPORT DATA ###################################################################################


tryTo(paste0('INFO: Reading input file ', opt$pgfile),{
    dat <- fread(opt$pgfile)
}, paste0('ERROR: problem trying to load ', opt$pgfile, ', does it exist?'))

tryTo(paste0('INFO: Massaging data from ', opt$pgfile, ' into a common style format for processing'), {
    dat <- standardize_format(dat)
}, 'ERROR: failed! Check for missing/corrupt headers?')

tryTo(paste0('INFO: Converting to long format'), {
    dat.long <- melt_table(dat)
}, 'ERROR: failed! Check for missing/corrupt headers?')

tryTo(paste0('INFO: Applying Log[base', opt$log_base, '](value+1) transformation to intensities'),{
    dat.long <- dat.long[, Intensity := log((Intensity + 1), base=opt$log_base)]
},'ERROR: failed! Was your log base numeric and > 1?')

tryTo('INFO: Excluding all unquantified or zero intensities', {
    dat.long <- dat.long[! is.na(Intensity)][Intensity != 0]
}, 'ERROR: failed!')



#### QC ############################################################################################

## Filtering
increasing_levels <- dat.long[, list('median'=median(Intensity)), by=Sample][order(median), Sample]
# decreasing_levels <- dat.long[, list('median'=median(Intensity)), by=Sample][order(-median), Sample]
dat.long[, Sample := factor(Sample, levels=increasing_levels)]

tryTo(paste0('INFO: Applying Filter Log[', opt$log_base, '](Intensity) > ',intensity_min_threshold),{
    dat.long <- filter_intensity(dat.long, intensity_min_threshold)
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
    dat.correlations <- get_correlations(dat)
    ezwrite(dat.correlations, QC_dir, 'sample_correlation.tsv')
    plot_correlation_heatmap(dat.correlations, QC_dir, 'sample_correlation.png')
}, 'ERROR: failed!')



tryTo('INFO: Importing experimental design',{
    design <- fread(opt$design)
    print(design[])
    conditions <- unique(design$condition)
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
        cat('\n')
        print(pgcounts[Sample %in% low_count_samples])
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


quit()



#### CLUSTERING ####################################################################################

prcomp(t(as.matrix(dat[,-c(1:3)])), center = TRUE, scale. = TRUE)

t(as.matrix(dat[,-c(1:3)]))

## PCA

## UMAP

## CLUSTER





tryTo('INFO: running PCA',{
    pca <- get_PCs(dat)
    plot_PCs(pca, cluster_dir, 'PCA.png')
}, 'ERROR: failed!')



tryTo('INFO: running UMAP ',{

}, 'ERROR: failed!')

tryTo('INFO: performing clustering',{

}, 'ERROR: failed!')



##cluster data(na=0)
cluster_data=pro[,grep('mzML',colnames(pro))]
cluster_data[is.na(cluster_data)]=0
cluster_data=cluster_data[which(rowSums(cluster_data)>0),]
log2_cluster_data=log2(cluster_data+1)

##PCA and plot
pca_data=t(log2_cluster_data)
pca=prcomp(pca_data, center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
pca_df = as.data.frame(pca$x)
pca_df$Condition=gsub('_[1-9]*.mzML$','',rownames(pca_df))
percentage=round(summary(pca)$importance[2,]*100, digits = 2)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#1E90FF",
            "#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE",
            "#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5",
            "#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C",
            "#FFFFE0","#EE82EE","#FF6347","#6A5ACD",
            "#9932CC","#8B008B","#8B4513","#DEB887")
p=ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+
  xlab(paste0("PC1","(",percentage[1],"%)")) +
  ylab(paste0("PC2","(",percentage[2],"%)"))+
  scale_color_manual(values = allcolour)+
  theme_classic()+
  stat_ellipse(level=0.95)
ggsave(plot = p,filename = paste0(out_dir,opt$prefix,"_PCA_circle.pdf"),height = 5,width = 7)
p=ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)+
  xlab(paste0("PC1","(",percentage[1],"%)")) +
  ylab(paste0("PC2","(",percentage[2],"%)"))+
  scale_color_manual(values = allcolour)+
  theme_classic()
ggsave(plot = p,filename = paste0(out_dir,opt$prefix,"_PCA.pdf"),height = 5,width = 7)


##hc
dist_mat <- dist(t(log2_cluster_data)) #
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- gsub('.mzML$','',colnames(log2_cluster_data))
pdf(file =paste0(out_dir,opt$prefix,"_hc_cluster_log2.pdf"))
plot(hc_cluster,cex=0.8,col="dark red",labels = samplesname,main="HC Cluster")
dev.off()
##mean to one
data_in1 <- data.frame(row.names = rownames(log2_cluster_data))
for (i in unique(gsub("_[1-9]*.mzML$",'',colnames(log2_cluster_data)))) {
  tmp <- log2_cluster_data[,grep(i, colnames(log2_cluster_data)),]
  tmp_i <- data.frame(rowMeans(tmp),row.names = rownames(log2_cluster_data))
  colnames(tmp_i) <-i
  data_in1 <- cbind(data_in1,tmp_i)
}
dist_mat <- dist(t(data_in1)) #
hc_cluster <- hclust(dist_mat,method = "complete")
samplesname<- colnames(data_in1)
pdf(file = paste0(out_dir,opt$prefix,"_hc_cluster_log2_mean.pdf"))
plot(hc_cluster,cex=0.8,col="dark red",labels = samplesname,main="HC Cluster")
dev.off()

##umap
if (length(unique(gsub('_[1-9]*.mzML$','',colnames(log2_cluster_data)))) > 6) {
  umap_data <- t(log2_cluster_data)
  set.seed(100)
  umap_out <- umap(umap_data)
  umap_df <- as.data.frame(umap_out$layout) 
  colnames(umap_df) <-c("UMAP1","UMAP2")
  umap_df$Condition <- gsub('_[1-9]*.mzML$$','',rownames(umap_df))
  
  p=ggplot(umap_df) +
    geom_point(aes(x=UMAP1, y=UMAP2, color=Condition),size=4)+
    scale_color_manual(values = allcolour)+
    theme_classic()
  ggsave(plot = p,filename = paste0(out_dir,opt$prefix,"umap.pdf"),height = 5,width = 7)
}

#DE analysis--------------------------------------------------------------------
##mkdir

if (!dir.exists(paste0(opt$outdir,"/DE_analysis/"))){
  dir.create(paste0(opt$outdir,"/DE_analysis/"),recursive = T)
}
out_dir=paste0(opt$outdir,"/DE_analysis/")



if (!is.null(opt$design_matrix)){
  ##read files
  design_matrix=read.csv(opt$design_matrix)
  ##parameters 
  fdr_cutoff=0.05
  lfc_cutoff=0.585
  ##ttest
  condition=unique(design_matrix$condition)
  for (i in condition) {
    data_i=pro[,grep(i,colnames(pro))]
    data_control=pro[,grep(unique(design_matrix$control[which(design_matrix$condition == i)]),colnames(pro))]
    data=cbind(data_i,data_control)
    data[is.na(data)]=0
    rownames(data)=pro$Protein.Group
    rm=apply(data, 1, function(x){
      sum(x == 0) > ncol(data)/2
      })
    df=data[!rm,]
    df=df[apply(df,1, var) != 0, ]
    df=df[apply(df[,grep(i,colnames(df))],1, var) != 0, ]
    df=df[apply(df[,grep(unique(design_matrix$control[which(design_matrix$condition == i)]),colnames(df))],1, var) != 0, ]
    pvalue=apply(df, 1, function(x){
      a =factor(c(rep('treat',ncol(data_i)),
                rep("control",ncol(data_control))),
              levels = c('treat',"control"))
      fvalue=var.test(x~a)
      if (!is.na(fvalue$p.value)){ 
        if (fvalue$p.value > 0.05){
        t.test(x~a, var.equal = T)
          }else{
            t.test(x~a, var.equal = F)
            }}
      })
    result_ttest=data.frame(ID=names(pvalue), 
                          Pvalue = as.numeric(unlist(lapply(pvalue,function(x) x$p.value))),
                          log2FC = log2(as.numeric(unlist(lapply(pvalue,function(x) x$estimate[1]/(x$estimate[2]+1))))))
    result_ttest$adj.Pvalue=p.adjust(result_ttest$Pvalue, method = 'BH', n = length(result_ttest$Pvalue))
    result_ttest= merge(df,result_ttest,by.x =0,by.y =1,all=F)
    result_ttest =merge(pro[,-grep("mzML",colnames(pro))],result_ttest,by.x ='Protein.Group',by.y=1,all=F)
    result_ttest=result_ttest[order(result_ttest$log2FC,decreasing = T),]
    write.csv(result_ttest,file = paste0(out_dir,i,'_',unique(design_matrix$control[which(design_matrix$condition==i)]),'_ttest.csv'),row.names = F)
    result_ttest <- na.omit(result_ttest)
    options(ggrepel.max.overlaps=Inf)
    vol_plot=result_ttest
    vol_plot$Group <- "Others"
    vol_plot$Group[which(vol_plot$log2FC >= lfc_cutoff)] <-"UP"
    vol_plot$Group[which(vol_plot$log2FC <= -lfc_cutoff)] <-"DOWN"
    vol_plot$Group[which(vol_plot$adj.Pvalue >= fdr_cutoff)]<- "Others"
    up_gene_5 <- vol_plot[which(vol_plot$Group=="UP"),]
    up_gene_5 <- up_gene_5[order(up_gene_5$log2FC,decreasing = T)[1:5],]
    down_gene_5 <- vol_plot[which(vol_plot$Group=="DOWN"),]
    down_gene_5 <- down_gene_5[order(down_gene_5$log2FC,decreasing = F)[1:5],]
    top5_gene <- rbind(up_gene_5,down_gene_5)
    top5_gene <- top5_gene[!duplicated(top5_gene$Genes),]
    p=ggplot(vol_plot, aes(x = log2FC, y = -log10(adj.Pvalue))) +
      geom_point(aes(color = Group)) +
      scale_color_manual(values = c("blue", "grey","red"))  +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") +
      geom_label_repel(
        data = subset(top5_gene),
        aes(label = Genes),
        size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"))+
      geom_hline(yintercept=-log10(fdr_cutoff), linetype="dashed")+ 
      geom_vline(xintercept=lfc_cutoff, linetype="dashed")+ 
      geom_vline(xintercept=-lfc_cutoff, linetype="dashed")+
      theme_classic()
  
      ggsave(file = paste0(out_dir,i,"_",unique(design_matrix$control[which(design_matrix$condition==i)]),"_top10_fdr0.05_fc1.5_vocal.pdf"),plot = p,width = 8,height = 8)
      
  }
  }

