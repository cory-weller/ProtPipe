
#### FUNCTIONS #####################################################################################
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



replace_NAs <- function(DT, sdcols, newvalue) {
    # Within data.table `DT`, 
    # for `sdcols` specified columns, 
    # replaces all NA with `newvalue`
    DT[, (sdcols) := lapply(.SD, function(x) {ifelse(is.na(x),newvalue,x)}), .SDcols=sdcols]
}

replace_values <- function(DT, sdcols, oldvalue, newvalue) {
    # Within data.table `DT`, 
    # for `sdcols` specified columns, 
    # replaces all NA with `newvalue`
    DT[, (sdcols) := lapply(.SD, function(x) {ifelse(x==oldvalue,newvalue,x)}), .SDcols=sdcols]
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
    # Accepts an input protein group intensity data.table, whether spectronaut or DIA-NN format,
    # and restructures into one consistent style for downstream processing
    DT <- copy(DT.original)
    if("Protein.Ids" %in% colnames(DT)) {
        DT[, 'Protein.Ids' := NULL]
        DT[, 'Protein.Names' := NULL]
        setnames(DT, 'Protein.Group', 'Protein_Group')
        setnames(DT, 'First.Protein.Description', 'First_Protein_Description')
    } else if('PG.ProteinGroups' %in% colnames(DT)) {
        setnames(DT, 'PG.ProteinGroups', 'Protein_Group')
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

trim_colnames <- function(DT) {
    colnames_out <- gsub(pattern="\\[.*\\] ", replacement='', x=colnames(DT))   # trim leading [N] 
    colnames_out <- gsub(pattern="\\..*\\.PG\\.Quantity", replacement='', x=colnames_out)   # remove suffix
    return(colnames_out)
}

melt_intensity_table <- function(DT) {
    # Converts intensity data.table to long format
    # info_cols <- c('Protein_Group', 'Genes', 'First_Protein_Description')
    DT.long <- melt(DT, 
    measure.vars=colnames(DT[,-c(1:3)]),
    variable.name='Sample',
    value.name='Intensity')
    return(DT.long)
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
    # Currently UNUSED as the beeswarm function serves the purpose well
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

convert_log_to_raw <- function(DT.original, log_base) {
    DT <- copy(DT.original)
    samplenames <- colnames(DT[,-c(1:3)])
    DT[, (samplenames) := lapply(.SD, function(x) log_base^x), .SDcols=samplenames]
    return(DT[])
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


get_spearman <- function(DT.original) {
    DT <- copy(DT.original)
    #### Pairwise correlations between sample columns
    dt.samples <- DT[,-c(1:3)]     # Ignore info columns (subset to only intensity values)
    dt.corrs <- cor(as.matrix(na.omit(dt.samples)), method='spearman')  

    # Convert to lower triangle only
    # (no need for full correlation matrix with diagonal or repeated comparisons)
    dt.corrs[lower.tri(dt.corrs, diag=T)] <- NA
    dt.corrs <- as.data.table(dt.corrs, keep.rownames=T)
    lvls <- dt.corrs[,rn]
    dt.corrs <- melt(dt.corrs, measure.vars=dt.corrs[,rn], value.name='Spearman')
    dt.corrs <- dt.corrs[! is.na(Spearman)]
    setnames(dt.corrs, c('rn', 'variable'), c('SampleA','SampleB'))

    # Format correlations as 3 digits
    dt.corrs[, Spearman := as.numeric(format(Spearman, digits=3))]

    # Adjust levels such that both axes plot samples in the same order
    dt.corrs[, SampleA := factor(SampleA, levels=lvls)]
    dt.corrs[, SampleB := factor(SampleB, levels=lvls)]
    return(dt.corrs[])
}

get_pearson_matrix <- function(DT.original) {
    DT <- copy(DT.original)
    #### Pairwise correlations between sample columns
    dt.samples <- DT[,-c(1:3)]     # Ignore info columns (subset to only intensity values)
    dt.corrs <- cor(as.matrix(na.omit(dt.samples)), method='pearson')  

    dt.corrs <- as.data.table(dt.corrs, keep.rownames=T)
    #dt.corrs <- melt(dt.corrs, measure.vars=dt.corrs[,rn], value.name='Spearman')
    #dt.corrs <- dt.corrs[! is.na(Spearman)]
    #setnames(dt.corrs, c('rn', 'variable'), c('SampleA','SampleB'))

    # Format correlations as 3 digits
    dt.corrs[, Spearman := as.numeric(format(Spearman, digits=3))]

    # Adjust levels such that both axes plot samples in the same order
    dt.corrs.levels <- sort(as.character(unique(dat.long$Sample)))
    dt.corrs[, SampleA := factor(SampleA, levels=dt.corrs.levels)]
    dt.corrs[, SampleB := factor(SampleB, levels=dt.corrs.levels)]
    return(dt.corrs[])
}


get_PCs <- function(DT.intensity, DT.design) {
    out <- list()
    # Formatting prior to running PCA
    cluster.dat <- copy(DT.intensity[,-c(1:3)]) 
    samplenames <- copy(colnames(cluster.dat))
    replace_NAs(cluster.dat, samplenames, 0)    # Replace NA values with 0
    cluster.dat[, 'N_notzero' := apply(.SD, 1, function(x) sum(x!=0)), .SDcols=samplenames]
    cluster.dat <- cluster.dat[N_notzero!=0][,samplenames, with=F]

    cluster.dat <- t(cluster.dat)                       # Transpose before PCA
    pca <- prcomp(cluster.dat, center = TRUE, scale. = TRUE)

    # Format summary output
    out$summary <- as.data.table(t(summary(pca)$importance), keep.rownames=T)
    setnames(out$summary, c('component','stdv','percent','cumulative'))

    # Format PCA output
    pca <- as.data.frame(pca$x)
    pca <- setDT(pca, keep.rownames=T)[]
    setnames(pca, 'rn', 'Sample')

    # Merge in design matrix to add 'condition' column
    out$components <- merge(pca, DT.design, by.x= 'Sample', by.y='sample_name')
    return(out)
}


plot_PCs <- function(PCA, output_dir, output_filename) {
    # Get % explained from PCA$summary table for PC1 and PC2
    pc1_label <- as.character(format(100*PCA$summary[component=='PC1', 'percent'], digits=3))
    pc1_label <- paste0('PC1: ', pc1_label, '%')

    pc2_label <- as.character(format(100*PCA$summary[component=='PC2', 'percent'], digits=3))
    pc2_label <- paste0('PC2: ', pc2_label, '%')

    # Plot with 95% CI ellipse
    g <- ggplot(PCA$components, aes(x=PC1, y=PC2, color=condition)) +
        geom_point() +
        stat_ellipse(level=0.95) +
        theme_few() +
        labs(x=pc1_label, y=pc2_label)

    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g,filename=paste0(output_dir, output_filename), width=18, height=18, units='cm')
}


plot_hierarchical_cluster <- function(DT, output_dir, output_filename) {
    dist_mat <- dist(t(DT[,-c(1:3)]))
    n <- ncol(DT)-3
    hc <- hclust(dist_mat, method = "complete")
    g <- ggdendrogram(hc, rotate=TRUE) + labs(title='Hierarchical clustering')

    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g, filename=paste0(output_dir, output_filename), height=n, width=15, units='cm')
}


get_umap <- function(DT.original, neighbors) {
    DT <- t(DT.original[,-c(1:3)])
    set.seed(100)
    DT.umap <- umap(DT, n_neighbors=neighbors)
    DT.out <- as.data.table(DT.umap$layout, keep.rownames=TRUE)
    setnames(DT.out, c('Sample', 'UMAP1', 'UMAP2'))
    DT.out <- merge(DT.out, design, by.x='Sample', by.y='sample_name')
    return(DT.out[])
}


plot_umap <- function(DT, output_dir, output_filename) {
    g <- ggplot(DT, aes(x=UMAP1, y=UMAP2, color=condition)) +
    geom_point() +
    theme_few() 
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g, filename=paste0(output_dir, output_filename), height=15, width=20, units='cm')
}

# filter_pgs <- function(DT.all, treatment_sample_names, treatment, control)  {
#     control_sample_names <- colnames(DT.all)[colnames(DT.all) %like% control]
#     n_treatment <- length(treatment_sample_names)
#     n_controls <- length(control_sample_names)

#     info_cols <- colnames(DT.all)[1:3]
#     DT <- DT.all[, c(info_cols, treatment_sample_names, control_sample_names), with=F]

#     # Filter: protein groups with at least half of control and treatment samples having non-zero value
#     DT[, 'N_missing_treatment' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(treatment_sample_names)]
#     DT[, 'N_missing_control' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(control_sample_names)]
#     DT <- DT[N_missing_treatment < (n_treatment/2)]
#     DT <- DT[N_missing_control < (n_controls/2)]
#     DT[, 'N_missing_treatment' := NULL]
#     DT[, 'N_missing_control' := NULL]
#     return(DT[])
# }


do_t_test <- function(DT, treatment_samples, control_samples) {

    dat <- copy(DT)
    n_treatment <- length(treatment_samples)
    n_control <- length(control_samples)

    # Retain first three columns plus all treatment and control columns
    dat <- dat[,c(colnames(DT)[1:3], treatment_samples, control_samples), with=F]
    
    # Convert NA to 0
    dat[is.na(dat)] <- 0

    # Drop rows (protein groups) with > 50% missingness in treatment OR control group
    dat[, 'missing_treatment' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=treatment_samples]
    dat[, 'missing_control' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=control_samples]
    dat <- dat[missing_treatment <= (n_treatment/2)]
    dat <- dat[missing_control <= (n_control/2)]
    dat[, c('missing_treatment','missing_control') := NULL]


    # Drop rows (protein groups) with 0 variance in treatment OR control group
    dat[, 'control_var' := apply(.SD, 1, var), .SDcols=c(control_samples)]
    dat[, 'treatment_var' := apply(.SD, 1, var), .SDcols=c(treatment_samples)]
    dat <- dat[control_var != 0]
    dat <- dat[treatment_var != 0]
    dat[, c('control_var','treatment_var') := NULL]

    # Perform t-test  on treatment and control columns
    t_test <- apply(dat[,-c(1:3)], 1, function(x){
        a =factor(c(rep('treatment',n_treatment),
                    rep("control",n_control)),
                    levels = c('treatment',"control"))
        fvalue=var.test(x~a)
        if (!is.na(fvalue$p.value)){ 
            if (fvalue$p.value > 0.05){
                result <- t.test(x~a, var.equal = T)
            }else{
                result <- t.test(x~a, var.equal = F)
            }
        }
        treatment_estimate <- as.numeric(unlist(result$estimate[1]))
        control_estimate <- as.numeric(unlist(result$estimate[2]))
        return(data.table('P'=result$p.value,
                    'treatment_estimate'=treatment_estimate,
                    'control_estimate'=control_estimate)
        )
    })


    t_test <- rbindlist(t_test)
    t_test <- cbind(dat[,c(1:3)], t_test)   # add back protein group / gene info cols
    t_test[, ratio := treatment_estimate / control_estimate]
    t_test[, p.adj := p.adjust(P, method='BH')]
    return(t_test[])
}


# run_contrast <- function(DT.all, treatment_sample_names, treatment, control)  {
#     control_sample_names <- colnames(DT.all)[colnames(DT.all) %like% control]
#     n_treatment <- length(treatment_sample_names)
#     n_controls <- length(control_sample_names)

#     info_cols <- colnames(DT.all)[1:3]
#     DT <- DT.all[, c(info_cols, treatment_sample_names, control_sample_names), with=F]

#     # Filter: protein groups with at least half of control and treatment samples having non-zero value
#     DT[, 'N_missing_treatment' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(treatment_sample_names)]
#     DT[, 'N_missing_control' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(control_sample_names)]
#     DT <- DT[N_missing_treatment < (n_treatment/2)]
#     DT <- DT[N_missing_control < (n_controls/2)]

#     DT[, 'control_var' := apply(.SD, 1, var), .SDcols=c(control_sample_names)]
#     DT[, 'treatment_var' := apply(.SD, 1, var), .SDcols=c(control_sample_names)]
#     DT[, 'overall_var' := apply(.SD, 1, var), .SDcols=c(control_sample_names, treatment_sample_names)]
#     DT <- DT[control_var != 0][treatment_var != 0][overall_var != 0]
#     DT.out <- copy(DT[, info_cols, with=F])
#     DT <- DT[, c(treatment_sample_names, control_sample_names), with=F]
#     # Run contrast
#     pvalue <- apply(DT, 1, function(x) {
#                                 a = factor(c(rep(treatment, n_treatment), rep(control, n_controls)),
#                                     levels=c(treatment, control))
#                                 fvalue = var.test(x~a)
#                                 if (!is.na(fvalue$p.value)){
#                                     if (fvalue$p.value > 0.05) {
#                                         t.test(x~a, var.equal = TRUE)
#                                     } else {
#                                         t.test(x~a, var.equal = FALSE)
#                                     }
#                                 }
#         }
#     )

#     DT.out[, p := as.numeric(unlist(lapply(pvalue,function(x) x$p.value)))]
#     DT.out[, q := p.adjust(p, method='BH')]       # BH method to convert P to FDR (q) value
#     DT.out[, ratio := as.numeric(unlist(lapply(pvalue,function(x) x$estimate[1]/(x$estimate[2])))) ]
#     return(DT.out[])
# }

# run_contrast_on_raw_count <- function(DT.all, treatment_sample_names, treatment, control)  {
#     control_sample_names <- colnames(DT.all)[colnames(DT.all) %like% control]
#     n_treatment <- length(treatment_sample_names)
#     n_controls <- length(control_sample_names)

#     info_cols <- colnames(DT.all)[1:3]
#     DT <- DT.all[, c(info_cols, treatment_sample_names, control_sample_names), with=F]
#     DT[, lapply(.SD, function(x) log_base^x), .SDcols=c(control_sample_names, treatment_sample_names)]

#     # Filter: protein groups with at least half of control and treatment samples having non-zero value
#     DT[, 'N_missing_treatment' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(treatment_sample_names)]
#     DT[, 'N_missing_control' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(control_sample_names)]
#     DT <- DT[N_missing_treatment < (n_treatment/2)]
#     DT <- DT[N_missing_control < (n_controls/2)]
#     DT.out <- copy(DT[, info_cols, with=F])
#     DT <- DT[, c(treatment_sample_names, control_sample_names), with=F]

#     # Run contrast
#     pvalue <- apply(DT, 1, function(x) {
#                                 a = factor(c(rep(treatment, n_treatment), rep(control, n_controls)),
#                                     levels=c(treatment, control))
#                                 fvalue = var.test(x~a)
#                                 if (!is.na(fvalue$p.value)){
#                                     if (fvalue$p.value > 0.05) {
#                                         t.test(x~a, var.equal = TRUE)
#                                     } else {
#                                         t.test(x~a, var.equal = FALSE)
#                                     }
#                                 }
#         }
#     )

#     DT.out[, p := as.numeric(unlist(lapply(pvalue,function(x) x$p.value)))]
#     DT.out[, q := p.adjust(p, method='BH', n=.N)]       # BH method to convert P to FDR (q) value
#     DT.out[, ratio := as.numeric(unlist(lapply(pvalue,function(x) x$estimate[1]/(x$estimate[2])))) ]
#     return(DT.out[])
# }


plot_volcano <- function(DT.original, treatment, control, log_base, lfc_threshold, fdr_threshold, out_dir, labelgene) {
    DT <- copy(DT.original)
    DT[, log_foldchange := log(ratio, base=log_base)]
    log_fc_threshold <- log(lfc_threshold, base=log_base)
    log_fdr_threshold <- -1*log10(fdr_threshold)
    DT[, labeltext := '']
    DT[, logq := -1*log10(p.adj) ]
    DT[, 'signif' := FALSE]
    DT[logq >= log_fdr_threshold, 'signif' := TRUE]
    DT[, 'magnitude' := 'LOW']
    DT[abs(log_foldchange) >= log_fc_threshold, 'magnitude' := 'HIGH']
    DT[, candidate := FALSE]
    DT[magnitude == 'HIGH' & signif == TRUE, candidate := TRUE]
    DT[candidate == TRUE, labeltext := Protein_Group]
    if (!is.null(labelgene)) {
        DT[Genes %like% labelgene, labeltext := Genes]
    }


    g <- ggplot(DT, aes(x=log_foldchange, y=logq, color=candidate)) +
        geom_point() +
        theme_few() +
        geom_text_repel(aes(label=labeltext), max.overlaps=50) +
        geom_hline(yintercept=log_fdr_threshold, linetype='dashed', alpha=0.5) +
        geom_vline(xintercept=log_fc_threshold, linetype='dashed', alpha=0.5) +
        geom_vline(xintercept=(-1*log_fc_threshold), linetype='dashed', alpha=0.5) +
        labs(x=paste0('Log[', opt$log_base, '](Intensity) fold change'),
                y='-Log[10](q)',
                title=paste0(treatment, ' vs ', control)
        ) +
        scale_color_manual(values=c('black','red')) +
        theme(legend.position = "none")
    output_filename <- paste0(treatment, '-vs-', control, '.png')
    cat(paste0('   -> ', out_dir, output_filename, '\n'))
    ggsave(g, filename=paste0(out_dir, output_filename), height=25, width=25, units='cm')
}

plot_volcano_from_raw <- function(DT.original, treatment, control, log_base, lfc_threshold, fdr_threshold, out_dir) {
    DT <- copy(DT.original)
    DT[, log_foldchange := log(ratio, base=log_base)]
    log_lfc_threshold <- log(lfc_threshold, base=log_base)
    log_fdr_threshold <- -1*log(fdr_threshold, base=log_base)
    DT[, labeltext := '']
    g <- ggplot(DT, aes(x=log_foldchange, y=-1*log10(q))) +
        geom_point() +
        theme_few() +
        geom_label_repel(aes(label=labeltext)) +
        geom_hline(yintercept=log_fdr_threshold, linetype='dashed', alpha=0.5) +
        geom_vline(xintercept=log_lfc_threshold, linetype='dashed', alpha=0.5) +
        geom_vline(xintercept=(-1*log_lfc_threshold), linetype='dashed', alpha=0.5) +
        labs(x=paste0('Log[', opt$log_base, '](Intensity) fold change'),
                y='-Log[10](q)',
                title=paste0(treatment, ' vs ', control)
        )
    output_filename <- paste0(treatment, '-vs-', control, '.png')
    cat(paste0('   -> ', out_dir, output_filename, '\n'))
    ggsave(g, filename=paste0(out_dir, output_filename), height=16, width=16, units='cm')
}
