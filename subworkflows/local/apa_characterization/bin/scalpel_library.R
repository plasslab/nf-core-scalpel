# [Script for Function needed in Scalpel analysis]
# ------------------------------------------------
# Differential isoform usage functions, verbatim from plasslab/SCALPEL src/scalpel_library.R
# (plotting functions & their Gviz/GenomicRanges requirements not needed here)
require(stringr)
require(stringr)


FindIsoforms <- function (seurat_obj, group_by, split_by = NULL, assay = "RNA",
                           threshold_pval = 0.05, threshold_abund = 0.1, threshold_var=0.05,
                           extended=F) {
    
    #check presence of group.by into seurat object
    if (!(group_by %in% colnames(seurat_obj[[]]))) {
        
        stop(paste0("Error: ", group.by, "not found in Seurat object !"))
        
    } else {
        #process seurat object if split.by is NULL
        if ( is.null(split_by) ) {
            message("process seurat object...")
            
            #get tab.counts (1)
            message("Aggregating counts in conditions...")
            counts.tab = Extract_aggTab(seurat_obj, group = group_by, assay)
            
            #process filtering ops by gene (2)
            message("Performing filtering and comparison test...")
            all_results = pbapply::pblapply(split(counts.tab, counts.tab$gene_name), function(geneTab) {
                #apply filtering checks & X2 test...
                if(nrow(geneTab)>1){
                    return(Perform_xTest(geneTab, threshold_abund, threshold_var))
                } else {
                    tmp = geneTab %>% dplyr::mutate(thr.abund=T, thr.var=T, p.value=NA, statistic=NA)
                    return(tmp)
                }
            }) %>% rbindlist()
            
            #calculate adjusted P.values (3)
            all_results = all_results %>%
                dplyr::mutate(p.value.adjusted = ifelse(thr.abund==F & thr.var==F & (!(is.na(p.value))), 0, NA))
            
            tmpA = dplyr::filter(all_results, p.value.adjusted==0)
            tmpA$p.value.adjusted = p.adjust(tmpA$p.value, method = "BH")
            tmpB = dplyr::filter(all_results, is.na(p.value.adjusted))
            all_results = rbind(tmpA, tmpB) %>% dplyr::arrange(p.value.adjusted, gene_name, transcript_name) %>%
                dplyr::distinct() %>%
                data.table::data.table()
            
            #apply p.value threshold (4)
            if(extended==F){
                all_results = dplyr::filter(all_results, p.value.adjusted<threshold_pval) %>%
                    select(!c("thr.abund", "thr.var"))
            } else {
                return(all_results)
            }
            
            #return
            return(all_results)
            
        } else {
            #check presence of split.by var into seurat object
            if (!(split_by %in% colnames(seurat_obj[[]]))) {
                
                stop(paste0("Error: ", group_by, "not found in Seurat object !"))
                
            } else {
                message("split seurat object and process with each object...")
                
                #splitting
                seurat.objs = Seurat::SplitObject(JoinLayers(seurat_obj), split.by = split_by)
                
                #processing
                message("Generation of aggregate counts...")
                lapply(names(seurat.objs), function(curr_obj){
                    message(curr_obj)
                    
                    #get tab.counts
                    counts.tab = Extract_aggTab(seurat.objs[[curr_obj]], group = group_by, assay, n=1) %>% 
                        dplyr::mutate(split.by = curr_obj)
                    
                })
            }
        }
    }
}

# Functions associateds..
#========================

Extract_aggTab <- function (seurat_obj, group, assay) {
    #Get aggregated counts tab
    tmp = Seurat::AggregateExpression(seurat_obj, group.by = group, assays = assay)[[1]]
    names.cols = colnames(tmp)
    tmp = tmp %>%
        tibble::as_tibble(rownames="gene.tr") %>%
        tidyr::separate(gene.tr, c("gene_name","transcript_name"), sep="\\*\\*\\*", remove = F) %>%
        dplyr::arrange(gene_name, transcript_name)
    
    #get genes with nb isoforms > n
    tmp.stats = dplyr::distinct(tmp, gene_name, transcript_name) %>%
        dplyr::group_by(gene_name) %>%
        dplyr::summarise(nb.isoforms = n_distinct(transcript_name))
    
    tmp = data.table::data.table(dplyr::filter(tmp, gene_name %in% tmp.stats$gene_name))
    colnames(tmp) = c("gene.tr", "gene_name", "transcript_name", names.cols)
    tmp = tmp %>% mutate(across(where(is.numeric), ~tidyr::replace_na(., 0)))
    return(tmp)
}

Perform_xTest <- function(tab, threshold.abund, threshold.var) {
    #process filtering ops
    tab = tab %>%
        rowwise() %>%
        mutate(
            row_sum = sum(c_across(where(is.numeric)))
        ) %>%
        ungroup() %>%
        mutate(
            thr.abund = if_all(where(is.numeric) & !matches("row_sum"), ~ (. / sum(.)) < threshold.abund),
            thr.var = all(across(where(is.numeric) & !matches("row_sum"), ~ sd(.) < 0.05))
        ) %>%
        select(-row_sum)
    
    #get numeric table
    numTab = tab %>% dplyr::filter(thr.abund==F & thr.var==F) %>%
        dplyr::select(where(is.numeric))
    
    #perform X2 test
    if(nrow(numTab)>1){
        #perform test if more than 1 isoforms following filtering ops
        chi2_res <- chisq.test(as.matrix(numTab)) %>% suppressWarnings()
        
        #check if X2 test condition are reunited (counts>=5)...
        if (all(round(chi2_res$expected,0) >= 5)) {
            #update tab
            tab = tab %>% dplyr::mutate(
                p.value = chi2_res$p.value,
                statistic = chi2_res$statistic
            )
        }else{
            #update tab
            tab = tab %>% dplyr::mutate(
                p.value = NA,
                statistic = NA
            )
        }
    } else {
        #update tab
        tab = tab %>% dplyr::mutate(
            p.value = NA,
            statistic = NA
        )
    }
    #return
    return(tab)
}



