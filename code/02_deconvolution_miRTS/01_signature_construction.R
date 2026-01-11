# Create a signature matrix that include
# marker miRNAs based on comparison with the avg of all other tissues,
# not any individual tissue (not top1 or top3)

library(edgeR)
df_input = "TA2_miR.raw_qc" #TA2_miR.3_new2024.sel  TA2_miR.qc_cpm

temp <- as.data.frame(t(as.data.frame( get(df_input) )))
temp <- temp[! gsub("__.*","",rownames(temp)) %in% c(
  #1
  "diaphragm", "glandular_breast_tissue",
  "prostate", "sclera", "urethra", "uterus"
  # 2
  , "adrenal_gland", "gallbladder", "spinal_cord", "tongue",
  "trachea"
),] # select the samples with n>2
rownames(temp) <- gsub("submandibular_gland", "salivary_gland", rownames(temp))
temp_input <- temp

sampleID_organ <- data.frame(
  sampleID=rownames(temp_input),
  organ=gsub("__.*","",rownames(temp_input))
)
sort(table(sampleID_organ$organ))
sampleID_organ.0 <- sampleID_organ

tissue_using = "bladder"

tissue_list <- c(
  "bladder", "lymph_node", "muscle", "salivary_gland", "pleurae", "spleen",
  "testis", "thyroid", "vein", "adipocyte", "artery", "bone",
  "esophagus", "pancreas", "skin", "stomach", "kidney",
  "liver", "nerve", "lung", "heart", "bowel", "brain"
)

marker_miR.compOthers <-
  data.frame()

for (tissue_using in tissue_list){

  sampleID_organ <- sampleID_organ.0 %>%
    mutate(organ=
             ifelse(organ==tissue_using, tissue_using, "Others"))

  d0 <- DGEList(t(temp_input))
  d0 <- calcNormFactors(d0)
  # d0$samples %>% arrange(norm.factors)

  mm <- model.matrix(~0 + organ, data = sampleID_organ)
  colnames(mm) <- gsub("^organ", "", colnames(mm))
  mm

  # keep <- filterByExpr(d0, mm )
  # sum(keep) # number of genes retained
  # d <- d0[keep,]
  d <- d0 #[keep,]   # Not using filterByExpr, since some miRs may be uniquely expressed in a tissue type.
  # plotMDS(d,labels = factor(gsub("__.*","", sampleID_organ$organ)), col = as.numeric(factor(sampleID_organ$organ)), cex = 1)

  y <- voom(d, mm, plot = T)
  fit <- lmFit(y, mm)
  head(coef(fit))

  make_all_contrasts <- function (group, delim="_vs_"){
    group <- sort(unique(as.character(group)))
    cb   <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])})
    contrasts<- limma::makeContrasts(contrasts=cb, levels=group)
    colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
    return(contrasts)
  }

  contr <- make_all_contrasts(group = unique(sampleID_organ$organ))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  contr_res.coef <- as.data.frame(tmp$coefficients) %>%
    rownames_to_column("miRNA") %>%
    melt(., id.var="miRNA",) %>% rename(contrast=variable)
  contr_res.coef.DE2 <- contr_res.coef #%>% filter((value)>log2(DE))

  contr_res.p <- as.data.frame(tmp$p.value) %>%
    rownames_to_column("miRNA") %>%
    melt(., id.var="miRNA",) %>% rename(contrast=variable)
  contr_res.p.lt.05 <- contr_res.p #%>% filter(value<P_using)

  TA2.Signature_miR.DE2.05 <-
    inner_join(contr_res.coef.DE2, contr_res.p.lt.05, by=c("miRNA", "contrast"))
  TA2.Signature_miR.DE2.05 <- TA2.Signature_miR.DE2.05 %>%
    mutate(enriched_organ=gsub("_vs_Others|Others_vs_", "", contrast))
  print(table(TA2.Signature_miR.DE2.05$enriched_organ))
  xx <- TA2.Signature_miR.DE2.05 %>%
    mutate(DE=value.x, DE_p=value.y) %>%
    select(miRNA, enriched_organ, DE, DE_p)
  marker_miR.compOthers <-
    rbind(marker_miR.compOthers, xx)
}

marker_miR.compOthers.sel <-
  marker_miR.compOthers %>%
  group_by(enriched_organ) %>%
  mutate(FDR=p.adjust(DE_p, method="BH")) %>%
  filter(DE>2 , FDR<0.2) # , DE_p<0.05)
marker_miR.compOthers.sel <- marker_miR.compOthers.sel %>%
  filter(miRNA %in% gsub("_", "-", Plasma_miR.all.1_intersect.2))
dim(marker_miR.compOthers.sel);  length(unique(marker_miR.compOthers.sel$miRNA))
aa <- names(table(marker_miR.compOthers.sel$enriched_organ))
tissue_list[!tissue_list %in% aa]

# marker_miR.compOthers.DE2.FDR.05__Plsm_union
# marker_miR.compOthers.DE2.FDR.2__Plsm_union
# marker_miR.compOthers.DE2.P.05__Plsm_union
# marker_miR.compOthers.DE2.FDR.05__Plsm_intersect
# marker_miR.compOthers.DE2.FDR.2__Plsm_intersect
# marker_miR.compOthers.DE2.P.05__Plsm_intersect
# marker_miR.compOthers.None__Plsm_intersect

assign( "marker_miR.compOthers.DE2.FDR.2__Plsm_intersect.8",
        list(Signature=
               TA2_tissue.TMM[unique(marker_miR.compOthers.sel$miRNA),
                              unique(marker_miR.compOthers.sel$enriched_organ)],
             miR_organ=marker_miR.compOthers.sel) )
all.equal(marker_miR.compOthers.DE2.FDR.2__Plsm_intersect.2,
          list(Signature=
                 TA2_tissue.TMM[unique(marker_miR.compOthers.sel$miRNA),
                                unique(marker_miR.compOthers.sel$enriched_organ)],
               miR_organ=marker_miR.compOthers.sel) )
