data_dir<-"~/chant2lab/ChanLab_Common_Drive/Prerana/projects/Binding_AUC_analysis/data/"
bms153_peptides<-read.delim(paste0(data_dir,"bms153_netMHCout.txt"),sep="\t",header=TRUE)
#bms<-read.delim(paste0(data_dir,"netMHCOutput/bms153_NetMHCpan.txt"),header=TRUE,skip = 1,sep="\t")
bms153_peptides$binder2<-bms153_peptides$Binder
bms153_peptides$tetramer2<-bms153_peptides$tetramer
bms153_peptides$bin_label<-ifelse(bms153_peptides$binder2=="REAL_BINDER",1,0)
bms153_peptides$tet_label<-ifelse(bms153_peptides$tetramer2=="TETRAMER+",1,0)
bms153_peptides$tet_label[which(bms153_peptides$binder2 == "NON_BINDER")]<-0
bms_tet<-bms153_peptides[which(!is.na(bms153_peptides$tetramer2)),]

roc_obj_binding <- roc(response = bms153_peptides$bin_label, predictor = bms153_peptides$BA.score)
roc_obj_immuno_binders_only <- roc(response = bms_tet$tet_label, predictor = bms_tet$BA.score)
roc_obj_immuno <- roc(response = bms153_peptides$tet_label, predictor = bms153_peptides$BA.score)

roc_list <- list(
  Binding = roc_obj_binding,
  Immunogenicity = roc_obj_immuno,
  Immunogenicity_binders  = roc_obj_immuno_binders_only
)

auc_vals <- sapply(roc_list, auc)

roc<-ggroc(roc_list) +
  theme_classic() +
  labs(
    color = "Model",
    subtitle = paste(
      names(auc_vals),
      sprintf("(AUC = %.3f)", auc_vals),
      collapse = "   "
    )
  )

pdf("./BMS153/ROC_Fig5.pdf",height = 6,width=12)
print(binder_roc)
dev.off()