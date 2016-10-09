#Bacterial Pathogenesis Tutorial - Skin microbiome and infection data analysis
library("phyloseq")
library("plyr")
library("ggplot2")
library("reshape")
#Set working directory to where your files are - YOU WILL NEED TO CHANGE THIS LINE RACHEL
setwd("~/Documents/UTS/metaSUB/Honours_2016/data/run3/R_analysis/")
#import data to phyloseq
map<-import_qiime_sample_data("mapping_metasub.txt")
biom_otu_table<- import_biom("otu_table_uparse_tax_goodalign_json.biom", "rep_set_tree.tre")
trains_otu <- merge_phyloseq(biom_otu_table,map)
colnames(tax_table(trains_otu)) = c("Domain", "Phylum", "Class", "Order", 
                                 "Family", "Genus", "Species")
#View sample metadata
trains_sample_data<-as.matrix(sample_data(trains_otu))
#There should now be an object called trains sample map that 
#you can click on to easily view your sample data.
#Five number summary of sequence coverage of samples
fivenum(colSums(otu_table(trains_otu)))
#View number of sequences per sample
seq_sums<-sort(sample_sums(trains_otu),decreasing=TRUE)
seq_sums<-as.matrix(seq_sums)
#Will create an obsject showing sequence counts per sample in decreasing order, 
#you can click on this from the environment pane to see your samples and counts
#Rarefy to 5000 seqs per sample
rare_trains_5k<-rarefy_even_depth(trains_otu, sample.size=5000, replace=FALSE, 
                                trimOTUs=TRUE, rngseed=711)
#View first 10 lines of the OTU table
otu_table(rare_trains_5k)[1:10,]
otus_sorted<-otu_table(rare_trains_5k)[order(rowSums(otu_table(rare_trains_5k)),
                              decreasing=T),]
otus_sorted[1:10,]
#Show overall proportion of top 20 OTUs to the total OTU table (rarefied)
barplot(sort(taxa_sums(rare_trains_5k),
             TRUE)[1:20]/sum(sample_sums(rare_trains_5k)),
        las=2, ylab="proportion")
#Alpha diversity estimates and plotting
phyloseq::plot_richness(rare_trains_5k, x="Surface", 
                        measures=c("Observed", "Chao1","Shannon"))+
  geom_boxplot()+xlab("Group")+ylab("Diversity")
#Test significance of alpha diversity
alpha_estimates<-estimate_richness(rare_trains_5k, split = TRUE, measures = c("Observed", "Chao1", "Shannon"))
#Sort samples in alpha_estimates by sample name
alpha_estimates_sorted<-alpha_estimates[order(rownames(alpha_estimates)),]
#sort samples in sample_data(rare_trains_5k) by sample name
sorted_sampledata<-sample_data(rare_trains_5k)[order(rownames(sample_data(rare_trains_5k))),]
#Create a vector of which group each sample belongs to from the 'Surface' column and add to alpha estimates file
group<-sorted_sampledata$Surface
alpha_estimates_sorted<-cbind(alpha_estimates_sorted, group)
rm(group)
#Sort sorted_alpha_estimates by group
alpha_estimates_sorted<-alpha_estimates_sorted[order(alpha_estimates_sorted$group),]

#Find test for among several groups, then do post-hoc testing in Wilcoxin rank sum test

#Beta_diversity
BC_ord<-ordinate(rare_trains_5k, method="NMDS", "wunifrac")
plot_ordination(rare_trains_5k, BC_ord, type = "samples", 
                shape="Surface", color="Station",label="X.SampleID",
                title="All Samples")
#N_1 is a clear outlier, will need to remove this from analysis and do the above again.

#Plotting data from OTU table
#Get data from OTU table in format for plotting
otu_rare_trains_otu<-t(otu_table(rare_trains_5k))
otu_rare_trains_otu<-(otu_rare_trains_otu/5000)*100
otu_rare_trains_otu<-otu_rare_trains_otu[,order(-colSums(otu_rare_trains_otu))]
surface<-sample_data(rare_trains_5k)[sample_names(rare_trains_5k),"Surface"]
otu_rare_trains_otu<-cbind(otu_rare_trains_otu, surface)
####UP to here
otu_rare_trains_otu_sort<-otu_rare_trains_otu[order(otu_rare_trains_otu$Surface),]
otu_rare_melt<-melt(otu_rare_trains_otu_sort, id="Surface")
otu_rare_melt_top10<-otu_rare_melt[1:540,]
#Plot stacked bar graph
ggplot(otu_rare_melt_top10, aes(x=variable, y=value, fill=Surface))+
  stat_summary(fun.y=mean, position="stack", geom="bar")+
  labs(y="Relative abundance %", x="OTU")+
  scale_y_continuous(breaks=seq(0,100,10))
#Plot side by side bar graph with standard error
ggplot(otu_rare_melt_top10, aes(x=variable, y=value, fill=Surface)) +stat_summary(fun.y=mean, geom="bar",position=position_dodge(1)) + 
  stat_summary(fun.ymin = function(x) mean(x)-sd(x), 
               fun.ymax = function(x) mean(x)+sd(x),
               geom="errorbar",position=position_dodge(1), width=.2)+
  labs(y="Relative abundance %", x="OTU")+
  scale_y_continuous(breaks=seq(0,100,10))
#Get identities of OTUs
pre_top10_taxa<-tax_table(rare_trains_5k)[colnames(otu_rare_trains_otu_sort)[1:10],]
pre_top10_taxa
#then show results of differential abundance testing in QIIME.

###FROM HERE DOWN I HAVEN'T ADJUSTED THE SCRIPT TO YOUR SAMPLES RACHEL!


#Differential abundance analysis with DESeq
#This chunk of code will not work on computer lab computers as DESeq2 could not be installed
#library("DESeq2")
#skin_pre_des<-subset_samples(skin_data, infection_outcome_s!="control")
#skin_pre_des<-subset_samples(skin_pre_des, day_of_study_s=="preinfection")
#skin_pre_des<-filter_taxa(skin_pre_des, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
#skin_pre_des<-phyloseq_to_deseq2(skin_pre, ~subject_health_status)
#skin_pre_des$subject_health_status <- relevel(skin_pre_des$subject_health_status, ref="resolver")
#skin_pre_des<-DESeq(skin_pre_des)
#plotDispEsts(skin_pre_des)
#res05<- results(skin_pre_des, alpha=0.05)
#resOrdered <- res05[order(res05$padj),]
#signum<-sum(res05$padj < 0.05, na.rm=TRUE)
#resOrdered[1:signum,]
#Look at abundance of sig diff OTUs between pustule formers and resolvers
#sig_otus<-rownames(resOrdered)[1:7]

#I ran the analysis on my machine and got the following results:
 # log2 fold change (MAP): subject_health_status pustule_former vs resolver 
#  Wald test p-value: subject_health_status pustule_former vs resolver 
# DataFrame with 7 rows and 6 columns
#            baseMean    log2FoldChange lfcSE      stat       pvalue        padj
#           <numeric>      <numeric> <numeric> <numeric>    <numeric>   <numeric>
#  OTU_7    165.7817390       7.307377  1.754827  4.164159 3.125022e-05 0.001656261
#  OTU_145    3.2377725       6.850261  1.894212  3.616416 2.987099e-04 0.007915813
#  OTU_1406  25.1749418       7.446709  2.231192  3.337548 8.452122e-04 0.014932083
#  OTU_14    44.0157143       3.769900  1.195741  3.152774 1.617269e-03 0.021428815
#  OTU_82     1.6727456       5.869333  2.008025  2.922938 3.467450e-03 0.030629142
#  OTU_96     2.7298277       6.215280  2.105857  2.951426 3.163105e-03 0.030629142
#  OTU_1233   0.9598385       5.300857  1.944421  2.726188 6.407051e-03 0.048510532


#Manually create a lit of the significantly different OTUs.
sig_otus<-c("OTU_7","OTU_145","OTU_1406","OTU_14","OTU_82","OTU_96","OTU_1233")
tax_table(skin_pre)[sig_otus,]
sig_otu_table<-data.frame(otu_pre_sort[,sig_otus])
sig_otu_table$subject_health_status<-otu_pre_sort$subject_health_status
sig_otu_table
sig_otu_table_melt<-melt(sig_otu_table, id="subject_health_status")
sig_dif_plot<-ggplot(sig_otu_table_melt, aes(x=variable, y=value, fill=subject_health_status)) +stat_summary(fun.y=mean, geom="bar",position=position_dodge(1)) + 
  stat_summary(fun.ymin = function(x) mean(x)-sd(x), 
               fun.ymax = function(x) mean(x)+sd(x),
               geom="errorbar",position=position_dodge(1), width=.2)+
  labs(y="Relative abundance %", x="OTU")+
  scale_y_continuous(breaks=seq(-10,100,10))
sig_dif_plot
sig_dif_plot+scale_y_log10()
tax_table(skin_pre)[sig_otus,]
