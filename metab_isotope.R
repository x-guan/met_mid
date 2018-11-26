### This R code calculates deterium-labeled metabolites
### Written by Xin Guan (github.com/x-guan)

library(xcms)
library(ecipex)
library(tidyverse)

rm(list=ls())
setwd("/Users/Guan/Xin/R/Github/metab_isotope")

#########################################################################################
### bia intermediates of interest
# if xcms is used
my_mz_name <- c("norcoclaurine", "coclaurine", "methylcoclaurine", "hydroxy-methylcoclaurine", "reticuline", 
                "acutumine", "acutumidine", "mz_382.14", "mz_440.15",
                "scoulerine?", "tetrahydrocolumbamine?", "canadine?", "stylopine?")   # metabolite names
my_mz_value <- c(272.12, 286.14, 300.16, 316.15, 330.16, 
                 398.14, 384.12, 382.14, 440.15,
                 328.15, 342.17, 340.15, 324.12)   # mz values
my_retime <- c(299, 622, 347, 746, 353,
               316, 226, 540, 846,
               680, 722, 867, 1020)   # retention time in sec
my_formulae <- c("C16H17NO3", "C17H19NO3", "N18H21NO3", "C18H21NO4", "C19H23NO4",
                 "C19H24NO6Cl", "C18H22NO6Cl", "C18H20NO6Cl", "C22H30NO6Cl",
                 "C19H21NO4", "C20H23NO4", "C20H21NO4", "C19H17NO4")   # formulae
my_deut <- c(13, 16, 19, 18, 21,
             23, 20, 19, 30,
             20, 22, 21, 17)   # number of non-exchangeable deuteriums
# if python script is used
my_mz_name <- c("norcoclaurine", "coclaurine", "methylcoclaurine", "reticuline", 
                "acutumine", "acutumidine")   # metabolite names
my_formulae <- c("C16H17NO3", "C17H19NO3", "N18H21NO3", "C19H23NO4",
                 "C19H24NO6Cl", "C18H22NO6Cl")   # formulae
my_deut <- c(13, 16, 19, 21,
             23, 20)   # number of non-exchangeable deuteriums


########################################################################################
### xcms analysis of mzdata files
# mzdata files
file_dir <- "/Users/Guan/Xin/R/Github/metab_isotope/6520_mzdata/180302_mca_d2o"
mzdt <- list.files(path = file_dir, recursive = F, full.names = T)
mzdt <- mzdt[1:16]
mzdt_name <- gsub(paste0(file_dir, "/"), "", gsub(".mzdata.xml", "", mzdt))
# xcms pick
for (i in 1:length(mzdt)){
  mzdt_tmp <- mzdt[i]
  xset_1 <- xcmsSet(mzdt_tmp, fwhm = 20)   # pick peaks
  xset_2 <- as.data.frame(peaks(xset_1))
  write.csv(xset_2, file = paste0("180302_pick/", mzdt_name[i], ".csv"), row.names = F)
}
# xcms align
xset <- xcmsSet(mzdt, fwhm = 20)
xsg <- group(xset)
xsg <- retcor(xsg)
xsg <- group(xsg)
xsg <- fillPeaks(xsg)
dat <- groupval(xsg, "medret", "into")
# clean
dat_2 <- as.data.frame(dat)
colnames(dat_2) <- gsub(".mzdata", "", colnames(dat_2))
mz_rt <- as.data.frame(str_split(row.names(dat_2), "/", simplify = T))
dat_3 <- data.frame(mz = mz_rt[,1], rt = mz_rt[,2], dat_2) %>%
  rename(mca_12_l_d2o = mca_12_l1) %>%
  rename(mca_12_l_h2o = mca_12_l2) %>%
  rename(mca_12_p_d2o = mca_12_p1) %>%
  rename(mca_12_p_h2o = mca_12_p2)
write.csv(dat_3, file = "mca_d2o_align_1.csv", row.names = F)
# add info of feeding exp
dat_4 <- read.csv("mca_d2o_align_1.csv", header = T) %>%
  gather(sample_name, into, -c(mz, rt))
sample_name <- as.data.frame(str_split(dat_4$sample_name, "_", simplify = T))
dat_4$plant_id <- sample_name[,2]
dat_4$day <- with(dat_4, ifelse(plant_id == 1, 5, ifelse(plant_id == 3, 10, ifelse(plant_id == 13, 15, 0))))
dat_4$organ <- sample_name[,3]
dat_4$feed <- sample_name[,4]
write.csv(dat_4, file = "mca_d2o_align_2.csv", row.names = F)
# normalize by dry weight (per mg dry weight in 10ul 80% MeOH)
weight_info <- read.csv("mca_d2o_weight.csv", header = T)
dat_5 <- read.csv("mca_d2o_align_2.csv", header = T) %>%
  left_join(weight_info, by = c("plant_id", "organ", "feed")) %>%
  mutate(into_norm = into * (volume_ml*1000) / ((tube_dry_g - tube_g)*1000*10)) %>%
  group_by(mz, rt, sample_name, day, organ, feed) %>%
  summarise(into_norm = sum(into_norm))   # sum if identical mz/rt/id_rep occurs
write.csv(dat_5, file = "mca_d2o_align_3_norm.csv", row.names = F)
# extract datasets of bia intermediates
sub_1 <- read.csv("mca_d2o_align_3_norm.csv", header = T)
my_mz_name
my_mz_value
my_retime
isotope_info <- data.frame(mz_name = my_mz_name, mz_value = my_mz_value, retime = my_retime) %>%
  mutate(M0 = mz_value + 0, M1 = mz_value + 1, M2 = mz_value + 2) %>%
  mutate(M3 = mz_value + 3, M4 = mz_value + 4, M5 = mz_value + 5) %>%
  gather(isotope_name, isotope_mz, -c(mz_name, mz_value, retime)) %>%
  mutate(mz_upper = isotope_mz + 0.1, mz_lower = isotope_mz - 0.1) %>%
  mutate(rt_upper = retime + 12, rt_lower = retime - 12)
sub_2 <- NULL
for (i in 1:nrow(isotope_info)){
  dat_tmp1 <- with(sub_1, sub_1[mz < isotope_info$mz_upper[i] & mz > isotope_info$mz_lower[i],])
  dat_tmp2 <- with(dat_tmp1, dat_tmp1[rt < isotope_info$rt_upper[i] & rt > isotope_info$rt_lower[i],])
  if (nrow(dat_tmp2)){
    dat_tmp2$metabolite <- isotope_info$mz_name[i]
    dat_tmp2$ion <- isotope_info$isotope_name[i]
    sub_2 <- rbind.data.frame(sub_2, dat_tmp2)
  }
}
write.csv(sub_2, file = "mca_d2o_align_4_bia.csv", row.names = F)

#############################################
### fit isotopic abundance to mid
# adapted from Curt Fischer's code

### calculate theoretical MID
my_formulae
mids <- ecipex(my_formulae, isoinfo = nistiso, id = T, limit = 1e-6)
for(df_name in names(mids)){
  mids[[df_name]]$formula <- df_name
}
all_mids <- bind_rows(mids)

### calculate experimental MID
# if xcms is used
mid_1 <- read.csv("mca_d2o_align_4_bia.csv", header = T, stringsAsFactors = F)

### plot EIC abundance
mid_1 <- mid_1 %>%
  filter(organ == "l") %>%
  filter(day %in% c(5,10,15)) %>%
  filter(metabolite %in% c("norcoclaurine", "coclaurine", "methylcoclaurine", 
                           "reticuline", "acutumine", "acutumidine")) %>%
  mutate(metabolite = factor(metabolite, levels = c("norcoclaurine", "coclaurine", "methylcoclaurine", 
                                                    "reticuline", "acutumine", "acutumidine")))
write.csv(mid_1, "plot_eic.csv", row.names=F)
 
### plot MID and NE in a loop
my_mz_name
plot_mid <- NULL
plot_delta_mid <- NULL
plot_ne <- NULL
for (i in 1:length(my_mz_name)){
  ### if xcms is used
  mid_2 <- mid_1 %>%
    filter(metabolite == my_mz_name[i]) %>%
    filter(organ == "l") %>%
    group_by(day, organ, feed) %>%
    mutate(mid = into_norm / sum(into_norm)) %>%
    select(day, organ, feed, ion, mid) %>%
    arrange(day, organ, feed, ion)

  ### plot MID
  ggplot(mid_2, aes(x = ion, y = mid, color = factor(feed), ymin = 0, ymax = mid)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(width = 0.1, position=position_dodge(width=0.2)) +
    theme_bw() +
    ggtitle(paste0("MID (", my_mz_name[i], ")")) +
    ylab('MID') +
    xlab('Isotopes') +
    # facet_grid(day ~ organ, scale='free_y') +
    facet_grid(day ~ ., scale='free_y') +
    ggsave(paste0(my_mz_name[i], "_f1.pdf"), height = 5, width = 4)
  
  ### plot 'Experimental MID - theoretical MID'
  all_mid_annot <- all_mids %>%
    filter(formula == my_formulae[i]) %>%
    mutate(int_mass = round(mass, 0)) %>%
    group_by(formula, int_mass) %>%
    summarize(total_abund = sum(abundance), mass = sum(mass*abundance)/sum(abundance)) %>%
    group_by(int_mass) %>%
    arrange(int_mass) %>%
    mutate(ppm = (mass - min(mass))/(0.5*(mass + min(mass)))*1e6)
  # calculate 'Experimental MID - theoretical MID'
  with.theoretical <- all_mid_annot %>%
    ungroup %>%
    mutate(ion = paste0('M', (int_mass - min(int_mass)))) %>%
    filter(ion %in% c('M0', 'M1', 'M2', 'M3', 'M4', 'M5')) %>%
    mutate(theoretical_mid = total_abund / sum(total_abund)) %>%
    select(ion, theoretical_mid) %>%
    inner_join(mid_2, by = c('ion')) %>%
    rename(experimental_mid = mid) %>%
    group_by(ion, day, organ) %>%
    mutate(diff = experimental_mid - theoretical_mid)
  plot_delta_mid <- rbind.data.frame(plot_delta_mid, with.theoretical %>% mutate(metabolite = my_mz_name[i]))
  # plot
  with.theoretical %>%
    # filter(organ == "l") %>%
    ggplot(aes(x = day %>% as.numeric, 
               y = diff, color = factor(feed), ymin=0, ymax=diff)) +
    geom_point(position = position_dodge(width=0.2)) +
    # facet_grid(ion ~ organ, scale='free_y') +
    facet_grid(ion ~ ., scale='free_y') +
    theme_bw() +
    ggtitle(paste0("MID ~ days (", my_mz_name[i], ")")) +
    # ylab('Experimental MID - theoretical MID') +
    scale_y_continuous('Experimental MID - theoretical MID', 
                       limits = c(min(with.theoretical$diff), max(with.theoretical$diff))) +
    stat_smooth(method=lm, formula=y~poly(x,2), se=F, size=0.2) +
    xlab('Timepoint (days)') +
    ggsave(paste0(my_mz_name[i], "_f2.pdf"), height = 5, width = 4)
  # plot
  with.theoretical %>%
    filter(organ == "l") %>%
    ggplot(aes(x = ion, y = diff, color = factor(feed), ymin=0, ymax=diff)) +
    geom_errorbar(width=0.1, position=position_dodge(width=0.2)) +
    # facet_grid(day ~ organ, scale='free_y') +
    facet_grid(day ~ ., scale='free_y') +
    theme_bw() +
    ggtitle(paste0("MID ~ isotopes (", my_mz_name[i], ")")) +
    # ylab('Experimental MID - theoretical MID') +
    scale_y_continuous('Experimental MID - theoretical MID', 
                       limits = c(min(with.theoretical$diff), max(with.theoretical$diff))) +
    stat_smooth(method=lm, formula=y~poly(x,2), se=F, size=0.2) +
    xlab('Isotopes') +
    ggsave(paste0(my_mz_name[i], "_f3.pdf"), height = 5, width = 4)
  
  ### plot neutron excess
  with.theoretical_2 <- with.theoretical %>%
    filter(organ == "l") %>%
    ungroup %>%
    group_by(day, organ, feed) %>%
    mutate(num_ion = str_sub(ion, 2, 2) %>% as.integer) %>%
    mutate(diff_ne = diff * num_ion) %>%
    summarize(neutron_excess = sum(diff_ne)/my_deut[i])
  plot_ne <- rbind.data.frame(plot_ne, with.theoretical_2 %>% mutate(metabolite = my_mz_name[i]))
  # plot
  ggplot(with.theoretical_2, aes(x = day, y = neutron_excess, color = factor(feed))) +
    geom_point() +
    geom_line() +
    theme_bw() +
    ggtitle(paste0("Neutron excess (", my_mz_name[i], ")")) +
    # facet_grid(organ ~ ., scale='free_y') +
    ylab('Experimental_NE - natural_NE') +
    xlab('Timepoint (days)') +
    ggsave(paste0(my_mz_name[i], "_f4.pdf"), height = 2.5, width = 4)
}
write.csv(plot_mid, file="plot_mid.csv", row.names=F)
write.csv(plot_delta_mid, file="plot_mid_delta.csv", row.names=F)
write.csv(plot_ne, file="plot_ne.csv", row.names=F)

#################################################################################################
### plot eic
dat_eic <- read.csv("plot_eic.csv", header=T) %>%
mutate(metabolite = factor(metabolite, levels = c("norcoclaurine", "coclaurine", "methylcoclaurine", 
                                                  "reticuline", "acutumine", "acutumidine"))) %>%
  mutate(feed = factor(feed, levels = c("h2o", "d2o")))
ggplot(dat_eic, aes(x = ion, y = into_norm, color = factor(feed))) +
  geom_bar(aes(fill = feed), color = "black", width = 0.7, size = 0.5, stat = "identity", position = "dodge") +
  theme_bw() +
  # theme(axis.text.x=element_text(angle = 35,hjust = 1,vjust = 1.0)) +
  scale_x_discrete("Days of feeding (5, 10 and 15) and isotopes (M0 - M6)") +
  scale_y_continuous("EIC abundance") +
  facet_grid(metabolite ~ day %>% as.numeric, scale="free_y") +
  ggtitle("EIC of BIA intermediates (leaves)") +
  ggsave("plot_eic.pdf", height = 7.5, width = 8)

dat_eic_2 <- read.csv("plot_eic.csv", header=T) %>%
  group_by(day, feed, metabolite) %>%
  summarize(sum = sum(into_norm)) %>%
  ungroup %>%
  group_by(metabolite) %>%
  summarize(mean = mean(sum), sd = sd(sum)) %>%
  mutate(metabolite = factor(metabolite, levels = c("norcoclaurine", "coclaurine", "methylcoclaurine", 
                                                    "reticuline", "acutumine", "acutumidine")))
write.csv(dat_eic_2, "plot_eic_2.csv", row.names=F)
ggplot(dat_eic_2, aes(x=metabolite, y=mean)) +
  geom_bar(width=0.5, stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0) +
  theme_bw() +
  scale_y_continuous("EIC") +
  scale_x_discrete("BIA pathway intermediates") +
  ggsave("plot_eic_2.pdf", height = 2, width = 6)

ggplot(dat_eic, aes(x = ion, y = into_norm, color = factor(feed))) +
  geom_bar(aes(fill = feed), color = "black", width = 0.7, size = 0.5, stat = "identity", position = "dodge") +
  theme_bw() +
  # theme(axis.text.x=element_text(angle = 35,hjust = 1,vjust = 1.0)) +
  scale_x_discrete("Days of feeding (5, 10 and 15) and isotopes (M0 - M6)") +
  scale_y_continuous("EIC abundance") +
  facet_grid(metabolite ~ day %>% as.numeric, scale="free_y") +
  ggtitle("EIC of BIA intermediates (leaves)") +
  ggsave("plot_eic.pdf", height = 7.5, width = 8)

# ### plot mid
# dat_mid <- read.csv("plot_mid.csv", header=T) %>%
#   mutate(metabolite = factor(metabolite, levels = c("norcoclaurine", "coclaurine", "methylcoclaurine", 
#                                                     "reticuline", "acutumine", "acutumidine"))) %>%
#   mutate(feed = factor(feed, levels = c("h2o", "d2o")))
# ggplot(dat_mid, aes(x = ion, y = mid, color = factor(feed), ymin = 0, ymax = mid)) +
#   geom_point(position = position_dodge(width = 0.2)) +
#   geom_errorbar(width = 0.1, position=position_dodge(width=0.2)) +
#   theme_bw() +
#   scale_x_discrete("Days of feeding (5, 10 and 15) and isotopes (M0 - M6)") +
#   scale_y_continuous("MID(experimental)") +
#   facet_grid(metabolite ~ day %>% as.numeric, scale="free_y") +
#   ggtitle("MID of BIA intermediates (leaves)") +
#   ggsave("plot_mid.pdf", height = 7.5, width = 8)

### plot mid_delta
dat_mid_delta <- read.csv("plot_mid_delta.csv", header=T) %>%
  mutate(metabolite = factor(metabolite, levels = c("norcoclaurine", "coclaurine", "methylcoclaurine", 
                                                    "reticuline", "acutumine", "acutumidine"))) %>%
  mutate(feed = factor(feed, levels = c("h2o", "d2o")))
ggplot(dat_mid_delta, aes(x=ion, y = diff, color = factor(feed), ymin = 0, ymax = diff)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(width = 0.1, position=position_dodge(width=0.2)) +
  theme_bw() +
  scale_x_discrete("Days of feeding (5, 10 and 15) and isotopes (M0 - M6)") +
  scale_y_continuous("MID(experimental) - MID(theoretical)") +
  facet_grid(metabolite ~ day, scale="free_y") +
  ggtitle("MID of BIA intermediates (leaves)") +
  ggsave("plot_mid_delta.pdf", height = 7.5, width = 8)

### plot ne
dat_ne <- read.csv("plot_ne.csv", header=T)
# plot_ne_1
dat_ne_1 <- dat_ne %>%
  mutate(metabolite = factor(metabolite, levels=c("norcoclaurine", "coclaurine", "methylcoclaurine", 
                                                  "reticuline", "acutumine", "acutumidine"))) %>%
  mutate(feed = factor(feed, levels=c("h2o", "d2o")))
ggplot(dat_ne_1, aes(x=day, y=neutron_excess, group=feed, fill=feed)) +
  geom_bar(width=4, stat="identity", position="dodge") +
  # geom_boxplot() +
  # geom_dotplot(binaxis="y", stackdir="center", dotsize=1) +
  # geom_jitter(aes(x=metabolite, y=neutron_excess)) +
  theme_bw() +
  facet_wrap(~ metabolite, scale="free", ncol=3, drop=T) +
  scale_x_continuous("BIA pathway intermediates and days of feeding (5, 10 and 15)", breaks=c(5,10,15)) +
  scale_y_continuous("Neutron excess", limits=c(-0.02, 0.08), breaks=seq(-0.02, 0.08, by=0.02)) +
  ggsave("plot_ne.pdf", height = 4, width = 6)
# plot_ne_2
dat_ne_2 <- dat_ne %>%
  spread(feed, neutron_excess) %>%
  mutate(diff_ne = d2o-h2o) %>%
  group_by(metabolite) %>%
  summarise(mean = mean(diff_ne), sd = sd(diff_ne)) %>%
  mutate(metabolite = factor(metabolite, levels=c("norcoclaurine", "coclaurine", "methylcoclaurine", 
                                                  "reticuline", "acutumine", "acutumidine")))
write.csv(dat_ne_2, "plot_ne_2.csv", row.names=F)
ggplot(dat_ne_2, aes(x=metabolite, y=mean)) +
  geom_bar(width=0.5, stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0) +
  theme_bw() +
  scale_y_continuous("NE_D2O - NE_H2O") +
  scale_x_discrete("BIA pathway intermediates") +
  ggsave("plot_ne_2.pdf", height = 2, width = 6)
# calculate p-value from paired 2-sample t-test
ttest_ext <- NULL
for (i in 1:length(unique(dat_ne_2$metabolite))){
  dat_ne_3 <- dat_ne_1 %>%
    filter(metabolite == unique(dat_ne_2$metabolite)[[i]]) %>%
    spread(feed, neutron_excess)
  ttest_re <- t.test(dat_ne_3$d2o, dat_ne_3$h2o, paired=T, alternative="greater")
  ttest_ext_tmp <- data.frame(metabolite=unique(dat_ne_2$metabolite)[[i]], 
                              mean_diff=ttest_re[[5]], pvalue=ttest_re[[3]])
  ttest_ext <- rbind.data.frame(ttest_ext, ttest_ext_tmp)
}
write.csv(ttest_ext, "plot_ttest.csv", row.names=F)






