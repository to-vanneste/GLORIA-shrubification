setwd("~/GLORIA/Data/Paper 1")

library(dplyr)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(car)
library(scales)
library(rstan)
library(brms)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(tidyr)
library(MCMCglmm)
library(tidybayes)
library(bayesplot)
library(ggpubr)
library(MuMIn)
library(broom.mixed)
library(ggeffects)
library(patchwork)
library(gmodels)
library(phyr)
library(rtrees)

# absolute cover changes (log-response ratios)

df <- read.csv("shrubs.csv", header = T)
test3 <- df %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(RR = log(cover_tot/lag(cover_tot))/(year - lag(year)), check = cover_tot + lag(cover_tot), check2 = cover_tot - lag(cover_tot)) %>%
  filter(check != 0)

df <- read.csv("shrubs.csv", header = T)
test4 <- df %>%
  filter(DE == "d") %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(RR = log(cover_tot/lag(cover_tot))/(year - lag(year)), check = cover_tot + lag(cover_tot)) %>%
  filter(check != 0)

df <- read.csv("shrubs.csv", header = T)
test5 <- df %>%
  filter(DE == "e") %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(RR = log(cover_tot/lag(cover_tot))/(year - lag(year)), check = cover_tot + lag(cover_tot)) %>%
  filter(check != 0)

p1<-ggplot(test3,aes(x=RR)) + 
  geom_histogram(bins = 100, aes(fill = RR >= 0, alpha = RR >= 0)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() +
  scale_fill_manual(values = c("black", "black")) +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_x_continuous(limits = c(-0.5, 1)) +
  labs(x = "", y = "No. of plots") +
  ggtitle("All shrubs") +
  theme(legend.position = "none", plot.title = element_text(face="bold"))

p2<-ggplot(test4,aes(x=RR)) + 
  geom_histogram(bins = 100, aes(fill = RR >= 0, alpha = RR >= 0)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() +
  scale_fill_manual(values = c("#9999CC", "#9999CC")) +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_x_continuous(limits = c(-0.5, 1)) +
  labs(x = "Cover change (log ratio/year)", y = "") +
  ggtitle("Deciduous shrubs") +
  theme(legend.position = "none", plot.title = element_text(face="bold"), axis.title.y = element_blank())

p3<-ggplot(test5,aes(x=RR)) + 
  geom_histogram(bins = 100, aes(fill = RR >= 0, alpha = RR >= 0)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() +
  scale_fill_manual(values = c("#66CC99", "#66CC99")) +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_x_continuous(limits = c(-0.5, 1)) +
  labs(x = "", y = "") +
  ggtitle("Evergreen shrubs") +
  theme(legend.position = "none", plot.title = element_text(face="bold"), axis.title.y = element_blank())

p <-p1 + p2 + p3 + plot_layout(ncol = 3, axes = "collect") 
png(file="Fig1_alternative.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

m <- lmer(RR~1+(1|region/summit/aspect), test3, na.action = na.omit, REML = T)
summary(m)
confint(m)

m <- lmer(RR~1+(1|region/summit/aspect), test4, na.action = na.omit, REML = T)
summary(m)
confint(m)

m <- lmer(RR~1+(1|region/summit/aspect), test5, na.action = na.omit, REML = T)
summary(m)
confint(m)

df <- read.csv("shrubs.csv", header = T)
df <- df  %>% filter(species_name == "Alyssum_spinosum" |
                       species_name == "Arctostaphylos_uva-ursi" |
                       species_name == "Betula_nana" |
                       species_name == "Calluna_vulgaris" |
                       species_name == "Dryas_octopetala" |
                       species_name == "Empetrum_nigrum" |
                       species_name == "Erica_herbacea" |
                       species_name == "Helianthemum_oelandicum" |
                       species_name == "Juniperus_communis" |
                       species_name == "Loiseleuria_procumbens" |
                       species_name == "Rhododendron_myrtifolium" |
                       species_name == "Salix_herbacea" |
                       species_name == "Salix_reticulata" |
                       species_name == "Thymus_nervosus" |
                       species_name == "Thymus_praecox" |
                       species_name == "Vaccinium_myrtillus" |
                       species_name == "Vaccinium_uliginosum" |
                       species_name == "Vaccinium_vitis-idaea")
df <- df  %>% select(year, region, summit, aspect, plot, species_name, cover)
df <- spread(df, species_name, cover)
df <- df %>% replace(is.na(.), 0)
df <- gather(df, species_name, cover, 6:23)

test <- df %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  ungroup() %>%
  arrange(region, summit, aspect, plot, year)

test2 <- df %>%
  arrange(region, summit, aspect, plot, species_name, year) %>%
  filter(year == min(year) | year == max(year)) %>%
  group_by(region, summit, aspect, plot, species_name) %>%
  mutate(RR = log((cover+0.001)/lag((cover+0.001)))/(year - lag(year)), check = cover + lag(cover)) %>%
  filter(check != 0)

comm <- test2 %>% select(region,summit,aspect,plot,species_name,RR)
colnames(comm) <- c("region","summit","aspect","plot","species","RR")
traits <- read.csv("shrub_traits2.csv", header = T)
df <- left_join(comm, traits, by = "species")

p<-ggplot(df, aes(reorder(POWO_accepted2, -RR, FUN=median), RR, color = DE)) +
  geom_boxplot() +
  geom_point(alpha=0.3) +
  #stat_summary(geom = "point",fun.y = "mean",col = "black",size = 3,shape = 24,fill = "black") +
  scale_color_manual(values = c("#9999CC", "#66CC99")) +
  theme_bw() +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.title = element_blank(), legend.box.margin = margin(0, 0, 0, 0), legend.text=element_text(size=7),
        legend.background = element_rect(fill="white",size=0.5,linetype="solid",colour ="gray"),
        legend.key.size = unit(0.2, "cm"), legend.spacing.y = unit(0, "pt"), legend.spacing.x = unit(0, "pt"), 
        axis.text.x = element_text(angle=50, vjust=1, hjust=1, face = "italic")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x="", y="Cover change (log-ratio/year)", color="")

png(file="Fig2_alternative.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

results <- df %>%
  group_by(POWO_accepted2, DE) %>%
  summarise(
    mean_RR = mean(RR),
    ci_lower = ci(RR)[2],
    ci_upper = ci(RR)[3])

p<-ggplot(results, aes(reorder(POWO_accepted2, -mean_RR, FUN=mean), mean_RR, color = DE)) +
  geom_point() +
  geom_errorbar(aes(x =reorder(POWO_accepted2, -mean_RR, FUN=mean), ymin = ci_lower, ymax = ci_upper), width = 0) +
  scale_color_manual(values = c("#9999CC", "#66CC99")) +
  theme_bw() +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.title = element_blank(), legend.box.margin = margin(0, 0, 0, 0), legend.text=element_text(size=7),
        legend.background = element_rect(fill="white",size=0.5,linetype="solid",colour ="gray"),
        legend.key.size = unit(0.2, "cm"), legend.spacing.y = unit(0, "pt"), legend.spacing.x = unit(0, "pt"), 
        axis.text.x = element_text(angle=50, vjust=1, hjust=1, face = "italic")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x="", y="Cover change (log-ratio/year)", color="")

png(file="Fig2_alternative2.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

#df <- left_join(results, traits, by = "POWO_accepted2")

m <- lmer(RR ~ scale(log(mean_Height)) + scale(mean_SLA) + scale(mean_SeedMass) + scale(log(mean_LeafN)) +
            (1|species) + (1|region/summit/aspect), data = df, na.action = na.omit, REML = T)

tb <- tidy(m, effects="fixed", conf.int=T, conf.method="boot")
tb$term <- c("Intercept","Height","Leaf N","Seed mass","SLA")
tb <- tb[-c(1), ]
tb$term <- factor(tb$term, level = c("Leaf N", "Seed mass", "SLA", "Height"))
p1<-ggplot(tb, aes(x=term,y=estimate, color = term)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_color_viridis_d() +
  #scale_color_manual(values = c("#482173FF", "#2D708EFF", "#51C56AFF","#85D54AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(x="", y="Effect size") +
  scale_y_continuous(breaks = c(-0.05,-0.025,0,0.025,0.05), labels = c("-0.050","-0.025","0.000","0.025","0.050")) +
  ggtitle("(a) MTV") +
  theme_bw() +
  theme(legend.position = "none",  plot.title=element_text(hjust=0), axis.text.x = element_text(size=7))

pa <- ggpredict(m, terms = c("median_Height")) %>% plot(show_data = T, show_ci = T, color = "#482173FF", dot_alpha = 0.1) + 
  theme_bw() + 
  labs(x="MTV height (m)", y = "Cover change (log ratio/year)") +
  ggtitle("(a)")

m <- lmer(RR ~ scale(COV_Height) + scale(COV_SLA) + scale(COV_SeedMass) + scale(COV_LeafN) +
            (1|species) + (1|region/summit/aspect), data = df, na.action = na.omit, REML = T)

tb <- tidy(m, effects="fixed", conf.int=T, conf.method="boot")
tb$term <- c("Intercept","Height","Leaf N","Seed mass","SLA")
tb <- tb[-c(1), ]
tb$term <- factor(tb$term, level = c("Leaf N", "Seed mass", "SLA", "Height"))
p2<-ggplot(tb, aes(x=term,y=estimate, color = term)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_color_viridis_d() +
  #scale_color_manual(values = c("#482173FF", "#2D708EFF", "#51C56AFF","#85D54AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(x="", y="Effect size") +
  scale_y_continuous(breaks = c(-0.05,-0.025,0,0.025,0.05), labels = c("-0.050","-0.025","0.000","0.025","0.050")) +
  ggtitle("(b) ITV") +
  theme_bw() +
  theme(legend.position = "none",  plot.title=element_text(hjust=0), axis.text.y = element_blank(), axis.text.x = element_text(size=7))

pb<- ggpredict(m, terms = c("SD_Height")) %>% plot(show_data = T, show_ci = T, color = "#2D708EFF", dot_alpha = 0.1) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(x = "ITV Height", y = "", color = "") +
  ggtitle("(b)")
pc<- ggpredict(m, terms = c("SD_SLA")) %>% plot(show_data = T, show_ci = T, color = "#2D708EFF", dot_alpha = 0.1) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(x = "ITV SLA", y = "", color = "") +
  #scale_x_continuous(breaks = c(0,5,10), labels = c("0","","10")) +
  ggtitle("(c)")
pd<- ggpredict(m, terms = c("SD_SeedMass")) %>% plot(show_data = T, show_ci = T, color = "#2D708EFF", dot_alpha = 0.1) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(x = "ITV Seed mass", y = "", color = "") +
  ggtitle("(d)")

comm <- test2 %>% select(region,summit,aspect,plot,species_name,RR)
colnames(comm) <- c("region","summit","aspect","plot","species","RR")
ind <- read.csv("shrubs_EIVE.csv", header = T)
df <- left_join(comm, ind, by = "species")

m <- lmer(RR ~ M + N + R + L + Te +
            (1|species) + (1|region/summit/aspect), data = df, na.action = na.omit, REML = T)
tb <- tidy(m, effects="fixed", conf.int=T, conf.method="boot")
tb$term <- c("Intercept","M","N","R","L","T")
tb <- tb[-c(1), ]
p3<-ggplot(tb, aes(x=term,y=estimate, color = term)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_color_viridis_d() +
  #scale_color_manual(values = c("#482173FF", "#2D708EFF", "#51C56AFF","#85D54AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(x="", y="Effect size") +
  scale_y_continuous(breaks = c(-0.05,-0.025,0,0.025,0.05), labels = c("-0.050","-0.025","0.000","0.025","0.050")) +
  ggtitle("(c) Niche position") +
  theme_bw() +
  theme(legend.position = "none",  plot.title=element_text(hjust=0), axis.text.x = element_text(size=7))

pe<- ggpredict(m, terms = c("R")) %>% plot(show_data = T, show_ci = T, color = "#51C56AFF", dot_alpha = 0.1) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(x = "R indicator value", y = "", color = "") +
  #scale_x_continuous(breaks = c(0,2.5,5), labels = c("0.0","2.5","5.0")) +
  ggtitle("(e)")

m <- lmer(RR ~ Mnw + Nnw + Rnw + Lnw + Tnw +
            (1|species) + (1|region/summit/aspect), data = df, na.action = na.omit, REML = T)
tb <- tidy(m, effects="fixed", conf.int=T, conf.method="boot")
tb$term <- c("Intercept","M","N","R ","L","T")
tb <- tb[-c(1), ]
p4<-ggplot(tb, aes(x=term,y=estimate, color = term)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_color_viridis_d() +
  #scale_color_manual(values = c("#482173FF", "#2D708EFF", "#51C56AFF","#85D54AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(x="", y="Effect size") +
  scale_y_continuous(breaks = c(-0.05,-0.025,0,0.025,0.05), labels = c("-0.050","-0.025","0.000","0.025","0.050")) +
  ggtitle("(d) Niche width") +
  theme_bw() +
  theme(legend.position = "none",  plot.title=element_text(hjust=0), axis.text.y = element_blank(), axis.text.x = element_text(size=7))

p <- ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, align = "hv")
png(file="Fig3_alternative.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

p <- p1 + p2 + plot_layout(axes = "collect", axis_titles = "collect") + p3 + p4 + plot_layout(ncol = 4) 

png(file="Fig3_alternative2.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

p <- pa + pb + pc + pd +pe + plot_layout(ncol = 5, axes = "collect") 

png(file="Fig4_alternative.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

## Shrub species in 5% of plots
df <- read.csv("shrubs.csv", header = T)
df <- df  %>% filter(species_name == "Alyssum_spinosum" |
                       species_name =="Arctostaphylos_alpinus" |
                       species_name =="Arctostaphylos_uva-ursi" |
                       species_name =="Betula_nana" |
                       species_name =="Calluna_vulgaris" |
                       species_name =="Cytisus_galianoi" |
                       species_name =="Daphne_striata" |
                       species_name =="Dryas_octopetala" |
                       species_name =="Empetrum_nigrum" |
                       species_name =="Erica_herbacea" |
                       species_name =="Euphorbia_acanthothamnos" |
                       species_name =="Helianthemum_oelandicum" |
                       species_name =="Juniperus_communis" |
                       species_name =="Loiseleuria_procumbens" |
                       species_name =="Polygala_chamaebuxus" |
                       species_name =="Prunus_prostrata" |
                       species_name =="Rhododendron_ferrugineum" |
                       species_name =="Rhododendron_myrtifolium" |
                       species_name =="Salix_herbacea" |
                       species_name =="Salix_reticulata" |
                       species_name =="Salix_retusa" |
                       species_name =="Saxifraga_retusa" |
                       species_name =="Sideritis_glacialis" |
                       species_name =="Thymus_capitatus" |
                       species_name =="Thymus_nervosus" |
                       species_name =="Thymus_nummularius" |
                       species_name =="Thymus_praecox" |
                       species_name =="Thymus_serpylloides" |
                       species_name =="Thymus_serpyllum" |
                       species_name =="Vaccinium_myrtillus" |
                       species_name =="Vaccinium_uliginosum" |
                       species_name =="Vaccinium_vitis-idaea")

test3 <- df %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(RR = log(cover_tot/lag(cover_tot))/(year - lag(year)), check = cover_tot + lag(cover_tot), check2 = cover_tot - lag(cover_tot)) %>%
  filter(check != 0)

test3 %>%
  group_by() %>%
  summarise(
    mean_RR = mean(RR),
    ci_lower = ci(RR)[2],
    ci_upper = ci(RR)[3])

test4 <- df %>%
  filter(DE == "d") %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(RR = log(cover_tot/lag(cover_tot))/(year - lag(year)), check = cover_tot + lag(cover_tot)) %>%
  filter(check != 0)

test4 %>%
  group_by() %>%
  summarise(
    mean_RR = mean(RR),
    ci_lower = ci(RR)[2],
    ci_upper = ci(RR)[3])

test5 <- df %>%
  filter(DE == "e") %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(RR = log(cover_tot/lag(cover_tot))/(year - lag(year)), check = cover_tot + lag(cover_tot)) %>%
  filter(check != 0)

test5 %>%
  group_by() %>%
  summarise(
    mean_RR = mean(RR),
    ci_lower = ci(RR)[2],
    ci_upper = ci(RR)[3])

p1<-ggplot(test3,aes(x=RR)) + 
  geom_histogram(bins = 100, aes(fill = RR >= 0, alpha = RR >= 0)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() +
  scale_fill_manual(values = c("black", "black")) +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "", y = "No. of plots") +
  ggtitle("All shrubs") +
  theme(legend.position = "none", plot.title = element_text(face="bold"))

p2<-ggplot(test4,aes(x=RR)) + 
  geom_histogram(bins = 100, aes(fill = RR >= 0, alpha = RR >= 0)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() +
  scale_fill_manual(values = c("#9999CC", "#9999CC")) +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Cover change (log ratio/year)", y = "") +
  ggtitle("Deciduous shrubs") +
  theme(legend.position = "none", plot.title = element_text(face="bold"), axis.title.y = element_blank())

p3<-ggplot(test5,aes(x=RR)) + 
  geom_histogram(bins = 100, aes(fill = RR >= 0, alpha = RR >= 0)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() +
  scale_fill_manual(values = c("#66CC99", "#66CC99")) +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "", y = "") +
  ggtitle("Evergreen shrubs") +
  theme(legend.position = "none", plot.title = element_text(face="bold"), axis.title.y = element_blank())

p <-p1 + p2 + p3 + plot_layout(ncol = 3, axes = "collect") 
png(file="Fig1_new.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

m <- lmer(RR~1+(1|region/summit/aspect), test3, na.action = na.omit, REML = T)
summary(m)
confint(m)

m <- lmer(RR~1+(1|region/summit/aspect), test4, na.action = na.omit, REML = T)
summary(m)
confint(m)

m <- lmer(RR~1+(1|region/summit/aspect), test5, na.action = na.omit, REML = T)
summary(m)
confint(m)

df <- read.csv("shrubs.csv", header = T)
df <- df  %>% filter(species_name == "Alyssum_spinosum" |
                       species_name =="Arctostaphylos_alpinus" |
                       species_name =="Arctostaphylos_uva-ursi" |
                       species_name =="Betula_nana" |
                       species_name =="Calluna_vulgaris" |
                       species_name =="Cytisus_galianoi" |
                       species_name =="Daphne_striata" |
                       species_name =="Dryas_octopetala" |
                       species_name =="Empetrum_nigrum" |
                       species_name =="Erica_herbacea" |
                       species_name =="Euphorbia_acanthothamnos" |
                       species_name =="Helianthemum_oelandicum" |
                       species_name =="Juniperus_communis" |
                       species_name =="Loiseleuria_procumbens" |
                       species_name =="Polygala_chamaebuxus" |
                       species_name =="Prunus_prostrata" |
                       species_name =="Rhododendron_ferrugineum" |
                       species_name =="Rhododendron_myrtifolium" |
                       species_name =="Salix_herbacea" |
                       species_name =="Salix_reticulata" |
                       species_name =="Salix_retusa" |
                       species_name =="Saxifraga_retusa" |
                       species_name =="Sideritis_glacialis" |
                       species_name =="Thymus_capitatus" |
                       species_name =="Thymus_nervosus" |
                       species_name =="Thymus_nummularius" |
                       species_name =="Thymus_praecox" |
                       species_name =="Thymus_serpylloides" |
                       species_name =="Thymus_serpyllum" |
                       species_name =="Vaccinium_myrtillus" |
                       species_name =="Vaccinium_uliginosum" |
                       species_name =="Vaccinium_vitis-idaea")
df <- df  %>% select(year, region, summit, aspect, plot, species_name, cover)
df <- spread(df, species_name, cover)
df <- df %>% replace(is.na(.), 0)
df <- gather(df, species_name, cover, 6:37)

test <- df %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  ungroup() %>%
  arrange(region, summit, aspect, plot, year)

test2 <- df %>%
  arrange(region, summit, aspect, plot, species_name, year) %>%
  filter(year == min(year) | year == max(year)) %>%
  group_by(region, summit, aspect, plot, species_name) %>%
  mutate(RR = log((cover+0.001)/lag((cover+0.001)))/(year - lag(year)), check = cover + lag(cover)) %>%
  filter(check != 0)

comm <- test2 %>% select(region,summit,aspect,plot,species_name,RR)
colnames(comm) <- c("region","summit","aspect","plot","species","RR")
traits <- read.csv("shrub_traits_ext.csv", header = T)
df <- left_join(comm, traits, by = "species")
phylo <- sp_list_df(sp_list = traits$species, taxon = "plant")
df <- left_join(df, phylo, by = "species")
write.csv(df, "shrub_traits.csv")

p<-ggplot(df, aes(reorder(POWO_accepted2, -RR, FUN=median), RR, color = DE)) +
  geom_boxplot() +
  geom_point(alpha=0.3) +
  #stat_summary(geom = "point",fun.y = "mean",col = "black",size = 3,shape = 24,fill = "black") +
  scale_color_manual(values = c("#9999CC", "#66CC99")) +
  theme_bw() +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.title = element_blank(), legend.box.margin = margin(0, 0, 0, 0), legend.text=element_text(size=7),
        legend.background = element_rect(fill="white",size=0.5,linetype="solid",colour ="gray"),
        legend.key.size = unit(0.2, "cm"), legend.spacing.y = unit(0, "pt"), legend.spacing.x = unit(0, "pt"), 
        axis.text.x = element_text(angle=50, vjust=1, hjust=1, face = "italic")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x="", y="Cover change (log-ratio/year)", color="")

png(file="Fig2_new.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

results <- df %>%
  group_by(POWO_accepted2, DE) %>%
  summarise(
    mean_RR = mean(RR),
    ci_lower = ci(RR)[2],
    ci_upper = ci(RR)[3])

p<-ggplot(results, aes(reorder(POWO_accepted2, -mean_RR, FUN=mean), mean_RR, color = DE)) +
  geom_point() +
  geom_errorbar(aes(x =reorder(POWO_accepted2, -mean_RR, FUN=mean), ymin = ci_lower, ymax = ci_upper), width = 0) +
  scale_color_manual(values = c("#9999CC", "#66CC99")) +
  theme_bw() +
  theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.title = element_blank(), legend.box.margin = margin(0, 0, 0, 0), legend.text=element_text(size=7),
        legend.background = element_rect(fill="white",size=0.5,linetype="solid",colour ="gray"),
        legend.key.size = unit(0.2, "cm"), legend.spacing.y = unit(0, "pt"), legend.spacing.x = unit(0, "pt"), 
        axis.text.x = element_text(angle=50, vjust=1, hjust=1, face = "italic")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x="", y="Cover change (log-ratio/year)", color="")

png(file="Fig2_new2.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

#df <- left_join(results, traits, by = "POWO_accepted2")

m <- lmer(RR ~ scale(log(mean_Height)) + scale(mean_SLA) + scale(log(mean_SeedMass)) + scale(mean_LeafN) + 
            (1|family/genus) + (1|region/summit/aspect), data = df, na.action = na.omit, REML = T)

tb <- tidy(m, effects="fixed", conf.int=T, conf.method="boot")
tb$term <- c("Intercept","Height","SLA","Seed mass","Leaf N")
tb <- tb[-c(1), ]
tb$term <- factor(tb$term, level = c("Leaf N", "Seed mass", "SLA", "Height"))
p1<-ggplot(tb, aes(x=term,y=estimate, color = term)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_color_viridis_d() +
  #scale_color_manual(values = c("#482173FF", "#2D708EFF", "#51C56AFF","#85D54AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(x="", y="Effect size") +
  scale_y_continuous(breaks = c(-0.05,0,0.05), labels = c("-0.05","0.00","0.05"), limits = c(-0.07,0.07)) +
  ggtitle("(a) Mean traits") +
  theme_bw() +
  theme(legend.position = "none",  plot.title=element_text(hjust=0), axis.text.x = element_text(size=7))

pa <- ggpredict(m, terms = c("mean_Height")) %>% plot(show_data = T, show_ci = T, color = "#FDE725FF", dot_alpha = 0.1) + 
  theme_bw() + 
  labs(x="Mean plant height (m)", y = "Cover change (log ratio/year)") +
  ggtitle("(a)")
pb <- ggpredict(m, terms = c("mean_LeafN")) %>% plot(show_data = T, show_ci = T, color = "#7AD151FF", dot_alpha = 0.1) + 
  theme_bw() + 
  theme(axis.title.y = element_blank()) +
  labs(x="Mean leaf N content (mg/g)", y = "") +
  ggtitle("(b)")

m <- lmer(RR ~ scale(COV_Height) + scale(log(COV_SLA)) + scale(log(COV_SeedMass)) + scale(COV_LeafN) +
            (1|family/genus) + (1|region/summit/aspect), data = df, na.action = na.omit, REML = T)

tb <- tidy(m, effects="fixed", conf.int=T, conf.method="boot")
tb$term <- c("Intercept","Height","SLA","Seed mass","Leaf N")
tb <- tb[-c(1), ]
tb$term <- factor(tb$term, level = c("Leaf N", "Seed mass", "SLA", "Height"))
p2<-ggplot(tb, aes(x=term,y=estimate, color = term)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_color_viridis_d() +
  #scale_color_manual(values = c("#482173FF", "#2D708EFF", "#51C56AFF","#85D54AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(x="", y="Effect size") +
  scale_y_continuous(breaks = c(-0.05,0,0.05), labels = c("-0.05","0.00","0.05"), limits = c(-0.07,0.07)) +
  ggtitle("(b) ITV") +
  theme_bw() +
  theme(legend.position = "none",  plot.title=element_text(hjust=0), axis.text.y = element_blank(), axis.text.x = element_text(size=7))


comm <- test2 %>% select(region,summit,aspect,plot,species_name,RR)
colnames(comm) <- c("region","summit","aspect","plot","species","RR")
ind <- read.csv("shrubs_EIVE_ext.csv", header = T)
df <- left_join(comm, ind, by = "species")
df <- left_join(df, phylo, by = "species")
write.csv(df,"shrubs_ind.csv")

m <- lmer(RR ~ M + N + R + L + Te +
            (1|family/genus) + (1|region/summit/aspect), data = df, na.action = na.omit, REML = T)

tb <- tidy(m, effects="fixed", conf.int=T, conf.method="boot")
tb$term <- c("Intercept","M","N","R","L","T")
tb <- tb[-c(1), ]
p3<-ggplot(tb, aes(x=term,y=estimate, color = term)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_color_viridis_d() +
  #scale_color_manual(values = c("#482173FF", "#2D708EFF", "#51C56AFF","#85D54AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(x="", y="Effect size") +
  scale_y_continuous(breaks = c(-0.05,0,0.05), labels = c("-0.05","0.00","0.05"), limits = c(-0.07,0.07)) +
  ggtitle("(c) Niche optimum") +
  theme_bw() +
  theme(legend.position = "none",  plot.title=element_text(hjust=0), axis.text.x = element_text(size=7))

pc<- ggpredict(m, terms = c("L")) %>% plot(show_data = T, show_ci = T, color = "#2D708EFF", dot_alpha = 0.1) + 
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(x = "Light optimum", y = "", color = "") +
  ggtitle("(c)")

m <- lmer(RR ~ Mnw + Nnw + Rnw + Lnw + Tnw +
            (1|family/genus) + (1|region/summit/aspect), data = df, na.action = na.omit, REML = T)

tb <- tidy(m, effects="fixed", conf.int=T, conf.method="boot")
tb$term <- c("Intercept","M","N","R ","L","T")
tb <- tb[-c(1), ]
p4<-ggplot(tb, aes(x=term,y=estimate, color = term)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  scale_color_viridis_d() +
  #scale_color_manual(values = c("#482173FF", "#2D708EFF", "#51C56AFF","#85D54AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() + 
  labs(x="", y="Effect size") +
  scale_y_continuous(breaks = c(-0.05,0,0.05), labels = c("-0.05","0.00","0.05"), limits = c(-0.07,0.07)) +
  ggtitle("(d) Niche width") +
  theme_bw() +
  theme(legend.position = "none",  plot.title=element_text(hjust=0), axis.text.y = element_blank(), axis.text.x = element_text(size=7))

p <- p1 + p2 + plot_layout(axes = "collect", axis_titles = "collect") + p3 + p4 + plot_layout(ncol = 4) 

png(file="Fig3_new.png",width=20,height=7,units="cm",res=300, pointsize=12)
p
dev.off()

p <- pa + pb + pc + plot_layout(axes = "collect", ncol = 3) 

png(file="Fig4_new.png",width=20,height=10,units="cm",res=300, pointsize=12)
p
dev.off()

## Shrub species in 5% of plots
df <- read.csv("shrubs.csv", header = T)
df <- df  %>% filter(species_name == "Alyssum_spinosum" |
                       species_name =="Arctostaphylos_alpinus" |
                       species_name =="Arctostaphylos_uva-ursi" |
                       species_name =="Betula_nana" |
                       species_name =="Calluna_vulgaris" |
                       species_name =="Cytisus_galianoi" |
                       species_name =="Daphne_striata" |
                       species_name =="Dryas_octopetala" |
                       species_name =="Empetrum_nigrum" |
                       species_name =="Erica_herbacea" |
                       species_name =="Euphorbia_acanthothamnos" |
                       species_name =="Helianthemum_oelandicum" |
                       species_name =="Juniperus_communis" |
                       species_name =="Loiseleuria_procumbens" |
                       species_name =="Polygala_chamaebuxus" |
                       species_name =="Prunus_prostrata" |
                       species_name =="Rhododendron_ferrugineum" |
                       species_name =="Rhododendron_myrtifolium" |
                       species_name =="Salix_herbacea" |
                       species_name =="Salix_reticulata" |
                       species_name =="Salix_retusa" |
                       species_name =="Saxifraga_retusa" |
                       species_name =="Sideritis_glacialis" |
                       species_name =="Thymus_capitatus" |
                       species_name =="Thymus_nervosus" |
                       species_name =="Thymus_nummularius" |
                       species_name =="Thymus_praecox" |
                       species_name =="Thymus_serpylloides" |
                       species_name =="Thymus_serpyllum" |
                       species_name =="Vaccinium_myrtillus" |
                       species_name =="Vaccinium_uliginosum" |
                       species_name =="Vaccinium_vitis-idaea")

test3 <- df %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(RR = (cover_tot - lag(cover_tot))/(year - lag(year)), check = cover_tot + lag(cover_tot), check2 = cover_tot - lag(cover_tot)) %>%
  filter(check != 0)

test3 %>%
  group_by() %>%
  summarise(
    mean_RR = mean(RR),
    ci_lower = ci(RR)[2],
    ci_upper = ci(RR)[3])

test4 <- df %>%
  filter(DE == "d") %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(RR = (cover_tot - lag(cover_tot))/(year - lag(year)), check = cover_tot + lag(cover_tot)) %>%
  filter(check != 0)

test4 %>%
  group_by() %>%
  summarise(
    mean_RR = mean(RR),
    ci_lower = ci(RR)[2],
    ci_upper = ci(RR)[3])

test5 <- df %>%
  filter(DE == "e") %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year)) %>%
  mutate(RR = (cover_tot - lag(cover_tot))/(year - lag(year)), check = cover_tot + lag(cover_tot)) %>%
  filter(check != 0)

test5 %>%
  group_by() %>%
  summarise(
    mean_RR = mean(RR),
    ci_lower = ci(RR)[2],
    ci_upper = ci(RR)[3])

## Shrub species in 5% of plots
df <- read.csv("shrubs.csv", header = T)
df <- df  %>% filter(species_name == "Alyssum_spinosum" |
                       species_name =="Arctostaphylos_alpinus" |
                       species_name =="Arctostaphylos_uva-ursi" |
                       species_name =="Betula_nana" |
                       species_name =="Calluna_vulgaris" |
                       species_name =="Cytisus_galianoi" |
                       species_name =="Daphne_striata" |
                       species_name =="Dryas_octopetala" |
                       species_name =="Empetrum_nigrum" |
                       species_name =="Erica_herbacea" |
                       species_name =="Euphorbia_acanthothamnos" |
                       species_name =="Helianthemum_oelandicum" |
                       species_name =="Juniperus_communis" |
                       species_name =="Loiseleuria_procumbens" |
                       species_name =="Polygala_chamaebuxus" |
                       species_name =="Prunus_prostrata" |
                       species_name =="Rhododendron_ferrugineum" |
                       species_name =="Rhododendron_myrtifolium" |
                       species_name =="Salix_herbacea" |
                       species_name =="Salix_reticulata" |
                       species_name =="Salix_retusa" |
                       species_name =="Saxifraga_retusa" |
                       species_name =="Sideritis_glacialis" |
                       species_name =="Thymus_capitatus" |
                       species_name =="Thymus_nervosus" |
                       species_name =="Thymus_nummularius" |
                       species_name =="Thymus_praecox" |
                       species_name =="Thymus_serpylloides" |
                       species_name =="Thymus_serpyllum" |
                       species_name =="Vaccinium_myrtillus" |
                       species_name =="Vaccinium_uliginosum" |
                       species_name =="Vaccinium_vitis-idaea")

test3 <- df %>%
  group_by(year, region, summit, aspect, plot) %>%
  summarise(cover_tot = sum(cover)) %>%
  arrange(region, summit, aspect, plot, year) %>%
  group_by(region, summit, aspect, plot) %>%
  filter(year == min(year) | year == max(year))
test3$year <- as.factor(test3$year)
ggplot(test3, aes(x=cover_tot, fill=year)) + geom_density(alpha=.3)
