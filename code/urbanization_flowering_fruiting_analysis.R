

library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)
library(cAIC4)
library(ggeffects)
library(sjPlot)
library(sf)
library(cowplot)


#----------------------------------------
fl1=read_csv(file="data/flowering_modeling_data_20221123.csv")
fr1=read_csv(file="data/fruiting_modeling_data_20221123.csv")

summary(fl1)
summary(fr1)
colnames(fl1)


# standardized data for modeling
fl1m=fl1 %>% select(doy, sc2, binomial_species, LT_tm_spr, LT_ppt_spr, tm_win_anm, tm_spr_anm, ppt_spr_anm, pop_den)
fl1m=cbind.data.frame(fl1m[, 1:3], scale(fl1m[, 4:8]), scale(log(fl1m[,9])))
summary(fl1m)



fr1m=fr1 %>% select(doy, sc2, binomial_species, LT_tm_spr, LT_ppt_spr, tm_win_anm, tm_spr_anm, ppt_spr_anm, pop_den)
fr1m=cbind.data.frame(fr1m[, 1:3], scale(fr1m[, 4:8]), scale(log(fr1m[,9])))
summary(fr1m)



#------- flowering - modeling -----------------
# check correlations
cor(na.omit(fl1[, c(1:3,8:11,16:23)]))
cor(na.omit(fr1[, c(1:3,8:11,16:23)]))



# no urbanization
# determine the optimal structure of fixed effects

fm1 <- lmer(doy ~ county.lat + alt + 
              tm_spr_anm + tm_win_anm*ppt_spr_anm + 
              (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(fm1)

fm2 <- lmer(doy ~ county.lat + alt + 
              tm_spr_anm + ppt_spr_anm + tm_win_anm:ppt_spr_anm+
              (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(fm2)

AICc(fm1, fm2)


fm3 <- lmer(doy ~ LT_MAT + LT_MAP + 
              tm_spr_anm + tm_win_anm*ppt_spr_anm + 
              (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(fm3)

fm4 <- lmer(doy ~ LT_tm_spr + LT_ppt_spr + 
              tm_spr_anm + tm_win_anm*ppt_spr_anm + 
              (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(fm4)

AICc(fm1, fm3, fm4)

# best model
fm.fl <- lmer(doy ~ LT_tm_spr + LT_ppt_spr + 
                tm_spr_anm + ppt_spr_anm + tm_win_anm:ppt_spr_anm+
                (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2), 
              data = fl1m, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(fm.fl)
r.squaredGLMM(fm.fl)  


fm.fl_coef=coef(fm.fl)
fm.fl_sum=summary(fm.fl)
fmfl_cffixed=as.data.frame(fm.fl_sum$coefficients)
fmfl_cffixed$variable=row.names(fmfl_cffixed)
fmfl_cffixed$variable[1]="Intercept"
fmfl_cffixed$model="No_urb"
fmfl_cfsc=fm.fl_coef$sc2
fmfl_cfspp=fm.fl_coef$binomial_species
fmfl_cfrdm=rbind.data.frame(fmfl_cfspp,fmfl_cfsc)
fmfl_cfrdm$var=c(rep("species", nrow(fmfl_cfspp)),rep("state_county",nrow(fmfl_cfsc)))
fmfl_cfrdm$group=row.names(fmfl_cfrdm)

cor.test(fmfl_cfspp$`(Intercept)`,fmfl_cfspp$tm_spr_anm)  
cor.test(fmfl_cfspp$`(Intercept)`,fmfl_cfspp$tm_win_anm) 



# urbanization model

tm1 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              pop_den*tm_spr_anm + pop_den*tm_win_anm + 
              pop_den*ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm1)

tm2 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              pop_den*tm_spr_anm + pop_den*tm_win_anm + 
              ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm2)

tm3 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              tm_spr_anm + pop_den*tm_win_anm + 
              ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm3)

tm4 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              tm_spr_anm + tm_win_anm + 
              ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm4)

tm5 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              tm_spr_anm + 
              ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm5)


AICc(tm1, tm2, tm3, tm4, tm5)

tm6 <- lmer(doy ~ LT_tm_spr+LT_ppt_spr + pop_den:LT_tm_spr + pop_den:LT_ppt_spr + 
              tm_spr_anm + ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fl1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm6)

AICc(tm1, tm5, tm6)


# best model
tm.fl <- lmer(doy ~ LT_tm_spr + LT_ppt_spr + 
                tm_spr_anm +LT_tm_spr:pop_den+ LT_ppt_spr:pop_den+
                ppt_spr_anm + tm_win_anm:ppt_spr_anm +
                (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
              data = fl1m, REML=T, control = lmerControl(optimizer="bobyqa")) 

summary(tm.fl)
r.squaredGLMM(tm.fl)  


tm.fl_coef=coef(tm.fl)
tm.fl_sum=summary(tm.fl)
tmfl_cffixed=as.data.frame(tm.fl_sum$coefficients)
tmfl_cffixed$variable=row.names(tmfl_cffixed)
tmfl_cffixed$variable[1]="Intercept"
tmfl_cffixed$model="Urb"
tmfl_cfsc=tm.fl_coef$sc2
tmfl_cfspp=tm.fl_coef$binomial_species
tmfl_cfrdm=rbind.data.frame(tmfl_cfspp,tmfl_cfsc)
tmfl_cfrdm$var=c(rep("species", nrow(tmfl_cfspp)),rep("state_county",nrow(tmfl_cfsc)))
tmfl_cfrdm$group=row.names(tmfl_cfrdm)


cor.test(tmfl_cfspp$`(Intercept)`,tmfl_cfspp$tm_spr_anm)  
cor.test(tmfl_cfspp$`(Intercept)`,tmfl_cfspp$tm_win_anm)  



#------- fruiting - modeling ------------------------

# no urbanization
# determine the optimal structure of fixed effects

fm1 <- lmer(doy ~ county.lat + alt + 
              tm_spr_anm + tm_win_anm*ppt_spr_anm + 
              (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2), 
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(fm1)

fm2 <- lmer(doy ~ LT_MAT + LT_MAP +
              tm_spr_anm + tm_win_anm*ppt_spr_anm + 
              (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2), 
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(fm2)

fm3 <- lmer(doy ~ LT_tm_spr + LT_ppt_spr +
              tm_spr_anm + tm_win_anm*ppt_spr_anm + 
              (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2), 
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(fm3)

AICc(fm1, fm2, fm3)
 

fm4 <- lmer(doy ~ LT_tm_spr*LT_ppt_spr +
              tm_spr_anm + tm_win_anm:ppt_spr_anm + ppt_spr_anm+
              (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2), 
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(fm4)



# best model
fm.fr <- lmer(doy ~ LT_tm_spr + LT_ppt_spr +
                tm_spr_anm + ppt_spr_anm + 
                tm_win_anm:ppt_spr_anm + 
                (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2),
              data = fr1m, REML=T, control = lmerControl(optimizer="bobyqa"))
summary(fm.fr)
r.squaredGLMM(fm.fr)  


fm.fr_coef=coef(fm.fr)
fm.fr_sum=summary(fm.fr)
#str(fm.fr_sum)
fmfr_cffixed=as.data.frame(fm.fr_sum$coefficients)
fmfr_cffixed$variable=row.names(fmfr_cffixed)
fmfr_cffixed$variable[1]="Intercept"
fmfr_cffixed$model="No_urb"
fmfr_cfsc=fm.fr_coef$sc2
fmfr_cfspp=fm.fr_coef$binomial_species
fmfr_cfrdm=rbind.data.frame(fmfr_cfspp,fmfr_cfsc)
fmfr_cfrdm$var=c(rep("species", nrow(fmfr_cfspp)),rep("state_county",nrow(fmfr_cfsc)))
fmfr_cfrdm$group=row.names(fmfr_cfrdm)

cor.test(fmfr_cfspp$`(Intercept)`,fmfr_cfspp$tm_spr_anm)  
cor.test(fmfr_cfspp$`(Intercept)`,fmfr_cfspp$tm_win_anm)  



# urbanization model
tm1 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              pop_den*tm_spr_anm + pop_den*tm_win_anm + 
              pop_den*ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm1)

tm2 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              pop_den*tm_spr_anm + pop_den*tm_win_anm + 
              ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm2)

tm3 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              pop_den*tm_spr_anm + tm_win_anm + 
              ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm3)

tm4 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              tm_spr_anm + tm_win_anm + 
              ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2), 
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm4)

tm5 <- lmer(doy ~ pop_den*LT_tm_spr + pop_den*LT_ppt_spr + 
              tm_spr_anm +  
              ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2),  
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm5)

tm6 <- lmer(doy ~ LT_tm_spr + pop_den*LT_ppt_spr +  
              tm_spr_anm + ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2),
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm6)

tm7 <- lmer(doy ~ LT_tm_spr + pop_den:LT_ppt_spr + pop_den +LT_tm_spr:LT_ppt_spr+
              tm_spr_anm + ppt_spr_anm + tm_win_anm:ppt_spr_anm +
              (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2),
            data = fr1m, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(tm7)

AICc(tm1, tm4, tm5, tm6, tm7)


# best model
tm.fr <- lmer(doy ~ LT_tm_spr + pop_den:LT_ppt_spr + pop_den +
                tm_spr_anm + ppt_spr_anm + tm_win_anm:ppt_spr_anm +
                (1+tm_spr_anm+tm_win_anm|binomial_species)+ (1|sc2),
              data = fr1m, REML=T, control = lmerControl(optimizer="bobyqa")) 

summary(tm.fr)
r.squaredGLMM(tm.fr)  


tm.fr_coef=coef(tm.fr)
tm.fr_sum=summary(tm.fr)
#str(tm.fr_sum)
tmfr_cffixed=as.data.frame(tm.fr_sum$coefficients)
tmfr_cffixed$variable=row.names(tmfr_cffixed)
tmfr_cffixed$variable[1]="Intercept"
tmfr_cffixed$model="Urb"
tmfr_cfsc=tm.fr_coef$sc2
tmfr_cfspp=tm.fr_coef$binomial_species
tmfr_cfrdm=rbind.data.frame(tmfr_cfspp,tmfr_cfsc)
tmfr_cfrdm$var=c(rep("species", nrow(tmfr_cfspp)),rep("state_county",nrow(tmfr_cfsc)))
tmfr_cfrdm$group=row.names(tmfr_cfrdm)


cor.test(tmfr_cfspp$`(Intercept)`,tmfr_cfspp$tm_spr_anm)  
cor.test(tmfr_cfspp$`(Intercept)`,tmfr_cfspp$tm_win_anm)  


# all models and model outputs
save(fm.fl, tm.fl, fm.fr, tm.fr, file="results/models_flowering_fruiting_rv2.Rdata")

cffixed=rbind.data.frame(fmfl_cffixed, tmfl_cffixed, fmfr_cffixed, tmfr_cffixed)
cffixed= cffixed %>% 
  rename(coef=Estimate, sd=`Std. Error`, t_value=`t value`, p_value=`Pr(>|t|)`) %>% 
  mutate(lwr=coef-sd, upr=coef+sd)
cffixed$phenology=c(rep("flowering",nrow(fmfl_cffixed)+nrow(tmfl_cffixed)), 
                    rep("fruiting",nrow(fmfr_cffixed)+nrow(tmfr_cffixed)))

cfrdm=list(fmfl_cfrdm, tmfl_cfrdm, fmfr_cfrdm, tmfr_cfrdm)

save(cffixed, cfrdm, file="results/modelcoef_flowering_fruiting_rv2.Rdata")





#- Figures ----- model output --------------------------------------------------------
library(cowplot)

load("results/models_flowering_fruiting_rv2.Rdata")
load("results/modelcoef_flowering_fruiting_rv2.Rdata")

cffixed = cffixed %>% filter(variable != "Intercept")

unique(cffixed$variable)

cffixed[cffixed$variable=="LT_ppt_spr:pop_den",6]="pop_den:LT_ppt_spr"
cffixed[cffixed$variable=="LT_tm_spr:pop_den",6]="pop_den:LT_tm_spr"
unique(cffixed$variable)

cffixed$variable=factor(cffixed$variable, 
                        levels=c("LT_tm_spr", "LT_ppt_spr", "tm_spr_anm", 
                                 "ppt_spr_anm", "ppt_spr_anm:tm_win_anm",
                                 "pop_den:LT_tm_spr", "pop_den:LT_ppt_spr", "pop_den"))

labels=c("LT_tm_spr"="LT_tm_sp","LT_ppt_spr"="LT_ppt_sp", 
         "tm_spr_anm"="Tm_sp_anm","ppt_spr_anm"="ppt_sp_anm",
         "ppt_spr_anm:tm_win_anm"="Tm_wn_anm:ppt_sp_anm",
         "pop_den:LT_tm_spr"="pop_den:LT_Tm_sp",
         "pop_den:LT_ppt_spr"="pop_den:LT_ppt_sp", 
         "pop_den"="pop_den")


# Figure 2
tiff(filename="figures_rv/flowering_fruiting_2model_coef_rv.tif", width= 2400, height= 2000, res=300, compression= "lzw")
ggplot(cffixed, aes(
  x = variable,
  y = coef, ymin = lwr, ymax = upr,
  #col = model
  shape = model
)) + 
  geom_hline(yintercept = 0, alpha = .5, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = -.5), size=1) + # minus to reverse factor levels, will lead to warning
  scale_x_discrete(limits = rev(levels(cffixed$variable)), 
                   labels=labels) +

  scale_shape_manual(values=c(21,22), labels=c("No urbanization","Urbanization")) +
  labs(x = "", y = "Coefficient value (estimate and standard error)", shape = "") +
  facet_wrap(~phenology) +
  coord_flip() + 
  theme_cowplot()+
  theme(legend.position="top")+
  guides(shape = guide_legend(override.aes = list(linetype = 0)))
dev.off()


# combine random slope plots
rdm2.fl = cfrdm[[2]] %>% filter(var=="species") %>%
  select(1,2,5,11) %>% 
  pivot_longer(col=c(-group, -`(Intercept)`), names_to = "variable", values_to="coef")


rdm2.fr = cfrdm[[4]] %>% filter(var=="species") %>%
  select(c(1,2,5,10)) %>% 
  pivot_longer(col=c(-group, -`(Intercept)`), names_to = "variable", values_to="coef")


rdm2=rbind.data.frame(rdm2.fl, rdm2.fr)
rdm2$pheno=c(rep("flowering",nrow(rdm2.fl)),rep("fruiting",nrow(rdm2.fr)))
rdm2$sig=1
rdm2$sig[rdm2$variable=="tm_win_anm"]=0

var_names=c("tm_spr_anm"="spring temp anomaly", "tm_win_anm"="winter temp anomaly",
            "flowering"="flowering", "fruiting"="fruiting")

# Figure 3
tiff(filename="figures_rv/flowering_fruiting_random_coef_Urb_rv.tif", width= 2100, height= 1800, res=300, compression= "lzw")
ggplot(rdm2, aes(
  x = `(Intercept)`,
  y = `coef`
)) +
  geom_point() +
  geom_smooth(method="lm", aes(col=factor(sig)), se=F)+
  scale_color_manual(values=c("blue","red"))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(variable~pheno, scales="free_y", labeller = as_labeller(var_names)) +
  labs(x="random intercept", y="phenological sensitivity across species")+
  theme_cowplot() +
  theme(legend.position="none")
dev.off()




#----------------- effects plot ---------------------------------------------------

load("results/models_flowering_fruiting_rv2.Rdata")


# Figure 4
line_df <- data.frame(x=fl1m$pop_den, y=155, xend = fl1m$pop_den, yend = 157)

p1=plot_model(tm.fl, type="pred", terms=c("pop_den", "LT_tm_spr[-0.9976, -0.0059, 1.0006]"), ci.lvl = F,
              title = "",
              colors = c("blue","purple","red")) +
  scale_color_manual(values=c("blue","purple","red"), labels=c(0.7, 7.4, 14.2),
                     name="long-term average \nof spr temp (\u00b0C)") +
  scale_x_continuous(breaks=c(  -2.72, -1.06, 0.6, 2.27, 3.93),
                     labels=c( 1, 10, 100, 1000, 10000), limits=c(-4.4, 6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(155, 210), expand = c(0, 0)) +
  geom_segment(data=line_df, aes(x=x, y=y, xend =xend, yend = yend),inherit.aes = FALSE)+
  labs(x=bquote("population "~"km"^-2), y="flowering date (day of year)")+
  theme_cowplot()+
  theme(legend.position=c(0.65, 0.85), legend.background=element_blank())


p2=plot_model(tm.fl, type="pred", terms=c("pop_den", "LT_ppt_spr[-1.0098, 0.0034, 1.0165]"), ci.lvl = F,
              title = "",
              colors = c("deepskyblue","dodgerblue2","blue")) +
  scale_color_manual(values=c("deepskyblue","dodgerblue2","blue"), labels=c(255, 295, 335),
                     name="long-term average \nof spring prcp (mm)") +
  scale_x_continuous(breaks=c(  -2.72, -1.06, 0.6, 2.27, 3.93),
                     labels=c( 1, 10, 100, 1000, 10000), limits=c(-4.4, 6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(155, 210), expand = c(0, 0)) +
  geom_segment(data=line_df, aes(x=x, y=y, xend =xend, yend = yend),inherit.aes = FALSE)+
  labs(x=bquote("population "~"km"^-2), y="flowering date (day of year)")+
  theme_cowplot()+
  theme(legend.position=c(0.55, 0.85),
        legend.background=element_blank())

line_df <- data.frame(x=fr1m$pop_den, y=180, xend = fr1m$pop_den, yend = 181.5)

p3=plot_model(tm.fr, type="pred", terms=c("pop_den", "LT_ppt_spr[-0.9957, -0.0001, 0.9954]"), ci.lvl = F,
              title = "",
              colors = c("deepskyblue","dodgerblue2","blue")) +
  scale_color_manual(values=c("deepskyblue","dodgerblue2","blue"), labels=c(255, 293, 331),
                     name="long-term average \nof spring prcp (mm)") +
  scale_x_continuous(breaks=c(  -2.72, -1.06, 0.6, 2.27, 3.93),
                     labels=c( 1, 10, 100, 1000, 10000), limits=c(-4.4, 6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(180, 220), expand = c(0, 0)) +
  geom_segment(data=line_df, aes(x=x, y=y, xend =xend, yend = yend),inherit.aes = FALSE)+
  labs(x=bquote("population "~"km"^-2), y="fruiting date (day of year)")+
  theme_cowplot()+
  theme(legend.position=c(0.55, 0.85),
        legend.background=element_blank())


plot.all=cowplot::plot_grid(p1, p2, p3, labels=c("(a)","(b)","(c)"), label_size=14, ncol=3, nrow=1, rel_widths=c(1, 1, 1))

ggsave(filename = "figures_rv/flowering_fruiting_effectplot_interactions3_rv.tif", plot.all, device="tiff",
       dpi=300, height= 5, width=15, compression="lzw")


# Figure S4

line_df <- data.frame(x=fl1m$ppt_spr_anm, y=175, xend = fl1m$ppt_spr_anm, yend = 175.5)

tiff(filename="figures_rv/flowering_effectplot_anm_interaction_rv.tif", width= 1800, height= 1600, res=300, compression= "lzw")
plot_model(tm.fl, type="pred", terms=c("ppt_spr_anm", "tm_win_anm[-0.993, 0.0021, 0.997]"), 
           ci.lvl = F, title = "",
           legend.title = paste("scaled winter temp"), colors = c("blue","purple","red")) +
  scale_color_manual(values=c("blue","purple","red"), labels=c(-2.2, -0.1, 2.0),
                     name="winter temp anomaly (\u00b0C)") +
  labs(x="spring precipitation proportional anomaly (%)", y="flowering date (day of year)")+
  scale_x_continuous(breaks=c(-3.31, -1.66, -0.0187, 1.63, 3.27, 4.91),
                     labels=c(-1, -0.5, 0, 0.5, 1, 1.5)) +
  scale_y_continuous(limits = c(175, 195), expand = c(0, 0)) +
  geom_segment(data=line_df, aes(x=x, y=y, xend =xend, yend = yend),inherit.aes = FALSE)+
  theme_cowplot()+
  theme(legend.position=c(0.55, 0.85), legend.background=element_blank())
dev.off()





facet.labs=c("pop_den=10", "pop_den=100", "pop_den=1000")
names(facet.labs)=c("2.3", "4.6", "6.9")

flpt.pred=ggpredict(mod, terms= c("ppt_spr","tm_win_anm","log_popden[2.3, 4.6, 6.9]"))
fltp.pred=ggpredict(mod, terms= c("tm_win_anm","ppt_spr","log_popden[2.3, 4.6, 6.9]"))

p1=ggplot(flpt.pred, aes(x=x, y=predicted, col=group)) +
  geom_line() +
  facet_wrap(~facet, labeller = labeller(facet=facet.labs))+
  scale_color_manual(values=c("blue","purple","red"), labels=c(-3.2, 3.8, 10.8),
                     name="winter temp (\u00b0C)") +
  scale_x_continuous(breaks=c(-3, -1, 1, 3, 5),
                     labels=c(0, 200, 400, 600, 800)) +
  labs(x="spring precipitation (mm)", y="flowering date (day of year)")+
  theme_cowplot()  

p2=ggplot(fltp.pred, aes(x=x, y=predicted, col=group)) +
  geom_line() +
  facet_wrap(~facet, labeller = labeller(facet=facet.labs))+
  scale_color_manual(values=c("deepskyblue","dodgerblue2","blue"), labels=c(200, 300, 400),
                     name="spring prcp (mm)") +
  scale_x_continuous(breaks=c(-2.24, -1.39, -0.538, 0.315, 1.168, 2.02, 2.87 ),
                     labels=c(-12, -6, 0, 6, 12, 18, 24)) +
  labs(x="winter temperature (\u00b0C)", y="flowering date (day of year)")+
  theme_cowplot()  

plot.all=cowplot::plot_grid(p1, p2, labels=c("a","b"), label_size=14, ncol=1, nrow=2, rel_heights=c(1, 1))

ggsave(filename = "figures/flowering_effectplot_interactions3.tif", plot.all, device="tiff",
       dpi=300, height= 10, width=10, compression="lzw")




# ------------------ frost risk ------------------------------

# select species with at least one specimen collected before last frost 
fl2=read_csv(file="data/flowering_frost_modeling_data_rv.csv")

f.spp = fl2 %>% filter(frostrisk>0) %>% select(binomial_species) %>% distinct() 

fl3=fl2 %>% filter(binomial_species %in% f.spp$binomial_species) %>% na.omit()

sppn=fl3 %>% group_by(binomial_species) %>% summarise(n=n())

fl3$id=1:nrow(fl3)

# separate training and testing data
set.seed(123)

sppn$sample_n=as.integer(sppn$n*0.9)

fl3_tr <- fl3 %>% group_split(binomial_species) %>% map2_dfr(sppn$sample_n, ~ slice_sample(.x, n = .y))
fl3_ts <- fl3 %>% filter(! id %in% fl3_tr$id)

# standardized data for modeling
fl3m_tr=fl3_tr %>% select(frostrisk, sc2, binomial_species, LT_tm_spr, LT_ppt_spr, tm_win_anm, tm_spr_anm, ppt_spr_anm, pop_den)
fl3m_tr=cbind.data.frame(fl3m_tr[, 1:3], scale(fl3m_tr[, 4:8]), scale(log(fl3m_tr[,9])))
summary(fl3m_tr)


fl3m_ts=fl3_ts %>% select(frostrisk, sc2, binomial_species)
fl3m_ts$LT_tm_spr=(fl3_ts$LT_tm_spr - mean(fl3_tr$LT_tm_spr))/sd(fl3_tr$LT_tm_spr)
fl3m_ts$LT_ppt_spr=(fl3_ts$LT_ppt_spr - mean(fl3_tr$LT_ppt_spr))/sd(fl3_tr$LT_ppt_spr)
fl3m_ts$tm_win_anm=(fl3_ts$tm_win_anm - mean(fl3_tr$tm_win_anm))/sd(fl3_tr$tm_win_anm)
fl3m_ts$tm_spr_anm=(fl3_ts$tm_spr_anm - mean(fl3_tr$tm_spr_anm))/sd(fl3_tr$tm_spr_anm)
fl3m_ts$ppt_spr_anm=(fl3_ts$ppt_spr_anm - mean(fl3_tr$ppt_spr_anm))/sd(fl3_tr$ppt_spr_anm)
fl3m_ts$pop_den=(log(fl3_ts$pop_den) - mean(log(fl3_tr$pop_den)))/sd(log(fl3_tr$pop_den))


cor(fl3_tr[, c(1, 8:11, 19:24)])



# frost risk model
frisk.m1 <- lmer(frostrisk ~ LT_tm_spr + LT_ppt_spr +
                   LT_tm_spr*pop_den+ LT_ppt_spr:pop_den+
                   tm_spr_anm*pop_den + tm_win_anm*ppt_spr_anm +ppt_spr_anm*pop_den+
                   (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2),
                 data = fl3m_tr, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(frisk.m1)

frisk.m2 <- lmer(frostrisk ~ LT_tm_spr + LT_ppt_spr +
                   LT_tm_spr*pop_den+ LT_ppt_spr:pop_den+
                   tm_spr_anm*pop_den + tm_win_anm:ppt_spr_anm +ppt_spr_anm*pop_den+
                   (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2),
                 data = fl3m_tr, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(frisk.m2)

frisk.m3 <- lmer(frostrisk ~ LT_tm_spr + LT_ppt_spr +
                   LT_tm_spr*pop_den+ LT_ppt_spr:pop_den+
                   tm_spr_anm + tm_win_anm:ppt_spr_anm +ppt_spr_anm*pop_den+
                   (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2),
                 data = fl3m_tr, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(frisk.m3)

frisk.m4 <- lmer(frostrisk ~ LT_tm_spr + LT_ppt_spr +
                   LT_tm_spr:pop_den+ LT_ppt_spr:pop_den+
                   tm_spr_anm + tm_win_anm:ppt_spr_anm +ppt_spr_anm:pop_den+ppt_spr_anm+
                   (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2),
                 data = fl3m_tr, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(frisk.m4)

frisk.m5 <- lmer(frostrisk ~ LT_tm_spr + LT_ppt_spr +
                   LT_tm_spr:pop_den+ LT_ppt_spr:pop_den+
                   tm_spr_anm + tm_win_anm:ppt_spr_anm +ppt_spr_anm+
                   (1+tm_spr_anm+ tm_win_anm|binomial_species)+ (1|sc2),
                 data = fl3m_tr, REML=F, control = lmerControl(optimizer="bobyqa"))
summary(frisk.m5)

AICc(frisk.m1, frisk.m2, frisk.m3, frisk.m4, frisk.m5)

frisk.mod <- lmer(frostrisk ~ LT_tm_spr + LT_ppt_spr +
                    LT_tm_spr:pop_den+ LT_ppt_spr:pop_den+
                    tm_spr_anm + tm_win_anm:ppt_spr_anm +
                    (1+tm_spr_anm+ tm_win_anm|binomial_species),
                 data = fl3m_tr, REML=T, control = lmerControl(optimizer="bobyqa"))

summary(frisk.mod)
coef(frisk.mod)
r.squaredGLMM(frisk.mod)

plot_model(frisk.mod, type="pred", terms=c("pop_den","LT_tm_spr"))


fl3m_ts$pred=predict(frisk.mod, newdata=fl3m_ts)

cor.test(fl3m_ts$frostrisk, fl3m_ts$pred)
sqrt(sum((fl3m_ts$pred-fl3m_ts$frostrisk)^2)/length(fl3m_ts$pred))  
mean(abs(fl3m_ts$pred-fl3m_ts$frostrisk)) 


# Figure S2
tiff(filename="figures_rv/frostrisk_model_validation.tif", width= 1500, height= 1500, res=300, compression= "lzw")
ggplot(fl3m_ts, aes(x=frostrisk, y=pred))+
  geom_point(pch=".")+
  #geom_smooth(method="lm", se=F)+
  geom_abline(slope = 1, intercept = 0, color="red", linetype=2, linewidth=1)+
  xlim(-350,100)+ylim(-350,100)+
  #annotate("text", x = 0, y = -300, label = bquote("r"^2~"= 0.67"))+
  annotate("text", x = 0, y = -300, label = "r=0.8")+
  annotate("text", x = 4, y = -320, label = "p<0.001")+
  labs(x="frost risk of flowering (days)", y="predicted frost risk (days)")+
  theme_cowplot()
dev.off()




# frost risk model effect plot

# Figure 5
line_df <- data.frame(x=fl3m_tr$pop_den, y=-85, xend = fl3m_tr$pop_den, yend = -83)

p1=plot_model(frisk.mod, type="pred", terms=c("pop_den", "LT_tm_spr[-0.9988, 0.0031, 1.005]"), ci.lvl = F,
              title = "", 
              colors = c("blue","purple","red")) +
  
  scale_color_manual(values=c("blue","purple","red"), labels=c(0.4, 7.2, 14.0),
                     name="long-term average \nof spr temp (\u00b0C)") +
  scale_x_continuous(breaks=c(  -2.72, -1.06, 0.6, 2.27, 3.93),
                     labels=c( 1, 10, 100, 1000, 10000), limits=c(-4.4, 6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-85, -30), expand = c(0, 0)) +
  geom_segment(data=line_df, aes(x=x, y=y, xend =xend, yend = yend),inherit.aes = FALSE)+
  
  labs(x=bquote("population "~"km"^-2), y="frost risk (days)")+
  theme_cowplot()+
  theme(legend.position=c(0.65, 0.85), legend.background=element_blank())



p2=plot_model(frisk.mod, type="pred", terms=c("pop_den", "LT_ppt_spr[-0.9959, 0.0035, 1.003]"), ci.lvl = F,
              title = "",
              colors = c("deepskyblue","dodgerblue2","blue")) +
  scale_color_manual(values=c("deepskyblue","dodgerblue2","blue"), labels=c(257, 297, 337),
                     name="long-term average \nof spring prcp (mm)") +
  scale_x_continuous(breaks=c(  -2.72, -1.06, 0.6, 2.27, 3.93),
                     labels=c( 1, 10, 100, 1000, 10000), limits=c(-4.4, 6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-85, -30), expand = c(0, 0)) +
  geom_segment(data=line_df, aes(x=x, y=y, xend =xend, yend = yend),inherit.aes = FALSE)+
  
  labs(x=bquote("population "~"km"^-2), y="frost risk (days)")+
  theme_cowplot()+
  theme(legend.position=c(0.60, 0.85), legend.background=element_blank())


plot.all=cowplot::plot_grid(p1, p2, labels=c("(a)","(b)"), label_size=14, ncol=2, nrow=1, rel_widths=c(1, 1))

ggsave(filename = "figures_rv/frostrisk_effectplot_interactions_rv.tif", plot.all, device="tiff",
       dpi=300, height= 5, width=10, compression="lzw")



  

#------------- prediction (1981-2010 normal) ---------------------------
clim=read_csv(file="data/climate_data_StateCountyYear_20221123.csv")
LT = clim %>% group_by(sc) %>% summarise(LT_MAT=mean(MAT, na.rm=T),
                                         LT_MAP=mean(MAP, na.rm=T),
                                         LT_tm_win=mean(tm_win, na.rm=T),
                                         LT_tm_spr=mean(tm_spr, na.rm=T),
                                         LT_ppt_spr=mean(ppt_spr, na.rm=T))
colnames(LT)[1]="sc2"


# prepare new data
# species and sc2
spp_sc2=fl3 %>% select(sc2, binomial_species, LT_tm_spr, LT_ppt_spr) %>% distinct()

# spatial
spv=fl3 %>% select(sc2) %>% distinct() 

# climate
clmv=clim %>% filter(year>=1981, year<=2010) %>%
  mutate(sc2=sc) %>% left_join(., LT) %>%
  mutate(tm_win_anm=tm_win-LT_tm_win,
         tm_spr_anm=tm_spr-LT_tm_spr,
         ppt_spr_anm=(ppt_spr-LT_ppt_spr)/LT_ppt_spr) %>%
  group_by(sc2) %>%
  summarise(twn=mean(tm_win_anm, na.rm=T), tsp=mean(tm_spr_anm, na.rm=T), psp=mean(ppt_spr_anm, na.rm=T),
            twn.sd=sd(tm_win_anm, na.rm=T), tsp.sd=sd(tm_spr_anm, na.rm=T), psp.sd=sd(ppt_spr_anm, na.rm=T))

# population density
pd=data.frame(sc2=rep(unique(fl3$sc2),each=100),pop_den=rep(seq(1,24000,240),length(unique(fl3$sc2))))


# new data - 30y normal scenario 
ndata=left_join(spp_sc2, pd[, c(1,5)]) 

ndata= clmv %>% select(sc2, twn, tsp, psp) %>%
  rename(tm_win_anm=twn, tm_spr_anm=tsp, ppt_spr_anm=psp) %>%
  right_join(ndata)


ndata_m=ndata %>% select(sc2, binomial_species)
ndata_m$LT_tm_spr=(ndata$LT_tm_spr - mean(fl3_tr$LT_tm_spr))/sd(fl3_tr$LT_tm_spr)
ndata_m$LT_ppt_spr=(ndata$LT_ppt_spr - mean(fl3_tr$LT_ppt_spr))/sd(fl3_tr$LT_ppt_spr)
ndata_m$tm_win_anm=(ndata$tm_win_anm - mean(fl3_tr$tm_win_anm))/sd(fl3_tr$tm_win_anm)
ndata_m$tm_spr_anm=(ndata$tm_spr_anm - mean(fl3_tr$tm_spr_anm))/sd(fl3_tr$tm_spr_anm)
ndata_m$ppt_spr_anm=(ndata$ppt_spr_anm - mean(fl3_tr$ppt_spr_anm))/sd(fl3_tr$ppt_spr_anm)
ndata_m$pop_den=(log(ndata$pop_den) - mean(log(fl3_tr$pop_den)))/sd(log(fl3_tr$pop_den))


# predict flowering date and LFD
ndata_m$frostrisk=predict(frisk.mod, newdata=ndata_m)

summary(ndata_m$frostrisk)

summary(ndata_m)

save(ndata, ndata_m, file="results/frostrisk_prediction_rv.Rdata")
load("results/frostrisk_prediction_rv.Rdata")



# calculate slope of frostrisk~log(pop_den) for each species at each county
frisktrend <- function(data)
{
  lm1sum <- summary(lm(frostrisk~pop_den, data=data))
  slp1=lm1sum$coefficients[2,1]
  return(data.frame(sc2=unique(data$sc2), slp.norm=slp1))
}

frisktrend1 = ndata_m %>% na.omit() %>% 
  group_by(sc2) %>% do(frisktrend(.))

summary(frisktrend1)

write.csv(frisktrend1, gzfile("results/predition_frisktrend_rv2.gz"))
frisktrend1=read.csv("results/predition_frisktrend_rv2.gz")
frisktrend1=frisktrend1[, -1]


# plot counties of increased or reduced frost risk (Figure 6)
county_area=read_csv(file="data/DC.EC_county_area_updated.csv")
county_area$sc2=substr(county_area$scd, 1, nchar(county_area$scd)-5)
head(county_area)

county_area=distinct(county_area)

frisktrend1= left_join(distinct(frisktrend1), distinct(county_area[,10:11]))

county=st_read(dsn ="C:/Data/shapefiles/cb_2018_us_county_5m/cb_2018_us_county_5m.shp")

county1=left_join(county, frisktrend1)



# Figure 6
p=ggplot(county1, aes(fill=slp.norm)) +
  geom_sf()+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name="change of\nfrost risk") +
  coord_sf(ylim=c(24,48), xlim=c(-87, -66))+
  theme_bw(14) +
  theme(legend.position = c(0.72, 0.22))


county=left_join(county, county_area[, 10:11]) %>% left_join(., LT)

p1=ggplot(county, aes(fill=LT_tm_spr, color=LT_tm_spr)) +
  geom_sf()+
  scale_fill_viridis(name="Long-term\nspring temp (\u00b0C)") +
  scale_color_viridis(guide="none") +
  coord_sf(ylim=c(24,48), xlim=c(-87, -66))+
  theme_bw(14) +
  theme(legend.position = c(0.72, 0.22))

p2=ggplot(county, aes(fill=LT_ppt_spr, color=LT_ppt_spr)) +
  geom_sf()+
  scale_fill_viridis(name="Long-term\nspring prcp (mm)") +
  scale_color_viridis(guide="none") +
  coord_sf(ylim=c(24,48), xlim=c(-87, -66))+
  theme_bw(14) +
  theme(legend.position = c(0.72, 0.22))

maps=cowplot::plot_grid(p, p1, p2, labels=c("(a)","(b)","(c)"), label_size=14, ncol=3, nrow=1, 
                        rel_widths=c(1, 1, 1))

ggsave(filename = "figures_rv/frostrisk_LTclimate.tif", maps, device="tiff",
       dpi=300, height= 5, width=12, compression="lzw")





