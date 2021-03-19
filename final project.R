library(SUMMER)
library(survey)
library(tidyverse)
library(haven)
library(rgdal)
library(foreign)
library(stringr)
library(mapproj)
library(cowplot)
library(lme4)
library(corrplot)
library(table1)

##county-level (child$sregion)

birth<-read.dta("data/KEBR72FL.DTA")
household<-read.dta("data/KEHR72FL.DTA")
individual<-read.dta("data/KEIR72FL.DTA")

child<-read.dta("data/KEKR72FL.DTA")
summary(child$m14)
# Antenatal care visits (m14) are presented in groups of none, 1-3, 4+.
# binary ANC4
child$ANC4<-with(child,ifelse(is.na(m14),NA,
                              ifelse(m14<4,0,1)))


#data merge

## Included demographic variable
# maternal education(v106), birth order(bord), household wealth ((hv270, v190, mv190)), 
# household residence type ( (hv025, v025, v102, mv025)), marital status(v502), ethnicity(v131), 
# parity(v219,v220),age at first marriage/cohabitation (v511), place of delivery (m15), 
# sex of household head (v151), religion(v130), maternal age (int((b3 - v011)/12))

individual.sub<-individual%>%select(v106,v502,v131,v130,v220,v511,v011,v151,caseid,sregion)
individual.sub$v131<-with(individual.sub,ifelse(v131 %in% c("kalenjin","kikuyu","kamba","kisii","luhya","luo"),v131,"other tribes"))

birth.sub<-birth%>%filter(bidx==1)%>%select(bord,b3,m15,caseid)
birth.sub$m15<-with(birth.sub,ifelse(m15 %in% c("government hospital","govt health center",
                                                "mission hospital/clinic","private hospital/clinic"),"health facility","non-health facility"))

child.sub<-child%>%filter(midx==1)%>%select(ANC4,caseid)
child.sub$hhid<-str_sub(child.sub$caseid,1,-4)
household.sub<-household%>%select(hv270,hv025,hhid)
merged<- individual.sub %>% inner_join(birth.sub, by = "caseid")%>%
  inner_join(child.sub,by="caseid")%>%
  inner_join(household.sub,by="hhid")
merged$maternal_age<-round((merged$b3 - merged$v011)/12,0)
merged<-merged%>%select(-c("b3","v011"))


colnames(merged)<-c("maternal education","marital status","ethnicity","religion","parity", 
                    "age of first marriage","sex of household head","caseid","sregion",
                    "birth order","place of delivery", "ANC4","hhid","household wealth",
                    "household residence type","maternal age")

merged$`age of first marriage`<-with(merged,ifelse(`age of first marriage`>=18,">=18","<18"))
merged$`maternal age`<-with(merged,ifelse(`maternal age`>=18,">=18","<18")) 
merged$parity<-factor(merged$parity)

# Table 1
table1(~.,data=merged[,-c(3,8,9,13)])

# demographic variable and ANC4
merged$parity<-as.numeric(merged$parity)
model_maternal_educ<-glmer(ANC4~`maternal education`+(1|sregion),data=merged,family=binomial)
model_marital_statusc<-glmer(ANC4~`marital status`+(1|sregion),data=merged,family=binomial)
model_eth<-glmer(ANC4~ethnicity+(1|sregion),data=merged,family=binomial)
model_religion<-glmer(ANC4~religion+(1|sregion),data=merged,family=binomial)
model_parity<-glmer(ANC4~parity+(1|sregion),data=merged,family=binomial)
model_age1m<-glmer(ANC4~`age of first marriage`+(1|sregion),data=merged,family=binomial)
model_sex_hh<-glmer(ANC4~`sex of household head`+(1|sregion),data=merged,family=binomial)
model_border<-glmer(ANC4~`birth order`+(1|sregion),data=merged,family=binomial)
model_place_deliver<-glmer(ANC4~`place of delivery`+(1|sregion),data=merged,family=binomial)
model_hh_wealth<-glmer(ANC4~`household wealth`+(1|sregion),data=merged,family=binomial)
model_hh_redi<-glmer(ANC4~`household residence type`+(1|sregion),data=merged,family=binomial)
model_maternal_age<-glmer(ANC4~`maternal age`+(1|sregion),data=merged,family=binomial)
model_null<-glmer(ANC4~1+(1|sregion),data=merged,family=binomial)

anova(model_hh_redi,model_null)

rownames(coef)<-NULL

coef<-coef(summary(model_hh_redi))
coef[,1]-1.96*coef[,2]
coef[,1]+1.96*coef[,2]
coef[,1]

#The individual weight for women (v005) is the household weight (hv005) multiplied by the inverse of the individual response rate for women in the stratum.
# Primary sampling unit v021 (PSU/EAS): smallest unit with correspondent GPS coordinate
child<-na.omit(subset(child,select=c("ANC4","v021","v023","v005","sregion")))

geo_county <- rgdal::readOGR(dsn="shps_subcounty",layer="ken_admbnda_adm1_iebc_20191031")
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))

ggplot(data = world) +
  geom_sf() +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "darkblue", fontface = "bold", check_overlap = FALSE) +
  annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", 
           fontface = "italic", color = "grey22", size = 6) +
  coord_sf(xlim = c(-102.15, -74.12), ylim = c(7.65, 33.97), expand = FALSE)


#harmonize the coding of county name between shapefile and dhs data sregion
geo_county$county_name[47]<-str_sort(unique(child$sregion))[31]
geo_county$county_name[-47]<-str_sort(unique(child$sregion))[-31]
# generates the spatial adjacency matrix Amat using the function getAmat()
Amat<-getAmat(geo,geo$REGNAME)
Amat_county<-getAmat(geo_county,geo_county$county_name)

DHSdesign<-svydesign(id=child$v021, strata=child$v023, weights=child$v005, data=child)
direct <- svyby(~ANC4, ~sregion, DHSdesign, svymean)


# polygon information is input, and the neighbors in the Amat argument - this is required for the ICAR.
# spatial model, acknowledging the sampling, design
smoothed <- fitGeneric(data = child, geo = geo_county,
                       Amat = Amat_county, responseType = "binary", responseVar = "ANC4", strataVar = "v023", weightVar = "v005", regionVar = "sregion", clusterVar = "~v021", CI = 0.95)

## map the posterior mean estimates
toplot <- smoothed$smooth
mapPlot(data = toplot, geo = geo_county, variables = c("mean"),
        labels = c("Posterior Mean"), by.data = "region", by.geo = "county_name")

## map 2.5% and 97.5% Quantile estimate
mapPlot(data = toplot, geo = geo_county, variables = c("lower", "upper"), labels = c("Lower", "Upper"), by.data = "region", by.geo = "county_name")

## map the interval width
toplot$width <- toplot$upper - toplot$lower
mapPlot(data = toplot, geo = geo_county, variables = c("width"),
        labels = c("Width"), by.data = "region", by.geo = "county_name")

data<-smoothed$smooth
 data[data$region=="mandera",]
 data[data$region=="west pokot",]
 data[data$region=="nairobi",]
 
 
 # map direct estimate
 toplot$HTest <- smoothed$HT$HT.est
 mapPlot(data = toplot, geo = geo_county, variables = c("HTest"),
         labels = c("Direct Estimates"), by.data = "region", by.geo = "county_name")
 # map lower and upper endpoints of 95% CI for direct estimate
 lo <- smoothed$HT$HT.est - 1.96 * sqrt(smoothed$HT$HT.var)
 hi <- smoothed$HT$HT.est + 1.96 * sqrt(smoothed$HT$HT.var) 
 toplot$HTlower <- lo
 toplot$HTupper <- hi
 mapPlot(data = toplot, geo = geo_county, variables = c("HTlower", "HTupper"), 
         labels = c("Direct Lower", "Direct Upper"), by.data = "region", by.geo = "county_name")
 
#compare the estimate variance between design-based direct estimate and smoothed estimate
est <- data.frame(direct = smoothed$HT$HT.est, 
                  smooth = smoothed$smooth$mean)
var <- data.frame(direct = smoothed$HT$HT.var, 
                  smooth = smoothed$smooth$var)
l1 <- range(est)
l2 <- range(var)
g1 <- ggplot(est, aes(x = direct, y = smooth)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  ggtitle("Estimates") + xlab("Direct estimate") + 
  ylab("Smoothed estimate") + xlim(l1) + ylim(l1)+
  theme(text=element_text(size=20)) 
g2 <- ggplot(var, aes(x = direct, y = smooth)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  ggtitle("Variance") + xlab("Direct estimate") + 
  ylab("Smoothed estimate") + xlim(l2) + ylim(l2)+
  theme(text=element_text(size=20)) 
plot_grid(g1,g2)

# adding area level covariate
# Covariates at the area level could be included in the smoothing model by specifying the X argument.
# The first column of the covariate matrix should be the names of regions that match the column and row names of the adjacency matrix.

## posterior estimate and interval for each county, comparing with direct estimation
data.plot<-data.frame(est=c(toplot$mean,toplot$HTest),
                      lower=c(toplot$lower,toplot$HTlower),
                      upper=c(toplot$upper,toplot$HTupper),
                      region=rep(toplot$region,by=2),
                      method=rep(c("smooth","direct"),each=47))

data.plot%>%
  mutate(region= fct_reorder(region, est))%>%
  ggplot(aes(x=region,y=est,color=method)) + 
  geom_point()+
  geom_errorbar(aes(x = region, ymin = lower, ymax = upper))+ylab("ANC4 estimate")+
  coord_flip()













