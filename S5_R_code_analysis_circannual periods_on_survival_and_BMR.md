###Appendix S7
## Testing for the correlation of free-running circannual period with survival using MCMCglmm and mulTree

Supplementary online material  to Julia Karagicheva, Eldar Rakhimberdiev, Anatoly Saveliev & Theunis Piersma 2017 _Circannual rhythms functionally link life histories and life cycles in birds._ - Journal of Functional Ecology 000: 000-000.

Here we provide a code for the analysis presented in the paper, assessing the relationships between deviations of circannual periods from 365 days (circannual deviations) and species-specific annual survival rates. 

In the analyses we used data from published studies on 14 passerine species and 2 shorebird species (Circannual_cycles_birds_data.csv)

Then we extended each glmmMCMC model by incorporating phylogenetic variance in it, using the package mulTree. From the output, we calculated phylogenetic signal and an analogue of Pagel lambda.

Install required packages

```{r, eval = F}

install.packages('snow')
if(!require(devtools)) install.packages("devtools")
library(devtools)
install_github("TGuillerme/mulTree", ref = "release")

```
you might need to install some packages manually

```{r, eval = F}

install.packages(c("MCMCglmm", "coda", "hdrcde", "snow", "ape", "corpcor", "curl"))

```
Load required packages

```{r, eval = F}

library(snow)
library(mulTree)

```
### PREPARE THE DATA 

Download data from GitHub

load table with circannual cycles and deviations of circannual period from 365 days (main table "Circannual_cycles_birds_data")

```{r, eval = F}

Circannual_cycles_birds_data<-read.csv('https://git.io/v5GIF',stringsAsFactors=F)

```
create unique ID 

```{r, eval = F}

Circannual_cycles_birds_data$uniqueid<-paste(Circannual_cycles_birds_data$species,
     Circannual_cycles_birds_data$Reference,
     Circannual_cycles_birds_data$BirdID_char,sep="_")

```
make cycle number character

```{r, eval = F}

Circannual_cycles_birds_data$cyclenr_char<-
     as.character(Circannual_cycles_birds_data$cyclenr_char)

```
load table "Survival_birds_data" with species-specific survival and calculate median per species

```{r, eval = F}

Survival_birds_data<-read.csv('https://git.io/v5Gtq',stringsAsFactors=F)
surv_sp<-aggregate(Survival_birds_data$survival,
      by=list(Survival_birds_data$species),median)

```
add survival values to the main table

```{r, eval = F}

Circannual_cycles_birds_data$surv<-NA
for(i in 1:nrow(Circannual_cycles_birds_data)){
Circannual_cycles_birds_data$surv[i]<-
     unique(surv_sp[surv_sp$Group.1 %in%
     Circannual_cycles_birds_data$species[i],]$x)
}

```
load table "BMR_birds_data.csv" with species-specific BMR  and calculate median per species

```{r, eval = F}

BMR_birds_data<-read.csv('https://git.io/v5ZSb',stringsAsFactors=F)
bmr_sp_median<-aggregate(BMR_birds_data$BMR_W,
                 by=list(BMR_birds_data$species),median)

```
add BMR values to the main table

```{r, eval = F}
Circannual_cycles_birds_data$bmr<-NA
for(i in 1:nrow(Circannual_cycles_birds_data)){
Circannual_cycles_birds_data$bmr[i]<-
     unique(bmr_sp_median[bmr_sp_median$Group.1==
            Circannual_cycles_birds_data$species[i],]$x)
}

```

Calculate residuals of survival on BMR

```{r, eval = F}

Circannual_cycles_birds_data$resbmr<-
residuals(lm(log(surv)~log(bmr),data=Circannual_cycles_birds_data))

```
Calculate residuals of circannual deviation on BMR

```{r, eval = F}

Circannual_cycles_birds_data$res_dev_bmr<-
     residuals(lm(deviation~log(bmr),data=Circannual_cycles_birds_data))

```

### RUN THE MODELS

### Set priors

```{r, eval = F}

p.var<-var(Circannual_cycles_birds_data$res_dev_bmr,na.rm=TRUE)
pri5_1 <- list(R = list(V = matrix(p.var/2), nu = 1), 
     G = list(G1 = list(V = matrix(p.var/2), nu = 1),
     G2 = list(V = matrix(p.var/2), nu = 1), 
     G3=list(V = matrix(p.var/2), nu = 1),
     G4=list(V = matrix(p.var/2), nu = 1),
     G5=list(V = matrix(p.var/2), nu = 1)))
```
### Get phylogenetic data

download phylogenetic trees (containing shorebirds), directly from http://birdtree.org/ or just download the 100 trees used in the presented analysis from https://git.io/v5ZxG :

```{r, eval = F}

download.file("https://git.io/v5ZxG", "with_shorebirds_tree_subset.RData",mode="wb")
# load the phylogenetic trees subset
load("with_shorebirds_tree_subset.RData")
tr_full<-with_shorebirds_tree_subset

# select one phylogenetic tree for model selection in MCMCglmm
mulTree_data <- as.mulTree(data = trait_data_with_shorebirds, tree = tr_full,    taxa = "animal",  rand.term=~uniqueid+Trait+Reference+species+animal)

ped<-mulTree_data$phy[[1]]
```

### run MCMCglmm 

```{r, eval = F}
model1_pedigree <- MCMCglmm(res_dev_bmr~resbmr:cyclenr_char+cyclenr_char,random= ~uniqueid+Trait+Reference+species+animal,data= trait_data_with_shorebirds, pedigree=ped,prior=pri5_1, nitt=130000*10,thin=10*5,burnin=3000*10,pr=TRUE,verbose=F,nodes="ALL")

model2_pedigree <- MCMCglmm(res_dev_bmr~resbmr+cyclenr_char,random= ~uniqueid+Trait+Reference+species+animal,data= trait_data_with_shorebirds, pedigree=ped,prior=pri5_1, nitt=130000*10,thin=10*5,burnin=3000*10,pr=TRUE,verbose=F,nodes="ALL")

model1_pedigree$DIC
# 4000.948
model2_pedigree$DIC
# 4009.138
```
The best model

```{r, eval = F}

summary(model1_pedigree)

```
control for autocorrelation

```{r, eval = F}

autocorr(model1_pedigree$VCV)

```

### Run the best model with phylogeny using mulTree

set the directory, where you wish to store the outputs

### set the model structure and parameters

```{r, eval = F}

mulTree_data <- 
as.mulTree(data = Circannual_cycles_birds_data, tree = tr_full,
     taxa = "animal", rand.term=~uniqueid+Trait+Reference+animal+species)
formula_fintest<-res_dev_bmr~resbmr:cyclenr_char+cyclenr_char
mulTree.parameters<-c(5000000, 5000, 100000)

```
### run mulTree

```{r, eval = F}

surv_w_shore<-mulTree(mulTree.data = mulTree_data, 
     formula = formula_fintest, priors = pri5_1, pr=TRUE, 
     parameters = mulTree.parameters, output = "surv_w_shore",
     chains=8, parallel='SOCK')

```

### EXTRACT OUTPUTS 

it is very likely that you will have to do it in a new R window. Therefore, again:

```{r, eval = F}
# load the packages
library(snow)
library(mulTree)

#... and set the work directory, where the model output is stored
```
### Check whether all the chains have converged

this function gets convergence parameters from 'conv' files, automatically saved in the same directory as model outputs

```{r, eval = F}

extract_conv<-function(model_name) {
  Files<-list.files( pattern=paste0(model_name,'.*_conv'))
  cat(length(Files), 'files found\n')
  Res<-c()
  for(i in 1:length(Files)){
  mod<-load(Files[i])
  treename<-Files[i]
  res<-unlist(converge.test)<1.1
Res1<-c(treename,res)
Res<-rbind(Res,Res1)
  }
  return(Res)
  }
  
  ```
### control for convergence

```{r, eval = F}

conv_mother<-extract_conv('surv_w_shore')

```
## Collect mulTree outputs

you can first load one of the models' output, to see its structure

we need p-values (pMCMC) for the fixed parameters: 
pMCMC_cycle: intercept between the transitory (baseline) and full cycle (cyclenr_char1)
pMCMC_surv_cycle_0: slope for the transitory cycle
pMCMC_surv_cycle_1: slope for the full cycle

### function extracting pMCMC

```{r, eval = F}

extract_pMCMC<-function(model_name) {
  Files<-list.files( pattern=paste0(model_name,'.*_chain'))
  cat(length(Files), 'files found\n')
	pMCMC_cycle<-c()
	pMCMC_slope_cycle_0<-c()
	pMCMC_slope_cycle_1<-c()

     for(i in 1:length(Files)){
  mod<-load(Files[i])
  pMCMC_cycle <-c(pMCMC_cycle, summary(model)$solutions[18])
  pMCMC_slope_cycle_0 <-c(pMCMC_slope_cycle_0, summary(model)$solutions[19])
  pMCMC_slope_cycle_1 <-c(pMCMC_slope_cycle_1, summary(model)$solutions[20])

  }
  Res<-list(pMCMC_cycle=pMCMC_cycle,pMCMC_slope_cycle_0=pMCMC_slope_cycle_0,
       pMCMC_slope_cycle_1=pMCMC_slope_cycle_1)
  return(Res)
}

```
### extract pMCMC for the model

```{r, eval = F}

pMCMC<-extract_pMCMC('surv_w_shore')
summary(pMCMC$pMCMC_cycle)
summary(pMCMC$pMCMC_slope_cycle_0)
summary(pMCMC$pMCMC_slope_cycle_1)

```
### function extracting fixed effects

```{r, eval = F}

extract_fixed<-function(model_name) {
  Files<-list.files( pattern=paste0(model_name,'.*_chain'))
  cat(length(Files), 'files found\n')
	Intercept<-c()
	cycle_nr_1<-c()
	slope_cycle_0<-c()
	slope_cycle_1<-c()
  
     for(i in 1:length(Files)){
  mod<-load(Files[i])
  Intercept <-c(Intercept, model$Sol[,1])
  cycle_nr_1 <-c(cycle_nr_1, model$Sol[,2])
  slope_cycle_0 <-c(slope_cycle_0, model$Sol[,3])
  slope_cycle_1 <-c(slope_cycle_1, model$Sol[,4])
  }
  Res<-list(Intercept=Intercept, cycle_nr_1=cycle_nr_1,
       slope_cycle_0=slope_cycle_0,slope_cycle_1=slope_cycle_1)
  return(Res)
}

```
### extract fixed effects

```{r, eval = F}

fixed_effects<-extract_fixed('surv_w_shore')
summary(fixed_effects$Intercept)
summary(fixed_effects$cycle_nr_1)
summary(fixed_effects$slope_cycle_0)
summary(fixed_effects$slope_cycle_1)

quantile(fixed_effects$Intercept,c(0.025,0.975))
quantile(fixed_effects$cycle_nr_1,c(0.025,0.975))
quantile(fixed_effects$slope_cycle_0,c(0.025,0.975))
quantile(fixed_effects$slope_cycle_1,c(0.025,0.975))

```

### function extracting random effects

```{r, eval = F}

extract_random<-function(model_name) {
  Files<-list.files( pattern=paste0(model_name,'.*_chain'))
  cat(length(Files), 'files found\n')
  uniqueid_post<-c()
  Trait_post<-c()
  Reference_post<-c()
  animal_post<-c()
  species_post<-c()
  units_post<-c()
     for(i in 1:length(Files)){
  mod<-load(Files[i])
  uniqueid_post <-c(uniqueid_post, model$VCV[,1])
  Trait_post <-c(Trait_post, model$VCV[,2])
  Reference_post <-c(Reference_post, model$VCV[,3])
  animal_post <-c(animal_post, model$VCV[,4])
  species_post <-c(species_post, model$VCV[,5])
  units_post <-c(units_post, model$VCV[,6])
   
  }
  Res<-list(uniqueid_post=uniqueid_post, Trait_post=Trait_post,
       Reference_post=Reference_post,animal_post=animal_post, 
       species_post=species_post,units_post=units_post)
  return(Res)
}

```
### extract random effects

```{r, eval = F}

random_effects<-extract_random('surv_w_shore')
summary(random_effects$uniqueid_post)
summary(random_effects$Trait_post)
summary(random_effects$Reference_post)
summary(random_effects$animal_post)
summary(random_effects$species_post)
summary(random_effects$units_post)

quantile(random_effects$uniqueid_post,c(0.025,0.975))
quantile(random_effects$Trait_post,c(0.025,0.975))
quantile(random_effects$Reference_post,c(0.025,0.975))
quantile(random_effects$animal_post,c(0.025,0.975))
quantile(random_effects$species_post,c(0.025,0.975))
quantile(random_effects$units_post,c(0.025,0.975))

```

### CALCULATE PHYLOGENETIC SIGNAL AND REPEATABILITIES

```{r, eval = F}


# phylogenetic signal and intra-specific repeatability

phyl_sig<-var(random_effects$animal_post)/(var(random_effects$uniqueid_post)+var(random_effects$Trait_post)+var(random_effects$Reference_post)+var(random_effects$animal_post)+var(random_effects$species_post)+var(random_effects$units_post))

# phyl_sig = 0.423506

lambda<-var(random_effects$animal_post)/(var(random_effects$animal_post)+var(random_effects$species_post))

#lambda = 0.812362

#Intra-specific repeatability

repeatab<-var(random_effects$species_post)/(var(random_effects$uniqueid_post)+var(random_effects$Trait_post)+var(random_effects$Reference_post)+
var(random_effects$animal_post)+var(random_effects$species_post)+var(random_effects$units_post))

# repeatab = 0.09782072

repeatab_adjusted<-var(random_effects$species_post)/(var(random_effects$species_post)+var(random_effects$units_post))

# repeatab_adjusted = 0.9517922

```
