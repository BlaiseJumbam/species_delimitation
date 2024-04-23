source("./code/packages_for_gymc.R")

# import trees

ef1a_yule_tr <- read_nexus_phylo("data/ef1a_0723_yule3_new.nex")
ef1a_coal_tr <- read_nexus_phylo("data/ef1a_0723_constant_coalescent.nex")

elo2_yule_tr <- read_nexus_phylo("data/elo2_0723_yule.nex")
elo2_coal_tr <- read_nexus_phylo("data/elo2_0723_constant_coalescent.nex")

its_yule_tr <- read_nexus_phylo("data/its_0723_yule.nex")
its_coal_tr <- read_nexus_phylo("data/its_0723_constant_coalescent.nex")

lsu_yule_tr <- read_nexus_phylo("data/lsu_0723_yule.nex")
lsu_coal_tr <- read_nexus_phylo("data/lsu_0723_constant_coalescent.nex")

mtlsu_yule_tr <- read_nexus_phylo("data/mtlsu_0723_yule.nex")
mtlsu_coal_tr <- read_nexus_phylo("data/mtlsu_0723_constant_coalescent.nex")

mtssu_yule_tr <- read_nexus_phylo("data/mtssu_0723_yule.nex")
mtssu_coal_tr <- read_nexus_phylo("data/mtssu_0723_constant_coalescent.nex")

rpb2_yule_tr <- read_nexus_phylo("data/rpb2_0723_yule.nex")
rpb2_coal_tr <- read_nexus_phylo("data/rpb2_0723_constant_coalescent.nex")

tub_yule_tr <- read_nexus_phylo("data/tub_0723_yule.nex")
tub_coal_tr <- read_nexus_phylo("data/tub_0723_constant_coalescent.nex")

# run gymc

ef1a_yule_gmyc <- gmyc(ef1a_yule_tr,method="multiple",interval=c(0,5),quiet=F)
ef1a_coal_gmyc <- gmyc(ef1a_coal_tr,method="multiple",interval=c(0,5),quiet=F)

elo2_yule_gmyc <- gmyc(elo2_yule_tr,method="multiple",interval=c(0,5),quiet=F)
elo2_coal_gmyc <- gmyc(elo2_coal_tr,method="multiple",interval=c(0,5),quiet=F)

its_yule_gmyc <- gmyc(its_yule_tr,method="multiple",interval=c(0,5),quiet=F)
its_coal_gmyc <- gmyc(its_coal_tr,method="multiple",interval=c(0,5),quiet=F)

lsu_yule_gmyc <- gmyc(lsu_yule_tr,method="multiple",interval=c(0,5),quiet=F)
lsu_coal_gmyc <- gmyc(lsu_coal_tr,method="multiple",interval=c(0,5),quiet=F)

mtlsu_yule_gmyc <- gmyc(mtlsu_yule_tr,method="multiple",interval=c(0,5),quiet=F)
mtlsu_coal_gmyc <- gmyc(mtlsu_coal_tr,method="multiple",interval=c(0,5),quiet=F)

mtssu_yule_gmyc <- gmyc(mtssu_yule_tr,method="multiple",interval=c(0,5),quiet=F)
mtssu_coal_gmyc <- gmyc(mtssu_coal_tr,method="multiple",interval=c(0,5),quiet=F)

rpb2_yule_gmyc <- gmyc(rpb2_yule_tr,method="multiple",interval=c(0,5),quiet=F)
rpb2_coal_gmyc <- gmyc(rpb2_coal_tr,method="multiple",interval=c(0,5),quiet=F)

tub_yule_gmyc <- gmyc(tub_yule_tr,method="multiple",interval=c(0,5),quiet=F)
tub_coal_gmyc <- gmyc(tub_coal_tr,method="multiple",interval=c(0,5),quiet=F)


# how many species are found 

summary(ef1a_yule_gmyc)
summary(ef1a_coal_gmyc)

summary(elo2_yule_gmyc)
summary(elo2_coal_gmyc)

summary(its_yule_gmyc)
summary(its_coal_gmyc)

summary(lsu_yule_gmyc)
summary(lsu_coal_gmyc)

summary(mtlsu_yule_gmyc)
summary(mtlsu_coal_gmyc)

summary(mtssu_yule_gmyc)
summary(mtssu_coal_gmyc)

summary(rpb2_yule_gmyc)
summary(rpb2_coal_gmyc)

summary(tub_yule_gmyc)
summary(tub_coal_gmyc)

# Plotting
# plots show 
# 1. # of lineages through time, with a red vertical line showing the 
# inferred position of the threshold
# 2. profile of the likelihood through time
# 3. tree with the individual clusters highlighted in red

plot(ef1a_yule_gmyc) 
plot(ef1a_coal_gmyc) 

plot(elo2_yule_gmyc) 
plot(elo2_coal_gmyc) 

plot(its_yule_gmyc) 
plot(its_coal_gmyc)

plot(lsu_yule_gmyc) 
plot(lsu_coal_gmyc)

plot(mtlsu_yule_gmyc) 
plot(mtlsu_coal_gmyc)

plot(mtssu_yule_gmyc) 
plot(mtssu_coal_gmyc)

plot(rpb2_yule_gmyc) 
plot(rpb2_coal_gmyc)

plot(tub_yule_gmyc) 
plot(tub_coal_gmyc)

# To know which samples are assigned to which species, use

spec.list(ef1a_yule_gmyc) %>% 
  write_csv(., file="results/ef1a_yule_gymc_multiple.csv")
spec.list(ef1a_coal_gmyc) %>% 
  write_csv(., file="results/ef1a_coal_gymc_multiple.csv")

spec.list(elo2_yule_gmyc)%>% 
  write_csv(., file="results/elo2_yule_gymc_multiple.csv")
spec.list(elo2_coal_gmyc)%>% 
  write_csv(., file="results/elo2_coal_gymc_multiple.csv")

spec.list(its_yule_gmyc)%>% 
  write_csv(., file="results/its_yule_gymc_multiple.csv")
spec.list(its_coal_gmyc)%>% 
  write_csv(., file="results/its_coal_gymc_multiple.csv")

spec.list(lsu_yule_gmyc)%>% 
  write_csv(., file="results/lsu_yule_gymc_multiple.csv")
spec.list(lsu_coal_gmyc)%>% 
  write_csv(., file="results/lsu_coal_gymc_multiple.csv")

spec.list(mtlsu_yule_gmyc)%>% 
  write_csv(., file="results/mtlus_yule_gymc_multiple.csv")
spec.list(mtlsu_coal_gmyc)%>% 
  write_csv(., file="results/mtlsu_coal_gymc_multiple.csv")

spec.list(mtssu_yule_gmyc)%>% 
  write_csv(., file="results/mtssu_yule_gymc_multiple.csv")
spec.list(mtssu_coal_gmyc)%>% 
  write_csv(., file="results/mtssu_coal_gymc_multiple.csv")

spec.list(rpb2_yule_gmyc)%>% 
  write_csv(., file="results/rpb2_yule_gymc_multiple.csv")
spec.list(rpb2_coal_gmyc)%>% 
  write_csv(., file="results/rpb2_coal_gymc_multiple.csv")

spec.list(tub_yule_gmyc)%>% 
  write_csv(., file="results/tub_yule_gymc_multiple.csv")
spec.list(tub_coal_gmyc)%>% 
  write_csv(., file="results/tub_coal_gymc_multiple.csv")

# plot the “support” for the delineated species. It can give an indication on 
# whether you can trust the results or not.

# estimate support
ef1a_yule_support <- gmyc.support(ef1a_yule_gmyc) 
ef1a_coal_support <- gmyc.support(ef1a_coal_gmyc) 

elo2_yule_support <- gmyc.support(elo2_yule_gmyc) 
elo2_coal_support <- gmyc.support(elo2_coal_gmyc)

its_yule_support <- gmyc.support(its_yule_gmyc) 
its_coal_support <- gmyc.support(its_coal_gmyc)

lsu_yule_support <- gmyc.support(lsu_yule_gmyc) 
lsu_coal_support <- gmyc.support(lsu_coal_gmyc)

mtlsu_yule_support <- gmyc.support(mtlsu_yule_gmyc) 
mtlsu_coal_support <- gmyc.support(mtlsu_coal_gmyc)

mtssu_yule_support <- gmyc.support(mtssu_yule_gmyc) 
mtssu_coal_support <- gmyc.support(mtssu_coal_gmyc)

rpb2_yule_support <- gmyc.support(rpb2_yule_gmyc) 
rpb2_coal_support <- gmyc.support(rpb2_coal_gmyc)

tub_yule_support <- gmyc.support(tub_yule_gmyc) 
tub_coal_support <- gmyc.support(tub_coal_gmyc)

# show values for affected nodes
is.na(ef1a_yule_support[ef1a_yule_support == 0]) <- TRUE
is.na(ef1a_coal_support[ef1a_coal_support == 0]) <- TRUE

is.na(elo2_yule_support[elo2_yule_support == 0]) <- TRUE
is.na(elo2_coal_support[elo2_coal_support == 0]) <- TRUE

is.na(its_yule_support[its_yule_support == 0]) <- TRUE
is.na(its_coal_support[its_coal_support == 0]) <- TRUE

is.na(lsu_yule_support[lsu_yule_support == 0]) <- TRUE
is.na(lsu_coal_support[lsu_coal_support == 0]) <- TRUE

is.na(mtlsu_yule_support[mtlsu_yule_support == 0]) <- TRUE
is.na(mtlsu_coal_support[mtlsu_coal_support == 0]) <- TRUE

is.na(mtssu_yule_support[mtssu_yule_support == 0]) <- TRUE
is.na(mtssu_coal_support[mtssu_coal_support == 0]) <- TRUE

is.na(rpb2_yule_support[rpb2_yule_support == 0]) <- TRUE
is.na(rpb2_coal_support[rpb2_coal_support == 0]) <- TRUE

is.na(tub_yule_support[tub_yule_support == 0]) <- TRUE
is.na(tub_coal_support[tub_coal_support == 0]) <- TRUE

# plot the tree and support values
plot(ef1a_yule_tr, cex=.6, no.margin=TRUE)
nodelabels(round(ef1a_yule_support, 2), cex=.7)

plot(ef1a_coal_tr, cex=.6, no.margin=TRUE)
nodelabels(round(ef1a_coal_support, 2), cex=.7)

plot(elo2_yule_tr, cex=.6, no.margin=TRUE)
nodelabels(round(elo2_yule_support, 2), cex=.7) 

plot(elo2_coal_tr, cex=.6, no.margin=TRUE)
nodelabels(round(elo2_coal_support, 2), cex=.7)

plot(its_yule_tr, cex=.6, no.margin=TRUE)
nodelabels(round(its_yule_support, 2), cex=.7) 

plot(its_coal_tr, cex=.6, no.margin=TRUE)
nodelabels(round(its_coal_support, 2), cex=.7)

plot(lsu_yule_tr, cex=.6, no.margin=TRUE)
nodelabels(round(lsu_yule_support, 2), cex=.7) 

plot(lsu_coal_tr, cex=.6, no.margin=TRUE)
nodelabels(round(lsu_coal_support, 2), cex=.7)

plot(mtlsu_yule_tr, cex=.6, no.margin=TRUE)
nodelabels(round(mtlsu_yule_support, 2), cex=.7) 

plot(mtlsu_coal_tr, cex=.6, no.margin=TRUE)
nodelabels(round(mtlsu_coal_support, 2), cex=.7)

plot(mtssu_yule_tr, cex=.6, no.margin=TRUE)
nodelabels(round(mtssu_yule_support, 2), cex=.7) 

plot(mtssu_coal_tr, cex=.6, no.margin=TRUE)
nodelabels(round(mtssu_coal_support, 2), cex=.7)

plot(rpb2_yule_tr, cex=.6, no.margin=TRUE)
nodelabels(round(rpb2_yule_support, 2), cex=.7) 

plot(rpb2_coal_tr, cex=.6, no.margin=TRUE)
nodelabels(round(rpb2_coal_support, 2), cex=.7)

plot(tub_yule_tr, cex=.6, no.margin=TRUE)
nodelabels(round(tub_yule_support, 2), cex=.7) 

plot(tub_coal_tr, cex=.6, no.margin=TRUE)
nodelabels(round(tub_coal_support, 2), cex=.7)

# results indicate that gymc does not delineate the species with confidence 
# as node values are lees than 1








