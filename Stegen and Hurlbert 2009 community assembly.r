## code for simulating 'metacommunities' and assembling local communities
## and calculating taxonomic, functional and phylogenetic beta diversity
## in order to evaluate the performance of different beta diversity metrics
## this version evolves traits and assigns spatial locations based on the assumed trait-space correlation
## which is itself a function of dispersal limitation and environmental filtering 
## this version has two traits, two environmental variables and two spatial dimensions
## this version evaluates the variance partitioning across 11 niche and dispersal breadths, but only one env.-space correlation

date()

library(ape)
library(geiger)
library(apTreeshape)
library(picante)
library(lattice)
individuals = 10000 # number of individuals in a local community
communities = 20 # number of local communities
meta.rich = 500 # metacommunity richness

## start parameters to manipulate

es.sd = 0.9 ## standard deviation for the environment-space correlation
rsq.desired=0.5 ## desired space-environment r.sq

#for (niche.breadth in c(10^seq(0,1,0.5))) { # variance of gaussian niche curve

#for (disp.breadth in c(10^seq(-4,1,0.5))) { # variance of gaussian dispersal kernal

###########################
# Variables to manipulate #
###########################
param_values = c(10^-9, 10^-4, 10^1)

for (niche.breadth in param_values) {
    for (disp.breadth in param_values) {
        
        meta.abun = rlnorm(meta.rich, meanlog = 1, sdlog = 1) ## "global" species abundance distribution

        ## end parameters to manipulate

        phylo.comp = matrix(c(0),nrow=0,ncol=23)
        trait.comp = matrix(c(0),nrow=0,ncol=23)
        taxa.comp = matrix(c(0),nrow=0,ncol=23)

        ## start simulated metacommunity, which is the overall phylogeny

        ##windows(21,14)
        ##par(mfrow=c(2,3))

        my.d = 0 #death rate - when set to zero = pure birth model
        my.b = 0.1 #birth rate

        tree1 = birthdeath.tree(b=my.b, d=my.d, taxa.stop=(meta.rich+1), return.all.extinct=FALSE) #generate random ultrametric trees with given richness
        tree2 = as.treeshape(tree1) #intermediate step to estimate tree imbalance
        meta = rescaleTree(tree1, 1) #standardise tree root-to-tip length

        Ic = colless(tree2,norm = "yule") # tree imbalance
        PD = sum(meta$edge.length) # alpha phylodiversity
        gamma = gammaStat(meta) # related to stemminess
        meta = cophenetic(meta)/2;
        rownames(meta)[rownames(meta)=="s501"]="del"
        colnames(meta)[colnames(meta)=="s501"]="del"
        meta = subset(meta, select = -del) 
        meta = subset(t(meta), select = -del) 
        meta = t(meta); dim(meta); range(as.numeric(rownames(meta)));

        ##########################
        # PHYLOGENY 
        ##########################
        meta = as.phylo(hclust(as.dist(meta)))
        phylo.dist = round(cophenetic(meta)/2,4); phylo.dist[1:5,1:5];
        ##plot(meta,typ="fan",show.tip.label = F,main="Metacommunity Phylogeny")

        ## end simulation of metacommunity, which is now meta

        ## start assinging environmental conditions and spatial locations to communities

        ##windows(14,7); par(mfrow=c(1,2))

        # first spatial and envirnmental axes
        env.opt.one = matrix(c(runif(communities, min=-10, max=10)),ncol=1,nrow=communities) # randomly selecting the environmental optima for each trait for each local community
        for (i in 1:ncol(env.opt.one)) {env.opt.one[,i] = (env.opt.one[,i]-mean(env.opt.one[,i]))/sd(env.opt.one[,i])} # standardizing environmental optima into z-scores

        for (loop in 1:1000) {

            comm.space.one = rnorm(communities, mean = env.opt.one, sd = es.sd);
            comm.space.one = (comm.space.one-mean(comm.space.one))/sd(comm.space.one); comm.space.one = matrix(c(comm.space.one),ncol=1) # standardizing community spatial locations into z-scores
            space.env.r.sq = round(summary(lm(comm.space.one~env.opt.one[,1]))$r.sq,2)
            if ( abs(space.env.r.sq - rsq.desired) < 0.05 ) {break} else{} 

        }

        plot(comm.space.one~env.opt.one[,1],main="1st axes Space-Environment Correlation",ylab="Local Spatial Position",xlab="Local Environmental Value"); 
        legend(min(env.opt.one[,1]),max(comm.space.one),legend = paste("R.sq = ",space.env.r.sq,sep=""),box.lty = 0)

        space.env.cov.one = cov(comm.space.one,env.opt.one[,1])

        # second spatial and envirnmental axes
        env.opt.two = matrix(c(runif(communities, min=-10, max=10)),ncol=1,nrow=communities) # randomly selecting the environmental optima for each trait for each local community
        for (i in 1:ncol(env.opt.two)) {env.opt.two[,i] = (env.opt.two[,i]-mean(env.opt.two[,i]))/sd(env.opt.two[,i])} # standardizing environmental optima into z-scores

        for (loop in 1:1000) {

            comm.space.two = rnorm(communities, mean = env.opt.two, sd = es.sd);
            comm.space.two = (comm.space.two-mean(comm.space.two))/sd(comm.space.two); comm.space.two = matrix(c(comm.space.two),ncol=1) # standardizing community spatial locations into z-scores
            space.env.r.sq = round(summary(lm(comm.space.two~env.opt.two[,1]))$r.sq,2)
            if ( abs(space.env.r.sq - rsq.desired) < 0.05 ) {break} else{} 

        }

        plot(comm.space.two~env.opt.two[,1],main="2nd axes Space-Environment Correlation",ylab="Local Spatial Position",xlab="Local Environmental Value"); 
        legend(min(env.opt.two[,1]),max(comm.space.two),legend = paste("R.sq = ",space.env.r.sq,sep=""),box.lty = 0)

        space.env.cov.two = cov(comm.space.two,env.opt.two[,1])

        ###################################################
        # SITE X SPATIAL COORDINATES MATRIX
        ###################################################
        site.xy = data.frame(x = comm.space.one, y = comm.space.two)




        ## end assinging environmental conditions and spatial locations to communities

        ## start calculation of the trait-space covariance

        trait.space.cov.1.1 = (exp(-disp.breadth) + space.env.cov.one*exp(-niche.breadth))/3  ## the covariance between traits and space along the first axes
        trait.space.cov.2.2 = (exp(-disp.breadth) + space.env.cov.two*exp(-niche.breadth))/3  ## the covariance between traits and space along the second axes
        trait.space.cov.1.2 = ( exp(-disp.breadth) )/3  ## the covariance between first trait and second space axes
        trait.space.cov.2.1 = ( exp(-disp.breadth) )/3  ## the covariance between second trait and first space axes
        trait.trait.cov = 0 ## covariance between traits
        space.space.cov = 0 ## covariance between spatial dimensions
        all.var = 1 ## variance of space and traits

        ## end calculation of the trait-space covariance


        ## start trait/space evolution, brownian motion

        traits = 4 # number of space + trait dimensions
        var.cov = matrix(c(
        all.var,trait.trait.cov,trait.space.cov.1.1,trait.space.cov.1.2,
        trait.trait.cov,all.var,trait.space.cov.2.1,trait.space.cov.2.2,
        trait.space.cov.1.1,trait.space.cov.1.2,all.var,space.space.cov,
        trait.space.cov.2.1,trait.space.cov.2.2,space.space.cov,all.var)
        ,ncol=traits,nrow=traits); var.cov

        ##for (loop in 1:1000) {
        trait.brown = sim.char(meta,var.cov) ## matrix of traits

        t.s.matrix = matrix(c(trait.brown),nrow=meta.rich); rownames(t.s.matrix) = rownames(trait.brown)
        for (i in 1:ncol(t.s.matrix)) {t.s.matrix[,i] = (t.s.matrix[,i]-mean(t.s.matrix[,i]))/sd(t.s.matrix[,i])} # standardizing traits into z-scores

        cov(t.s.matrix)-var.cov

        ##real.trait.space.cov = cov(t.s.matrix[,1],t.s.matrix[,2])
        ##if ( abs(real.trait.space.cov - trait.space.cov) < 0.05 ) {break} else{}
        ##}

        trait.matrix = rbind(cbind(env.opt.one,env.opt.two),as.matrix(t.s.matrix[,c(1,2)],ncol=2))
        for (i in 1:communities) {rownames(trait.matrix)[i] = paste("c",i,sep="") }
        trait.matrix = as.matrix(dist(trait.matrix,upper=T,diag=T)); trait.matrix = trait.matrix/max(trait.matrix); trait.matrix[1:5,1:5]; dim(trait.matrix) # distance matrix made with the environmental optima

        space.matrix = rbind(cbind(comm.space.one,comm.space.two),as.matrix(t.s.matrix[,c(3,4)],ncol=2))
        for (i in 1:communities) {rownames(space.matrix)[i] = paste("c",i,sep="") }
        space.matrix = as.matrix(dist(space.matrix,upper=T,diag=T)); space.matrix = space.matrix/max(space.matrix); space.matrix[1:5,1:5]; dim(space.matrix) # distance matrix made with the environmental optima

        #K.trait = Kcalc(trait.brown,meta) ## trait phylogenetic signal
        #K.trait

        ##trait.dendro = as.phylo(hclust(as.dist(trait.matrix)))
        ##space.dendro = as.phylo(hclust(as.dist(space.matrix)))
        ##trait.dendro.fsor = as.phylo(hclust(as.dist(trait.matrix[(communities+1):nrow(trait.matrix),(communities+1):ncol(trait.matrix)])))
        ##plot(trait.dendro,typ="fan",show.tip.label = F,cex=0.5,main="Metacommunity Trait Dendrogram")
        ##plot(space.dendro,typ="fan",show.tip.label = F,cex=0.5,main="Metacommunity Space Dendrogram")

        ## end trait/space evolution, brownian motion

        ## start community assembly, abundance and environment

        species = rownames(trait.matrix)[(communities+1):nrow(trait.matrix)]

        pres.ab =  matrix(c(0),ncol=meta.rich,nrow=communities)
        rich = numeric()
        for(i in 1:nrow(pres.ab)) {

            ##local.species = as.numeric(sample(species,local.rich,replace=F,prob=c(exp(-((trait.matrix[(communities+1):nrow(trait.matrix),i])^2)/(2*niche.breadth))/(sqrt(2*niche.breadth*pi))))) # niche breadth multiplied by 10 to get more rare species
            ##rel.abun = sample(local.species,individuals,replace=T,prob=c(exp(-((trait.matrix[c(local.species+communities),i])^2)/(2*niche.breadth))/(sqrt(2*niche.breadth*pi))))
            
          # THE PROBLEM IS HERE: the expression passed to prob is just a matrix full of 0's,
          # resulting in "too few positive probabilities"
          rel.abun = as.numeric(sample(species,individuals,replace=T,
                prob=c(
	                meta.abun ## probaility from the global abundance distribution
	                *exp(-((trait.matrix[c(communities+1):nrow(trait.matrix),i])^2)/(2*niche.breadth))/(sqrt(2*niche.breadth*pi)) ## probability from environmental filtering
	                *exp(-((space.matrix[c(communities+1):nrow(space.matrix),i])^2)/(2*disp.breadth))/(sqrt(2*disp.breadth*pi)) ## probability from spatial dispersal limitation
            )))
            rel.abun = as.data.frame(rel.abun)
            rel.abun = cbind(rel.abun,1)
            rel.abun = tapply(rel.abun[,2],rel.abun[,1],FUN=sum)
            rel.abun = cbind(as.numeric(names(rel.abun)),rel.abun)
            rel.abun[,2] = rel.abun[,2]/sum(rel.abun[,2])
            rich[i] = nrow(rel.abun)
            pres.ab[i,c(rel.abun[,1])] = rel.abun[,2]

        }
        pres.ab[1:5,1:5]
        ##hist(log(rel.abun[,2]),breaks=20,xlab="Relative Abundance",ylab="# of Species in Local Community",main="Example Local SAD")

        coms = cbind(trait.matrix[1:communities,1],space.matrix[1:communities,1],pres.ab)
        colnames(coms) = c("Env","Spatial",c(1:meta.rich))
        coms[1:5,1:5]

        ###################################################
        # SITE X SPECIES MATRIX IN EITHER REL.ABUN OR COMS
        ###################################################

        ## end community assembly, abundance and environment
    }
}
