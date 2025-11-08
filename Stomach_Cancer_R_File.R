library(data.table)
library(ggplot2)

df <- fread("C:/Users/lenna/OneDrive/Documents/Polygence/phenocode-KoGES_GCA.tsv/phenocode-KoGES_GCA.tsv",header=T,stringsAsFactors=F)
df$location = paste(df$chrom, df$pos, sep=":")

df3 <- fread("C:/Users/lenna/Downloads/gwas_gastric_hg19pos_l4QpCFE.tsv", header=T,stringsAsFactors = F)
df3$effect.al <-tstrsplit(df3$riskAllele,"-")[2]
df3$rsids <-tstrsplit(df3$riskAllele,"-")[1]
df3$idk <-tstrsplit(df3$rsids,"r")[2]

df3$dash <- df3$locations == "-"
df3$locations[df3$dash] <- df3$idk[df3$dash]
df3$chrom <-tstrsplit(df3$locations,":")[1]
df3$pos <-tstrsplit(df3$locations,":")[2]


df3.1<-df3[!(df3$pubmedId %in% c(35325318,35121771,37340842,30189721)),]

df3.1$hg19.location <- paste(df3.1$hg19.chr, df3.1$hg19.pos, sep=":")
df4<-merge(df,df3.1,by.x="location",by.y="hg19.location")


#remapping---------------------------------------------------

out.df <- data.frame(Original.SNP=character(), Original.rsid=character(), Korean.pval=numeric(),
                     New.SNP=character(),New.rsid=character(),New.pval=numeric(), 
                     Original.Pos=character(), New.Pos=character(),pval.apart=numeric(),chrom=character())

for(i in 1:nrow(df4)){
  tempchr <- df4$chrom.x[i]
  chr.df <- df[df$chrom==tempchr,]
  pos.df <- chr.df[chr.df$pos %between% c(df4$pos.x[i]-50000,df4$pos.x[i]+50000)]
  oldsnp <- paste(df4$chrom.x[i], df4$pos.x[i], sep=":")
  oldrsid <-df4$rsids.x[i]
  oldpval <- df4$pval[i]
  oldpos <- df4$pos.x[i]
  top.df <- pos.df[pos.df$pval==min(pos.df$pval),]
  newsnp <- paste(top.df$chrom, top.df$pos, sep=":")
  newrsid <-top.df$rsids
  newpval <- top.df$pval
  newpos <- top.df$pos
  pval.apart <- oldpval-newpval
  chrom.out.df<-top.df$chrom

  tmp <- data.frame(Original.SNP=oldsnp, Original.rsid=oldrsid, Original.pval=oldpval,
                    New.SNP=newsnp,New.rsid=newrsid,New.pval=newpval, 
                    Original.Pos=oldpos, New.Pos=newpos,pval.apart=pval.apart,chrom=chrom.out.df)
  out.df <- rbind(out.df,tmp)
}

out.df$distance.apart <- (out.df$Original.Pos-out.df$New.Pos)


#----------------------------------------------

ea.df <- data.frame(GWASCAT.SNP = df4$location)
ea.df$Korean.SNP <- df4$location
ea.df$published.allele <- df4$effect.al
ea.df$alt.allele <- df4$alt
ea.df$ref.allele <- df4$ref
ea.df$beta.value <- df4$beta.x
ea.df$catelog.OR <- df4$orValue
 
ea.df <- subset(ea.df, catelog.OR!="-")
ea.df <- subset(ea.df, published.allele!="?")

ea.df$alt.match <- ea.df$alt.allele==ea.df$published.allele
ea.df$ref.match <-ea.df$ref.allele==ea.df$published.allele

both.al.mismatch <- ea.df$alt.match==F & ea.df$ref.match==F
alt.A <- ea.df$alt.allele=="A"
alt.C <- ea.df$alt.allele=="C"
alt.T <- ea.df$alt.allele=="T"
alt.G <- ea.df$alt.allele=="G"
ea.df$alt.allele[both.al.mismatch & alt.A] <- "T"
ea.df$alt.allele[both.al.mismatch & alt.C] <- "G"
ea.df$alt.allele[both.al.mismatch & alt.T] <- "A"
ea.df$alt.allele[both.al.mismatch & alt.G] <- "C"

ref.A <- ea.df$ref.allele=="A"
ref.C <- ea.df$ref.allele=="C"
ref.T <- ea.df$ref.allele=="T"
ref.G <- ea.df$ref.allele=="G"
ea.df$ref.allele[both.al.mismatch & ref.A] <- "T"
ea.df$ref.allele[both.al.mismatch & ref.C] <- "G"
ea.df$ref.allele[both.al.mismatch & ref.T] <- "A"
ea.df$ref.allele[both.al.mismatch & ref.G] <- "C"

ea.df$alt.match <- ea.df$alt.allele==ea.df$published.allele
ea.df$ref.match <-ea.df$ref.allele==ea.df$published.allele

beta.ltmt <- ea.df$beta.value<0
OR.ltmt<- ea.df$published.allele<1
ea.df$direction.match <- beta.ltmt==OR.ltmt

newbeta <- ifelse(ea.df$alt.match==F, ea.df$beta.value*-1,ea.df$beta.value)
newbeta.ltmt<- newbeta<0
ea.df$direction.match <-newbeta.ltmt==OR.ltmt
  
#------------------------------------------------
 #interpopulation allele frequencies

af.df <- fread("C:/Users/lenna/OneDrive/Documents/Polygence/gwas_cataloghg19_fq_all.txt",header=T,stringsAsFactors=F,col.names = c("CHROM","POS","AF_EUR_unrel","AF_EAS_unrel","AF_AFR_unrel","AF_SAS_unrel","AF_AMR_unrel"))
af.df$location = paste(af.df$CHROM, af.df$POS, sep=":")
bf.df<-merge(af.df,ea.df,by.x="location",by.y="Korean.SNP")

#polarization----------------------------------------
bf.df$flip <- F
bf.df$flip[bf.df$alt.match==T & bf.df$catelog.OR<1]<-T
bf.df$flip[bf.df$ref.match==T & bf.df$catelog.OR>1]<-T

bf.df$AF_EUR_unrel_new<-bf.df$AF_EUR_unrel
EUR_flipped <- 1-bf.df$AF_EUR_unrel
bf.df$AF_EUR_unrel_new[bf.df$flip==T]<- EUR_flipped[bf.df$flip==T]

bf.df$AF_EAS_unrel_new<-bf.df$AF_EAS_unrel
EAS_flipped <- 1-bf.df$AF_EAS_unrel
bf.df$AF_EAS_unrel_new[bf.df$flip==T]<- EAS_flipped[bf.df$flip==T]

bf.df$AF_AFR_unrel_new<-bf.df$AF_AFR_unrel
AFR_flipped <- 1-bf.df$AF_AFR_unrel
bf.df$AF_AFR_unrel_new[bf.df$flip==T]<- AFR_flipped[bf.df$flip==T]

bf.df$AF_SAS_unrel_new<-bf.df$AF_SAS_unrel
SAS_flipped <- 1-bf.df$AF_SAS_unrel
bf.df$AF_SAS_unrel_new[bf.df$flip==T]<- SAS_flipped[bf.df$flip==T]

bf.df$AF_AMR_unrel_new<-bf.df$AF_AMR_unrel
AMR_flipped <- 1-bf.df$AF_AMR_unrel
bf.df$AF_AMR_unrel_new[bf.df$flip==T]<- AMR_flipped[bf.df$flip==T]

#allele-frequency-comparison----------------------------

bf.df$EURDIF <- bf.df$AF_EAS_unrel_new-bf.df$AF_EUR_unrel_new
cf.df <- data.frame(EURDIFmean = mean(bf.df$EURDIF))

bf.df$AFRDIF <- bf.df$AF_EAS_unrel_new-bf.df$AF_AFR_unrel_new
cf.df$AFRDIFmean <- mean(bf.df$AFRDIF)

bf.df$SASDIF <- bf.df$AF_EAS_unrel_new-bf.df$AF_SAS_unrel_new
cf.df$SASDIFmean <- mean(bf.df$SASDIF)

bf.df$AMRDIF <- bf.df$AF_EAS_unrel_new-bf.df$AF_AMR_unrel_new
cf.df$AMRDIFmean <- mean(bf.df$AMRDIF)

#remapping-visualization---------------------------------------------------

out.df$log.newpval <- -log10(out.df$New.pval)
out.df$log.originalpval <- -log10(out.df$Original.pval)

out.df2 <- data.frame(out.df[1,])

ggplot(out.df) + geom_point(aes(x = distance.apart, y = log.newpval, color="New SNPS"))+ 
  geom_point(aes( x = 0, y = log.originalpval, color="Original SNPS"), alpha=0.5 ) + 
  xlab("Distance from original SNP (base pairs)")+ ylab("-log10(P-Value)") +
  geom_segment(aes(x = 0, y = log.originalpval, xend = distance.apart, yend = log.newpval),
              colour = "black", alpha=0.5, linetype='solid')+
  theme_classic(base_size = 15) + theme(legend.position = "none") 

ggplot(out.df2) + geom_point(aes(x = distance.apart, y = log.newpval, color="New SNPS", size=3))+ 
  geom_point(aes( x = 0, y = log.originalpval, color="Original SNPS", size=3), alpha=0.5 ) + 
  scale_x_continuous(name="Distance from original SNP (base pairs)", limits=c(-50000, 50000))+ ylab("-log10(P-Value)") +
  geom_segment(aes(x = 0, y = log.originalpval, xend = distance.apart, yend = log.newpval),
               colour = "black", alpha=0.5, linetype='solid')+
  theme_classic(base_size = 15) + theme(legend.position = "none")

#ggsave("C:/Users/lenna/OneDrive/Documents/Polygence/Science Fair Figures/remapping2.png", dpi = 300)

#_direction of effect analysis vizualization__________________________________________________________________-
n.true <- table(ea.df$direction.match)[2]
n.false <- table(ea.df$direction.match)[1]
n.total <- nrow(ea.df)

ea.df.2 <- data.frame(Match=c("TRUE", "FALSE"), Percent=c(n.true/n.total,n.false/n.total))

ggplot(ea.df.2, aes(x=Match, y=Percent, fill=Match)) +
  geom_bar(stat="identity", width=0.7)+ ylab("Proportion")+ xlab("GWAS association direction of effect match") +
  theme_classic(base_size = 20)+ theme(legend.position = "none") + scale_y_continuous(expand=c(0,0))

ggsave("C:/Users/lenna/OneDrive/Documents/Polygence/Science Fair Figures/directionofeffect.png", dpi = 300)

#binomial test----------------------------------------------------------

binom.test(n.true, n.total)

#paired-T-Test-----------------------------------------------------------

t.test(x=bf.df$AF_EAS_unrel_new, y=bf.df$AF_EUR_unrel_new, paired=T)
t.test(x=bf.df$AF_EAS_unrel_new, y=bf.df$AF_AFR_unrel_new, paired=T)
t.test(x=bf.df$AF_EAS_unrel_new, y=bf.df$AF_SAS_unrel_new, paired=T)
t.test(x=bf.df$AF_EAS_unrel_new, y=bf.df$AF_AMR_unrel_new, paired=T)

#allele-frequency-vizualization----------------------------

cf.df$EURmean <- mean(bf.df$AF_EUR_unrel_new)
cf.df$EASmean <- mean(bf.df$AF_EAS_unrel_new)
cf.df$AFRmean <- mean(bf.df$AF_AFR_unrel_new)
cf.df$SASmean <- mean(bf.df$AF_SAS_unrel_new)
cf.df$AMRmean <- mean(bf.df$AF_AMR_unrel_new)

df.df<- data.frame(Region=c("EUR", "EAS", "AFR", "SAS", "AMR"), AlleleFrequency=c(cf.df$EURmean,cf.df$EASmean,cf.df$AFRmean,cf.df$SASmean,cf.df$AMRmean))

ggplot(df.df, aes(x=Region, y=AlleleFrequency, fill=Region)) +
  geom_bar(stat="identity", width=0.7)+ ylab("Mean Allele Frequency")+ xlab("Region") +
  theme_classic(base_size = 20)+ theme(legend.position = "none") + scale_y_continuous(expand=c(0,0))

#ggsave("C:/Users/lenna/OneDrive/Documents/Polygence/Science Fair Figures/allelefrequencies.png", dpi = 300)

#allele-frequency-vizualization2--------------------------

ggplot(bf.df) + geom_point(aes(x = "East Asia", y = AF_EAS_unrel_new, color="brown1"))+ 
  geom_point(aes( x = "South Asia", y = AF_SAS_unrel_new, color="darkgreen")) + 
  xlab("Region")+ ylab("Risk allele frequencies") +
  geom_segment(aes(x = 1, y = AF_EAS_unrel_new, xend = 2, yend = AF_SAS_unrel_new),
               colour = "black", alpha=0.5, linetype='solid')+
  theme_classic(base_size = 15) + theme(legend.position = "none") 

ggsave("C:/Users/lenna/OneDrive/Documents/Polygence/Science Fair Figures/allelefrequenciesSAS.png", dpi = 300)

#allele-frequency-#2---------------------------------------------------------------------------

write.table(out.df[,c("chrom","New.Pos")],file="C:/Users/lenna/OneDrive/Documents/Polygence/Science Fair Figures/New.SNP.Locations",sep="\t",quote=F,row.names=F,col.names=F)

#------------------------------------------------
#interpopulation allele frequencies #2

af.2df <- fread("C:/Users/lenna/OneDrive/Documents/Polygence/gwas_remapped_fq_all.txt",header=T,stringsAsFactors=F,col.names = c("CHROM","POS","AF_EUR_unrel","AF_EAS_unrel","AF_AFR_unrel","AF_SAS_unrel","AF_AMR_unrel"))
af.2df$location = paste(af.2df$CHROM, af.2df$POS, sep=":")

df5<-merge(df,af.2df,by.x="location",by.y="location")
#polarization #2----------------------------------------
df5$flip <- F
df5$flip[df5$beta<0]<-T

df5$AF_EUR_unrel_new<-df5$AF_EUR_unrel
EUR_flipped2 <- 1-df5$AF_EUR_unrel
df5$AF_EUR_unrel_new[df5$flip==T]<- EUR_flipped2[df5$flip==T]

df5$AF_EAS_unrel_new<-df5$AF_EAS_unrel
EAS_flipped2 <- 1-df5$AF_EAS_unrel
df5$AF_EAS_unrel_new[df5$flip==T]<- EAS_flipped2[df5$flip==T]

df5$AF_AFR_unrel_new<-df5$AF_AFR_unrel
AFR_flipped2 <- 1-df5$AF_AFR_unrel
df5$AF_AFR_unrel_new[df5$flip==T]<- AFR_flipped2[df5$flip==T]

df5$AF_SAS_unrel_new<-df5$AF_SAS_unrel
SAS_flipped2 <- 1-df5$AF_SAS_unrel
df5$AF_SAS_unrel_new[df5$flip==T]<- SAS_flipped2[df5$flip==T]

df5$AF_AMR_unrel_new<-df5$AF_AMR_unrel
AMR_flipped2 <- 1-df5$AF_AMR_unrel
df5$AF_AMR_unrel_new[df5$flip==T]<- AMR_flipped2[df5$flip==T]

#allele-frequency-comparison #2----------------------------

df5$EURDIF <- df5$AF_EAS_unrel_new-df5$AF_EUR_unrel_new
cf.df2 <- data.frame(EURDIFmean = mean(df5$EURDIF))

df5$AFRDIF <- df5$AF_EAS_unrel_new-df5$AF_AFR_unrel_new
cf.df2$AFRDIFmean <- mean(df5$AFRDIF)

df5$SASDIF <- df5$AF_EAS_unrel_new-df5$AF_SAS_unrel_new
cf.df2$SASDIFmean <- mean(df5$SASDIF)

df5$AMRDIF <- df5$AF_EAS_unrel_new-df5$AF_AMR_unrel_new
cf.df2$AMRDIFmean <- mean(df5$AMRDIF)

cf.df2$EURmean <- mean(df5$AF_EUR_unrel_new)
cf.df2$EASmean <- mean(df5$AF_EAS_unrel_new)
cf.df2$AFRmean <- mean(df5$AF_AFR_unrel_new)
cf.df2$SASmean <- mean(df5$AF_SAS_unrel_new)
cf.df2$AMRmean <- mean(df5$AF_AMR_unrel_new)
#paired-T-Test2-----------------------------------------------------------

t.test(x=df5$AF_EAS_unrel_new, y=df5$AF_EUR_unrel_new, paired=T)
t.test(x=df5$AF_EAS_unrel_new, y=df5$AF_AFR_unrel_new, paired=T)
t.test(x=df5$AF_EAS_unrel_new, y=df5$AF_SAS_unrel_new, paired=T)
t.test(x=df5$AF_EAS_unrel_new, y=df5$AF_AMR_unrel_new, paired=T)

#allele-frequency-vizualization2 #2--------------------------

ggplot(df5) + geom_point(aes(x = "East Asia", y = AF_EAS_unrel_new, color="brown1"))+ 
  geom_point(aes( x = "Europe", y = AF_EUR_unrel_new, color="darkgreen")) + 
  xlab("Region")+ ylab("Risk allele frequencies") +
  geom_segment(aes(x = 1, y = AF_EAS_unrel_new, xend = 2, yend = AF_EUR_unrel_new),
               colour = "black", alpha=0.5, linetype='solid')+
  theme_classic(base_size = 15) + theme(legend.position = "none") 

#ggsave("C:/Users/lenna/OneDrive/Documents/Polygence/Science Fair Figures/allelefrequenciesEUR2.png", dpi = 300)

