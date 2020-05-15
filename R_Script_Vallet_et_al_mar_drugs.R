########################### SCRIPT ################################
# Laminariales host does impact lipid temperaturetrajectories of the fungal endophyteParadendryphiella salina(G.K. Sutherl.)
# by : Marine Vallet, Tarik Meziane, Najet Thiney, Soizic Prado, Cédric Hubas*
# Submitted to Marine Drugs
# *Corresponding author : cedric.hubas@mnhn.fr
# Credits : Script = C.H., Data = M.V., N.T., S.P.
########################### SCRIPT ################################

#####################
# PACKAGES
#####################
library(reshape)
library(ggplot2)
library(rstatix)
library(ade4)
library(factoextra)
library(cowplot)

#####################
# DATA UPLOAD
#####################

# Optional ####################
#setwd("/Users/cedric.hubas/Desktop/Working_Directory/") # change working directory address if needed

# file extraction and upload ####################
#path <- "/Users/cedric.hubas/Desktop/Working_Directory/" # change working directory address if needed 
files <- list.files(pattern="*.txt") # if path changed use argument path = path
List <- lapply(files, function(x) read.table(x,skip=1)[,c(2,4)])

# Modify list names ####################
pos.point <- regexpr(".",files,fix=T) # use regular expression to find position of symbol "."
new.file.name <- substr(files, 1, pos.point-1)
names(List) <- new.file.name 
List

#####################
# DATA PROCESSING
#####################

# Remove internal standard (IS) ####################
List.WithoutC23 <- lapply(List, function(x)
	{vector=x[,1]!="23:0"
	x[vector,]
	})

# Caculate percentages ####################
percent <- lapply(List.WithoutC23, function(x)
	{fa=x[,2]
	total=sum(fa)
	data.frame(FA=x[,1],percent=fa/total*100)
	})

# Filter: remove peaks with percentage < threshold ####################
threshold <- 0
V <- lapply(percent,function(x) # identify peaks < threshold, generate TRUE/FALSE vector "V"
		{vector=x[,2]>= threshold
		})
super.list <- mapply(cbind, List.WithoutC23, V, SIMPLIFY=FALSE) # merge both lists

List.WithoutC23.WithoutThreshold <- lapply(super.list, function(x)
	{vector=x[,3]== TRUE
	x[vector,]
	})

# Caculate percentages again ####################
percent2 <- lapply(List.WithoutC23.WithoutThreshold, function(x)
	{fa=x[,2]
	total=sum(fa)
	data.frame(FA=x[,1],percent=fa/total*100)
	})

# Filter verification ####################
unlist(lapply(percent2, function(x) return(sum(x$percent/sum(x$percent)*100)))) # must give 100% 

#####################
# PERCENTAGES
#####################

dat <- do.call(rbind,percent2) # equivalent to function unlist
pos.point2 <- regexpr(".",rownames(dat),fix=T) # use regular expression to find position of symbol "."
new.file.name2 <- substr(rownames(dat), 1, pos.point2-1) # retrieve sample names
dat2 <- data.frame(dat,group=new.file.name2)
data.m <- melt(dat2,id=c(1,3),measure=c(2)) 
table <- cast(data.m, group ~ FA,sum)
table[is.na(table)] <- 0 ; table
colnames(table) <- gsub("w","n-",colnames(table))
table$sp <- substr(as.vector(table$group), 3, 4)
table$temp <- substr(as.vector(table$group), 5, 7)
table$salinity <- substr(as.vector(table$group), 8, 9)

# Optional ####################
#write.table(table,paste(getwd(),"/tables/FA.table.percent.txt",sep=""))

#####################
# CONCENTRATIONS
#####################

C.FA <- do.call(rbind,List)
pos.point3 <- regexpr(".",rownames(C.FA),fix=T) # use regular expression to find position of symbol "."
new.CFA.name <- substr(rownames(C.FA), 1, pos.point3-1) # retrieve sample names
C.FA2 <- data.frame(C.FA,group=new.CFA.name) 
C.FA2.m <- melt(C.FA2,id=c(1,3),measure=c(2)) 
table.C.FA=cast(C.FA2.m, group ~ V2,sum) 
table.C.FA[is.na(table.C.FA)] <- 0 ; table.C.FA

# Import IS table ####################
imported.C23 <- read.table(paste(getwd(),"/tables/fill.C23copy.txt",sep=""),h=T) # importer le fichier C23 avec les poids

column <- colnames(table.C.FA)!="23:0" # remove IS column
inv.column <- colnames(table.C.FA)=="23:0" # retreive IS column
Final.table.C.FA <- table.C.FA[,column] 
VectorC23 <- table.C.FA[,inv.column] # extract IS column
Concentration <- (imported.C23$C23mg*Final.table.C.FA[,-1])*1000/(VectorC23*imported.C23$Echmg) # calculate concentrations
Concentration$group <- table$group
colnames(Concentration) <- gsub("w","n-",colnames(Concentration))
Concentration$sp <- substr(as.vector(Concentration$group), 3, 4)
Concentration$temp <- substr(as.vector(Concentration$group), 5, 7)
Concentration$salinity <- substr(as.vector(Concentration$group), 8, 9)

# Optional ####################
#write.table(Concentration,paste(getwd(),"/tables/FA.table.conc.txt",sep=""))

#####################
# AESTHETICS
#####################

My_Theme = theme(
  axis.text.x = element_text(angle = 90,hjust = 1,size = 10),
  axis.text.y = element_text(size = 12),
  axis.title = element_text(size = 12),
  legend.title=element_text(size = 12),
  legend.text=element_text(size = 12),
  strip.background =element_rect(fill="white"))


#####################
# FIGURE 1
#####################

# modify table ####################
Workin.Conc <- Concentration[,1:22] 
dim <- dim(Workin.Conc)
TABLE.CONC <- data.frame(
	conc=as.vector(as.matrix(Workin.Conc)),
	fatty.acids=rep(names(Workin.Conc),each=dim[1]),
	sp=rep(Concentration$sp,dim[2]),
	temp=rep(Concentration$temp,dim[2]),
	salinity=rep(Concentration$salinity),dim[2]
)
TABLE.CONC$salinity <- c("Salinity 23.5 PSU","Salinity 50 PSU","Salinity 70 PSU")[as.factor(TABLE.CONC$salinity)]
TABLE.CONC$temp <- gsub("T","",TABLE.CONC$temp)

# get dominant FA ####################
TABLE.CONC2 <- TABLE.CONC[TABLE.CONC$fatty.acids=="16:0" |
	TABLE.CONC$fatty.acids=="16:1n-7" |
	TABLE.CONC$fatty.acids=="18:0" |
	TABLE.CONC$fatty.acids=="18:1n-7" |
	TABLE.CONC$fatty.acids=="18:1n-9" |
	TABLE.CONC$fatty.acids=="18:2n-6" |
	TABLE.CONC$fatty.acids=="18:3n-3" ,]

# Welch ANOVA ####################
my.conc <- Concentration[,c(5:6,14:18)]
s.my.conc <- split(my.conc,paste(Concentration$sp, Concentration$salinity))

group<-	as.factor(c(rep("10",3),rep("18",3),rep("25",3)))
res.over=NULL
for(j in 1:6) {
res.norm=NULL
res.bartlett =NULL
res.anova=NULL
for(i in 1:7){
	y <-s.my.conc[[j]][,i]
	res.norm[i] <- shapiro.test(aov(y~group)$res)$p
	res.bartlett[i] <- bartlett.test(aov(y~group)$res~group)$p.value
	bart.condition <- bartlett.test(aov(y~group)$res~group)$p.value
	if (bart.condition < 0.5){
     res.anova[i] <-  oneway.test(y~group, var.equal = FALSE)$p.value
}else{
     res.anova[i] <-  oneway.test(y~group, var.equal = TRUE)$p.value
}
	}
res.over[[j]] <- data.frame(
	FA=c("16:0","16:1n-7","18:0","18:1n-7","18:1n-9","18:2n-6","18:3n-3"),
	Shapiro.p.value= res.norm,
	Bartlett.p.value= res.bartlett,
	ANOVA.p.value=res.anova)
}

names(res.over) <- levels(factor(paste(Concentration$sp, Concentration$salinity))) ; res.over # ANOVA results

make_stars <- function(x) {
	stars<-NULL
for (i in 1:7) {
	if (x[i,4]>=0.05){
	stars[i] <- ""
}else{
		if (x[i,4]<0.05 & x[i,4]>=0.01) {
		stars[i] <- "*"
	}else{
			if (x[i,4]<0.01 & x[i,4]>=0.001) {
			stars[i] <- "**"
		}else{
				if (x[i,4]<0.001) {
				stars[i] <- "***"
			}else{
				stars[i] <- "***"
		}
	}
}}}
return(stars)
}

res.stars <- do.call(rbind,lapply(res.over,make_stars))
colnames(res.stars) <- c("16:0","16:1n-7","18:0","18:1n-7","18:1n-9","18:2n-6","18:3n-3") ; res.stars

# plot ####################

ggplot(TABLE.CONC2,aes(y=conc,x=fatty.acids,col=temp)) +
		geom_boxplot() + 
		theme_bw() + 
		facet_wrap(~paste(TABLE.CONC2$sp, " - ", TABLE.CONC2$salinity,sep="")) +
		scale_color_manual("Temperature (°C)",values=c("#81BDCA","#896587","#E98766")) +
		ylab(expression(paste("Concentration (µg.",g^-1,")")))+
		xlab("")+
		My_Theme

#####################
# FIGURE 2
#####################

myplot2 <- function (data,sal.name) {
multi=dudi.pca(data[,-c(1,24:26)],scannf=F,nf=2)
rownames(multi$li)=paste(1:dim(multi$li)[1],data$sp,data$temp,sep="")
colfac <- as.factor(data$sp)
Fac=as.factor(paste(data$sp,data$temp,sep=""))
X=aggregate(multi$li[,1],by=list(Fac),mean)
Y=aggregate(multi$li[,2],by=list(Fac),mean)

fviz_pca_ind(multi, geom = "point", pointshape=19,habillage=Fac, pointsize = 5) +
	labs(title = sal.name)+
	theme_bw()+
	scale_color_manual(values=c(rep("#FCCA28",3), rep("#FF7363",3))) +
	My_Theme +
	theme(legend.position = "none")+
	geom_segment(aes(x = X$x[1], 
		y = Y$x[1], 
		xend = X$x[2], 
		yend = Y$x[2]),
		arrow = arrow(length = unit(0.6, "cm")),col="#FCCA28") +
    geom_segment(aes(x = X$x[2], 
		y = Y$x[2], 
		xend = X$x[3], 
		yend = Y$x[3]),
		arrow = arrow(length = unit(0.6, "cm")),col="#FCCA28") +
    geom_segment(aes(x = X$x[4], 
		y = Y$x[4], 
		xend = X$x[5], 
		yend = Y$x[5]),
		arrow = arrow(length = unit(0.6, "cm")),col="#FF7363") +
    geom_segment(aes(x = X$x[5], 
		y = Y$x[5], 
		xend = X$x[6], 
		yend = Y$x[6]),
		arrow = arrow(length = unit(0.6, "cm")),col="#FF7363")
}

plotS1 <- myplot2(table[table$salinity=="S1",],sal.name="PCA - Ind: Salinity 23.5 PSU")
plotS2 <- myplot2(table[table$salinity=="S2",],sal.name="PCA - Ind: Salinity 50 PSU")
plotS3 <- myplot2(table[table$salinity=="S3",],"PCA - Ind: Salinity 70 PSU")

flitdata <- table[,]
flitdata$salinity <- gsub("S1", "Salinity 23.5 PSU",flitdata$salinity)
flitdata$salinity <- gsub("S2", "Salinity 50 PSU",flitdata$salinity)
flitdata$salinity <- gsub("S3", "Salinity 70 PSU",flitdata$salinity)
flitdata$index <- flitdata$"18:2n-6"/flitdata $"18:1n-9"
flitdata$temperature <- as.numeric(gsub("T","",flitdata$temp))

plotlinearindex3 <- ggplot(flitdata,aes(y=index,x=temperature,col=sp)) +
	facet_grid(~salinity)+
	geom_point() + 
	geom_smooth(method="lm")+
	theme_bw()+
	xlab("Temperature (°C)")+
	ylab(expression("FAI"[1][8][C]))+
	scale_color_manual(values=c("#FCCA28","#FF7363"))+
	theme(legend.position = "none")+
	My_Theme

upper <- plot_grid(plotS1,plotS2,plotS3,labels=c("a)", "b)", "c)"),nrow=1)
plot_grid(upper, plotlinearindex3,labels=c("", "d)"),nrow=2)

#####################
# ANCOVA
#####################

# data preparation ####################
anc.data <- split(flitdata,flitdata$salinity)

# Normality ####################
norm_check <- function(x) {
	model <- lm(index ~ temperature + sp, data = x)
	shapiro_test(model$res)
}
lapply(anc.data, norm_check) # Normality checked

# Homogeneity of variance ####################
var_check <- function(x) {
	model <- lm(index ~ temperature + sp, data = x)
	levene_test(model$res~as.factor(sp),data=x)
}
lapply(anc.data, var_check) # Homoscedasticity checked

# Linearity ####################
# checked visualy => ok

# Homogeneity of regression slopes ####################

anc.data[[1]] %>% anova_test(index~temperature*sp) # Salinity 23.5 PSU: interaction significant -> homogeneity not met
anc.data[[2]] %>% anova_test(index~temperature*sp) # Salinity 50 PSU: interaction significant -> homogeneity not met
anc.data[[3]] %>% anova_test(index~temperature*sp) #Salinity 70 PSU: interaction non significant -> homogeneity met

# Outliers ####################

boxplot(augment(model)$.std.resid)$out # no outliers

# ANCOVA ####################

anc.data[[3]] %>% anova_test(index~temperature+sp) #Salinity 70 PSU: models not significantly different

#####################
# Supp. Mat A1
#####################

names(Concentration)
aggregate(apply(Concentration[,1:22],1,sum),by=list(Concentration$sp, Concentration$temp, Concentration$salinity),mean)
aggregate(apply(Concentration[,1:22],1,sum),by=list(table$sp,table$temp,table$salinity),sd)
par(las=2)

sum.FA <- data.frame(FA.Concentration=apply(Concentration[,1:22],1,sum),
						sp= Concentration$sp,
						temperature=gsub("T","", Concentration$temp),
						salinity=c("Salinity 23.5 PSU","Salinity 50 PSU","Salinity 70 PSU")[as.factor(Concentration$salinity)])
						
ggplot(sum.FA,aes(y=FA.Concentration,x=temperature,col=sp)) +
		geom_boxplot() + 
		theme_bw() + 
		facet_grid(~salinity) +
		scale_color_manual(values=c("#FCCA28","#FF7363"))+
		ylab(expression(paste("Total fatty acid concentration (µg.",g^-1,")")))+
		xlab(expression(paste("Temperature (°C)")))+
		theme(strip.background = element_rect(fill="white"))+
		My_Theme
