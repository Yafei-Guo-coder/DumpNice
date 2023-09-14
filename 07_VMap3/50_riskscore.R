#riskscore
library("PredictABEL")
setwd("/Users/guoyafei/Desktop/riskscore")
#selection
data <- read.table("all.input.txt", header=F, stringsAsFactors = F)
colnames(data) <- c("ID","Latitude","Longitude","Elevation","type","ID","prec1","prec2","prec3","prec4","prec5","prec6","prec7","prec8","soil1","soil10","soil11","soil12","soil13","soil2","soil3","soil4","soil5","soil6","soil7","soil8","soil9","solar1","solar2","solar3","temp1","temp10","temp11","temp2","temp3","temp4","temp5","temp6","temp7","temp8","temp9",paste("rs",1:85,sep="") )
sub <- data[,c(5,2,3,4,7:14,16:19,27:85)]
decline <- c(1,2,3,4,6,8,9,10,11,13,14,16,17,19,20,22,23)
nodecline <- c(5,7,12,15,18,21,24,25)
sub$outcome <- NA 
sub[which(sub$type %in% decline), 76] <- 1
sub[which(sub$type %in% nodecline), 76] <- 0
cOutcome <- 76
#cNonGenPred1 <- c(5:12)
#cNonGenPred1 <- c(13:17)
#cNonGenPred1 <- c(18:20)
#cNonGenPred1 <- c(21:31)
#cNonGenPred1 <- c(5:31)
cNonGenPred1 <- c(0)
cNonGenPredCat1 <- c(0)
#cGenPred1 <- c(0)
cGenPred1 <- c(32:75)
cGenPredsCat1 <- c(0)


#environment related
data <- read.table("all.input2.txt", header=F, stringsAsFactors = F)
colnames(data) <- c("ID","Latitude","Longitude","Elevation","type","ID","prec1","prec2","prec3","prec4","prec5","prec6","prec7","prec8","soil1","soil10","soil11","soil12","soil13","soil2","soil3","soil4","soil5","soil6","soil7","soil8","soil9","solar1","solar2","solar3","temp1","temp10","temp11","temp2","temp3","temp4","temp5","temp6","temp7","temp8","temp9",paste("rs",1:29467,sep="") )
sub <- data[,c(5,2,3,4,7:14,16:19,27:29508)]
decline <- c(1,2,3,4,6,8,9,10,11,13,14,16,17,19,20,22,23)
nodecline <- c(5,7,12,15,18,21,24,25)
sub$outcome <- NA 
sub[which(sub$type %in% decline), 29499] <- 1
sub[which(sub$type %in% nodecline), 29499] <- 0

cOutcome <- 29499
cNonGenPred1 <- c(0)
cNonGenPredCat1 <- c(0)
cGenPred1 <- c(sample(50:29467,85))
cGenPredsCat1 <- c(0)

#random
data <- read.table("random.all.input.txt", header=F, stringsAsFactors = F)
colnames(data) <- c("ID","Latitude","Longitude","Elevation","type","ID","prec1","prec2","prec3","prec4","prec5","prec6","prec7","prec8","soil1","soil10","soil11","soil12","soil13","soil2","soil3","soil4","soil5","soil6","soil7","soil8","soil9","solar1","solar2","solar3","temp1","temp10","temp11","temp2","temp3","temp4","temp5","temp6","temp7","temp8","temp9",paste("rs",1:283,sep="") )
sub <- data[,c(5,2,3,4,7:14,16:19,27:283)]
decline <- c(1,2,3,4,6,8,9,10,11,13,14,16,17,19,20,22,23)
nodecline <- c(5,7,12,15,18,21,24,25)
sub$outcome <- NA 
sub[which(sub$type %in% decline), 274] <- 1
sub[which(sub$type %in% nodecline), 274] <- 0
cOutcome <- 274
cNonGenPred1 <- c(0)
cNonGenPredCat1 <- c(0)
cGenPred1 <- c(sample(50:273,85))
cGenPredsCat1 <- c(0)


########################
riskmodel <- fitLogRegModel(data=sub, cOutcome=cOutcome, cNonGenPreds=cNonGenPred1, cNonGenPredsCat=cNonGenPredCat1, cGenPreds=cGenPred1, cGenPredsCat=cGenPredsCat1)
summary(riskmodel)
ORmultivariate(riskModel=riskmodel)
predRisk <- predRisk(riskmodel)
all <-as.data.frame(cbind(predRisk, sub$type,sub$outcome))
all$V2 <- factor(all$V2, levels = c(1:25))
all$V3 <- as.factor(all$V3)

pdf("test1-random3.pdf", width=5, height=8)

p <- ggplot(all, aes(x = predRisk,group = V2, fill=V3))+
  geom_boxplot() +
  theme_classic() + 
  scale_fill_manual(values = c("#215085","#D53B3D")) +
  xlab("Polygeneic risk score")+
  ylab("Population")+
  theme(
    legend.position="none",
    #panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x=element_text(size=15),
    axis.text.y=element_blank(),
    axis.title.y=element_text(size = 18),
    axis.title.x=element_text(size = 18),
) 
print(p)
dev.off()

rangeaxis <- c(0,1) 
groups <- 10
pdf("test2-random3.pdf")
plotCalibration(data=sub, cOutcome=cOutcome, predRisk=predRisk, groups=groups, rangeaxis=rangeaxis)
print(p)
dev.off()

labels <- c("no-decline","decline")
pdf("test3-random3.pdf")
plotDiscriminationBox(data=sub, cOutcome=cOutcome, predrisk=predRisk, labels=labels)
print(p)
dev.off()

#3.plotCalibration

predRisk <- predRisk(riskmodel)
rangeaxis <- c(0,1) 
groups <- 20
pdf("1.pdf")
plotCalibration(data=sub2, cOutcome=cOutcome, predRisk=predRisk, groups=groups, rangeaxis=rangeaxis)
print(p)
dev.off()

#4.plotDiscriminationBox

labels <- c("decline", "no-decline")
pdf("test2.pdf")
plotDiscriminationBox(data=sub2, cOutcome=cOutcome, predrisk=predRisk, labels=labels)
print(p)
dev.off()

#5.plotPredictivenessCurve
data(ExampleData)
riskmodel1 <- ExampleModels()$riskModel1 
riskmodel2 <- ExampleModels()$riskModel2
predRisk1 <- predRisk(riskmodel) 
predRisk2 <- predRisk(riskmodel2)
rangeyaxis <- c(0,1) 
labels <- c("without genetic factors", "with genetic factors")
pdf("test3.pdf")
plotPredictivenessCurve(predrisk=cbind(predRisk1,predRisk2), rangeyaxis=rangeyaxis, labels=labels)
print(p)
dev.off()

#6.plotPriorPosteriorRisk
data(ExampleData) 
cOutcome <- 2
riskmodel1 <- ExampleModels()$riskModel1 
riskmodel2 <- ExampleModels()$riskModel2
predRisk1 <- predRisk(riskmodel1) 
predRisk2 <- predRisk(riskmodel2)
xlabel <- "Prior risk" 
ylabel <- "Posterior risk" 
titleplot <- "Prior versus posterior risk" 
rangeaxis <- c(0,1) 
labels<- c("without outcome", "with outcome")
pdf("test4.pdf")
p <- plotPriorPosteriorRisk(data=ExampleData, priorrisk=predRisk1, posteriorrisk=predRisk2, cOutcome=cOutcome, xlabel=xlabel, ylabel=ylabel, rangeaxis=rangeaxis, plotAll=TRUE, plottitle=titleplot, labels=labels)
print(p)
dev.off()

#7.plotRiskDistribution
data(ExampleData) 
cOutcome <- 2
riskmodel <- ExampleModels()$riskModel2
predRisk <- predRisk(riskmodel)
interval <- .05 
xlabel <- "Predicted risk" 
ylabel <- "Percentage" 
xrange <- c(0,1) 
yrange <- c(0,40) 
maintitle <- "Distribution of predicted risks" 
labels <- c("Without outcome", "With outcome")
pdf("test5.pdf")
p <- plotRiskDistribution(data=ExampleData, cOutcome=cOutcome, risks=predRisk, interval=interval, plottitle=maintitle, rangexaxis=xrange, rangeyaxis=yrange, xlabel=xlabel, ylabel=ylabel, labels=labels)
print(p)
dev.off()

#8.plotRiskscorePredrisk
data(ExampleData)
riskmodel <- ExampleModels()$riskModel2
predRisk <- predRisk(riskmodel)
cGenPred <- c(11:16)
riskScore <- riskScore(weights=riskmodel, data=ExampleData, cGenPreds=cGenPred, Type="unweighted")
rangexaxis <- c(0,12) 
rangeyaxis <- c(0,1) 
xlabel <- "Risk score" 
ylabel <- "Predicted risk" 
plottitle <- "Risk score versus predicted risk"
pdf("test6.pdf")
p <- plotRiskscorePredrisk(data=ExampleData, riskScore=riskScore, predRisk=predRisk, plottitle=plottitle, xlabel=xlabel, ylabel=ylabel, rangexaxis=rangexaxis, rangeyaxis=rangeyaxis)
print(p)
dev.off()

#9.plotROC
data(ExampleData) 
cOutcome <- 2
riskmodel1 <- ExampleModels()$riskModel1 
riskmodel2 <- ExampleModels()$riskModel2
predRisk1 <- predRisk(riskmodel1) 
predRisk2 <- predRisk(riskmodel2)
labels <- c("without genetic factors", "with genetic factors")
pdf("test7.pdf")
p <- plotROC(data=ExampleData, cOutcome=cOutcome, predrisk=cbind(predRisk1,predRisk2), labels=labels)
print(p)
dev.off()

#10.predRisk
data(ExampleData) 
cOutcome <- 2 
cID <- 1 
cNonGenPred <- c(3:10) 
cNonGenPredCat <- c(6:8) 
cGenPred <- c(11,13:16) 
cGenPredCat <- c(0)
riskmodel <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome, cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat, cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
predRisk <- predRisk(riskModel=riskmodel)


#11.reclassification
data(ExampleData) 
cOutcome <- 2
riskmodel1 <- ExampleModels()$riskModel1 
riskmodel2 <- ExampleModels()$riskModel2
predRisk1 <- predRisk(riskmodel1) 
predRisk2 <- predRisk(riskmodel2) 
cutoff <- c(0,.10,.30,1)
reclassification(data=ExampleData, cOutcome=cOutcome, predrisk1=predRisk1, predrisk2=predRisk2, cutoff)

#12.riskScore
data(ExampleData) 
cGenPred <- c(11:16)
riskmodel <- ExampleModels()$riskModel2
riskScore <- riskScore(weights=riskmodel, data=ExampleData, cGenPreds=cGenPred, Type="unweighted")

#13.simulatedDataset
ORfreq<-cbind(c(1.35,1.20,1.24,1.16), rep(1,4), c(.41,.29,.28,.51),rep(1,4))
popRisk <- 0.3 
popSize <- 10000
Data <- simulatedDataset(ORfreq=ORfreq, poprisk=popRisk, popsize=popSize)
pdf("test8.pdf")
p <- plotROC(data=Data, cOutcome=4, predrisk=Data[,3])
print(p)
dev.off()

plotDiscriminationBox(data=sub2, cOutcome=cOutcome, predrisk=predRisk, labels=labels)

#准备输入文件---------------------------------------------------------
#Package ‘PredictABEL’
#1.fitLogRegModel
data(ExampleData) 
cOutcome <- 2 
cNonGenPred <- c(3:10) 
cNonGenPredCat <- c(6:8) 
cGenPred <- c(11,13:16) 
cGenPredCat <- c(0)
riskmodel <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome, cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat, cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
summary(riskmodel)

#2.ORmultivariate
data(ExampleData) 
cOutcome <- 2 
cNonGenPred <- c(3:10) 
cNonGenPredCat <- c(6:8)
cGenPred <- c(11,13:16) 
cGenPredCat <- c(0)
riskmodel <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome, cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat, cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
ORmultivariate(riskModel=riskmodel)

#3.plotCalibration
data(ExampleData) 
cOutcome <- 2
riskmodel <- ExampleModels()$riskModel2
predRisk <- predRisk(riskmodel)
rangeaxis <- c(0,1) 
groups <- 20
pdf("test1.pdf")
plotCalibration(data=ExampleData, cOutcome=cOutcome, predRisk=predRisk, groups=groups, rangeaxis=rangeaxis)
print(p)
dev.off()

#4.plotDiscriminationBox
data(ExampleData) 
cOutcome <- 2
riskmodel <- ExampleModels()$riskModel2
predRisk <- predRisk(riskmodel) 
labels <- c("Without disease", "With disease")
pdf("test2.pdf")
plotDiscriminationBox(data=ExampleData, cOutcome=cOutcome, predrisk=predRisk, labels=labels)
print(p)
dev.off()

#5.plotPredictivenessCurve
data(ExampleData)
riskmodel1 <- ExampleModels()$riskModel1 
riskmodel2 <- ExampleModels()$riskModel2
predRisk1 <- predRisk(riskmodel1) 
predRisk2 <- predRisk(riskmodel2)
rangeyaxis <- c(0,1) 
labels <- c("without genetic factors", "with genetic factors")
pdf("test3.pdf")
p <- plotPredictivenessCurve(predrisk=cbind(predRisk1,predRisk2), rangeyaxis=rangeyaxis, labels=labels)
print(p)
dev.off()

#6.plotPriorPosteriorRisk
data(ExampleData) 
cOutcome <- 2
riskmodel1 <- ExampleModels()$riskModel1 
riskmodel2 <- ExampleModels()$riskModel2
predRisk1 <- predRisk(riskmodel1) 
predRisk2 <- predRisk(riskmodel2)
xlabel <- "Prior risk" 
ylabel <- "Posterior risk" 
titleplot <- "Prior versus posterior risk" 
rangeaxis <- c(0,1) 
labels<- c("without outcome", "with outcome")
pdf("test4.pdf")
p <- plotPriorPosteriorRisk(data=ExampleData, priorrisk=predRisk1, posteriorrisk=predRisk2, cOutcome=cOutcome, xlabel=xlabel, ylabel=ylabel, rangeaxis=rangeaxis, plotAll=TRUE, plottitle=titleplot, labels=labels)
print(p)
dev.off()

#7.plotRiskDistribution
data(ExampleData) 
cOutcome <- 2
riskmodel <- ExampleModels()$riskModel2
predRisk <- predRisk(riskmodel)
interval <- .05 
xlabel <- "Predicted risk" 
ylabel <- "Percentage" 
xrange <- c(0,1) 
yrange <- c(0,40) 
maintitle <- "Distribution of predicted risks" 
labels <- c("Without outcome", "With outcome")
pdf("test5.pdf")
p <- plotRiskDistribution(data=ExampleData, cOutcome=cOutcome, risks=predRisk, interval=interval, plottitle=maintitle, rangexaxis=xrange, rangeyaxis=yrange, xlabel=xlabel, ylabel=ylabel, labels=labels)
print(p)
dev.off()

#8.plotRiskscorePredrisk
data(ExampleData)
riskmodel <- ExampleModels()$riskModel2
predRisk <- predRisk(riskmodel)
cGenPred <- c(11:16)
riskScore <- riskScore(weights=riskmodel, data=ExampleData, cGenPreds=cGenPred, Type="unweighted")
rangexaxis <- c(0,12) 
rangeyaxis <- c(0,1) 
xlabel <- "Risk score" 
ylabel <- "Predicted risk" 
plottitle <- "Risk score versus predicted risk"
pdf("test6.pdf")
p <- plotRiskscorePredrisk(data=ExampleData, riskScore=riskScore, predRisk=predRisk, plottitle=plottitle, xlabel=xlabel, ylabel=ylabel, rangexaxis=rangexaxis, rangeyaxis=rangeyaxis)
print(p)
dev.off()

#9.plotROC
data(ExampleData) 
cOutcome <- 2
riskmodel1 <- ExampleModels()$riskModel1 
riskmodel2 <- ExampleModels()$riskModel2
predRisk1 <- predRisk(riskmodel1) 
predRisk2 <- predRisk(riskmodel2)
labels <- c("without genetic factors", "with genetic factors")
pdf("test7.pdf")
p <- plotROC(data=ExampleData, cOutcome=cOutcome, predrisk=cbind(predRisk1,predRisk2), labels=labels)
print(p)
dev.off()

#10.predRisk
data(ExampleData) 
cOutcome <- 2 
cID <- 1 
cNonGenPred <- c(3:10) 
cNonGenPredCat <- c(6:8) 
cGenPred <- c(11,13:16) 
cGenPredCat <- c(0)
riskmodel <- fitLogRegModel(data=ExampleData, cOutcome=cOutcome, cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat, cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
predRisk <- predRisk(riskModel=riskmodel)


#11.reclassification
data(ExampleData) 
cOutcome <- 2
riskmodel1 <- ExampleModels()$riskModel1 
riskmodel2 <- ExampleModels()$riskModel2
predRisk1 <- predRisk(riskmodel1) 
predRisk2 <- predRisk(riskmodel2) 
cutoff <- c(0,.10,.30,1)
reclassification(data=ExampleData, cOutcome=cOutcome, predrisk1=predRisk1, predrisk2=predRisk2, cutoff)

#12.riskScore
data(ExampleData) 
cGenPred <- c(11:16)
riskmodel <- ExampleModels()$riskModel2
riskScore <- riskScore(weights=riskmodel, data=ExampleData, cGenPreds=cGenPred, Type="unweighted")

#13.simulatedDataset
ORfreq<-cbind(c(1.35,1.20,1.24,1.16), rep(1,4), c(.41,.29,.28,.51),rep(1,4))
popRisk <- 0.3 
popSize <- 10000
Data <- simulatedDataset(ORfreq=ORfreq, poprisk=popRisk, popsize=popSize)
pdf("test8.pdf")
p <- plotROC(data=Data, cOutcome=4, predrisk=Data[,3])
print(p)
dev.off()
