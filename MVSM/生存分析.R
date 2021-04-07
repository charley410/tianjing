GBM=read.csv("E:\\бшнд\\GBM.csv",head=T,sep=',', stringsAsFactors=TRUE,fileEncoding = "UTF-8")
#res.cox<-coxph(Surv(time,satute)~new5,data=KRCC)
#sur<-survdiff(Surv(time,satute)~new5,data=KRCC)
fit <- survfit(Surv(time,statue)~Cluster,data=GBM)
#ggsurvplot(fit,data=KRCC,pval=T,risk.table = T,risk.table.col="red",pval.method=T)

ggsurvplot(fit,data=GBM,surv.median.line = "hv",legend.title = " ",legend.labs = c("cluster1", "cluster2", "cluster3"),pval = TRUE,legend = c(0.9,0.9),title ="GBM"(hjust=0.5),xlab = " Time (Days)",censor.shape = 124,censor.size = 2,conf.int = FALSE,ggtheme = theme_bw())
