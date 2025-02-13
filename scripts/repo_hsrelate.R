setwd("../output/20240422_hs_scenarios/")

dSeg <- fread("summary_dom_discrete.hs_segr.txt",data.table=F) %>%
    filter(id!="id") %>%
    mutate_at(1:4,as.numeric) %>%
    mutate(par=as.numeric(gsub("(\\d+)(.+$)","\\1",.[[5]])),
           rep=as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",.[[5]])))


parcomb <- fread("parcomb_discrete.csv",data.table=F) %>%
    rename(par=V1) %>%
    mutate(N=paste(N))


dSeg <- inner_join(dSeg,parcomb,by='par')


hsbind2dPlot <- function(data,with.regression=FALSE) {
    predictedGeom <- predictedAnno <-  NULL
    if (with.regression) {
        fit <- nls(model_form,data,start=list(rate=100,int=0.5))
        fit_rate <- round(summary(fit)$coefficients["rate",1],3)
        fit_int <-  round(summary(fit)$coefficients["int",1],3)

        predicted <- data.frame(s=seq(-1,0,length.out=10000)) %>%
            mutate(h=1 / (1/fit_int - fit_rate * s))

        predictedGeom <- geom_line(data=predicted,aes(s,h),linewidth=1)
        predictedAnno <- annotate("label",x=-0.8,y=0.45,
                                  label=paste0("rate: ",fit_rate,", intercept: ",fit_int))
    }
    
    ggplot(data)+
        geom_bin2d(aes(s,h),
                   binwidth=c(0.091,0.046)
                   )+
        predictedGeom+
        ## predictedAnno+
        scale_fill_distiller(direction = 1,
                             palette="YlOrRd",
                             trans="log10",
                             ## breaks = c(1, 10, 30, 100, 300, 1000, 3000, 10000),
                             ## labels = c(1, 10, 30, 100, 300, 1000, 3000, 10000)
                             )+
        coord_cartesian(xlim=c(-1,0),ylim=c(0,0.5))+
        labs(fill='Count of\nalleles')+
        NULL
}


p44_a <- hsbind2dPlot(filter(dSeg,N==20))+
    labs(title="N=20")
p44_b <- hsbind2dPlot(filter(dSeg,N==100))+
    labs(title="N=100")
p44_c <- hsbind2dPlot(filter(dSeg,N==10000),with.regression=T)+
    labs(title="N=10000")

p44 <- (p44_c/p44_b/p44_a)+
    plot_annotation(tag_levels = 'A')
p44

saveplot(p44,"plots/p42_s_vs_h_discrete_bin2d",9,12)
