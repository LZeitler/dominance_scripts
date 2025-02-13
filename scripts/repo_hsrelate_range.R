setwd("../output/range/")
library(patchwork)

## FIG2

d0 <- fread("hs-range_summary.txt",data.table=F) %>%
    rename(id=1,ac=2,s=3,h=4,genorig=5,file=6) %>%
    mutate(stage=substr(file,1,4),
           par=as.numeric(gsub("(.+par)(\\d+)(.+$)","\\2",file)),
           rep=as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",file))) %>%
    select(-file) %>%
    mutate(popOrigin=ifelse(genorig>=2000,"Core","Expansion")) %>% 
    group_by(stage,par,rep,popOrigin) %>%
    mutate(approxgen=max(genorig)) %>%
    group_by(stage,par,rep) %>%
    mutate(approxgen=min(approxgen)) %>%
    ungroup %>%
    mutate(mutAge=ifelse(popOrigin=="Core",approxgen-genorig+100000,approxgen-genorig),
           isSingleton=ifelse(ac==1,"Singleton","No singleton"))

d <- d0 %>% mutate(h=h+0.001)

stagesAssigned <- c("t1_A", "t1_B", "t3_A", "t3_B", "t3_C", "t3_Z", "t4_C", "t4_Z" )
stagesnames <-
    c("Interior 1, t=1","Interior Edge, t=1","Interior","Interior 2",
      "Edge","Core","Edge t+1000","Core t+1000")
names(stagesAssigned) <- stagesnames
d <- d %>%         mutate(stage_names=forcats::fct_recode(stage,!!!stagesAssigned))


## vertical density chart - histogram-like
hsHistPlot <- function(data,additional.s.filters=c(-1,-0.10009,-0.00010)) {
    data <- data %>%
        group_by(stage_names,h,s) %>%
        count %>%
        group_by(stage_names) %>%
        mutate(relCount=n/sum(n)) %>%
        ungroup %>%
        filter(s%in%additional.s.filters)
    ggplot(data)+
        geom_histogram(aes(paste(h-0.001),relCount,fill=stage_names),
                       stat='identity',alpha=0.7,
                       position=position_dodge())+
        geom_hline(aes(yintercept=1/121,color="input/\nde novo"))+
        ## scale_y_log10()+
        scale_x_discrete(guide=guide_axis(angle = 45))+
        ## scale_y_sqrt(labels = scales::percent)+
        scale_y_continuous(labels = scales::percent)+
        scale_fill_manual(values=c("#417A00","#1F0F87"))+
        scale_color_manual(values="red")+
        facet_grid(forcats::fct_relabel(factor(s), ~ paste0("s = ", .x))~.,
                   switch="y",
                   scales="free")+
        labs(x="Dominance coefficient (h)",
             y="Percentage counts of mutations",
             )+
        theme(legend.title=element_blank())
}


p_hist1 <- hsHistPlot(filter(d,stage%in%c("t3_C","t3_Z"))
                      )


p_hist2 <- hsHistPlot(filter(d,stage%in%c("t3_C","t3_Z"),
                             mutAge==0)
                      )


p_hist3 <- hsHistPlot(filter(d,stage%in%c("t3_C","t3_Z"),
                             isSingleton=="Singleton"
                      ))

p_hist_all <- (p_hist1+
               theme(legend.position = 'none',
                     plot.tag.position=c(0,1))+
               labs(title="All Mutations"))+
    (p_hist2+theme(legend.position = 'none',
                   axis.title.y=element_blank(),
                   plot.tag.position=c(-0.05,1))+
     labs(title="De novo"))+
    (p_hist3+theme(axis.title.y=element_blank(),
                   plot.tag.position=c(-0.05,1))+
     labs(title="Singletons"))
p_hist_all <- p_hist_all+
    plot_annotation(tag_levels = 'A')
p_hist_all


saveplot(p_hist_all,"plots/histgram_svsh_corevsedge_freescales_all",12,9)



### FIG3
#### SHARED ALLELES BETWEEN CORE, INTERIOR AND EDGE
ds <- d %>%
    filter(stage%in%c("t4_Z","t3_Z","t3_B","t3_C","t4_C")) %>% # remove uninteresting stages
    mutate(stage_comparisons=ifelse(stage%in%c("t3_C","t4_C"),"Edge","Not edge")) %>%
    mutate(tp=ifelse(stage%in%c("t4_Z","t4_C"),"expanded+1000","expanded")) %>%
    select(id,s,h,genorig,par,rep,stage_comparisons,tp) %>% 
    group_by(par,rep,stage_comparisons,tp) %>%
    unique() %>% 
    group_by(par,rep,tp) %>%
    mutate(dup=duplicated(id)|
               duplicated(id,fromLast = T))  %>%
    ## filter(dup==T) %>%
    group_by(par,rep,id,h,s,tp) %>%
    count %>%
    group_by(par,rep,id,h,s) %>% 
    select(-tp) %>%
    unique %>% 
    ungroup

dsrel <- ds %>% mutate(private=ifelse(n==1,"Private","Shared")) %>%
    group_by(h,s,private) %>%
    summarise(countShared=n()) %>%
    group_by(h) %>%
    mutate(pcountShared=countShared/sum(countShared))

dsrel %>% filter(private!="Private",countShared>0) # number

nbrs <- dsrel %>% group_by(private) %>%  summarise(sum(countShared)) # numbers
nbrs;nbrs[2,2]/sum(nbrs[,2])

ps6_a_shared <- ds %>%
    filter(n>1) %>% # shared
    ggplot()+
    geom_bar(aes(x=h,fill=factor(s)),width=0.04)+
    scale_fill_viridis_d(end=0.9,direction=-1)+
    labs(x="Dominance coefficient (h)",y="Shared allele count",fill="s",
         title="Shared Mutations (absolute)")
ps6_a_shared

ps6_a_private <- ds %>%
    filter(n==1) %>% # private
    ggplot()+
    geom_bar(aes(x=h,fill=factor(s)),width=0.04)+
    scale_fill_viridis_d(end=0.9,direction=-1)+
    labs(x="Dominance coefficient (h)",y="Private allele count",fill="s",
         title="Private Mutations (absolute)")
ps6_a_private


ps6_p_shared <- dsrel %>%
    filter(private=="Shared") %>%
    ggplot()+
    geom_col(aes(x=h,y=pcountShared,fill=factor(s)))+
    scale_y_continuous(labels=scales::percent)+
    scale_fill_viridis_d(end=0.9,direction=-1)+
    labs(x="Dominance coefficient (h)",y="Percent shared alleles",fill="s",
         title="Shared Mutations (%)")
ps6_p_shared

ps6_p_private <- dsrel %>%
    filter(private!="Shared") %>%
    ggplot()+
    geom_col(aes(x=h,y=pcountShared,fill=factor(s)))+
    scale_y_continuous(labels=scales::percent)+
    scale_fill_viridis_d(end=0.9,direction=-1)+
    labs(x="Dominance coefficient (h)",y="Percent private alleles",fill="s",
         title="Private Mutations (%)")
ps6_p_private

ps6 <- ((ps6_p_shared+theme(legend.position='none'))+ps6_p_private)/
    ((ps6_a_shared+theme(legend.position='none'))+ps6_a_private)+
    plot_annotation(tag_levels = "A")+
    plot_layout(guides = 'collect')
ps6

saveplot(ps6,"plots/ps6_sharedallele_percentandabs_sharedandprivate",12,9)




## FIG S3
## Helper functions
hsbind2dPlot <- function(data,val.only=FALSE) {
    if(val.only) {
        p <- list(fit_int,fit_rate)
    } else {
        p <- ggplot(data)+
            geom_bin2d(aes(s,h),
                       binwidth=c(0.091,0.046)
                       )+
            scale_fill_distiller(direction = 1,
                                 palette="YlOrRd",
                                 trans="log10",
                                 )+
            coord_cartesian(xlim=c(-1,0),ylim=c(0,0.5))+
            NULL
    }
    return(p)
}
prepCountData <- function(data) {
    ## prepare / uncount mutation data, after filtering,
    ## results in a simple but long list of mutations
    t <- data %>%
        unnest(counts) %>%
        ungroup %>%
        select(h,s,z) %>% 
        uncount(z)
    gc()
    return(t)
}


sbins <- seq(-1,-0.0001,length=11)
c <- 1
muts <- data.frame()
for (s in seq_along(sbins)) {
    for (i in 0:10) {
        muts <- bind_rows(muts,data.frame(s=sbins[s],h=i*0.1/2,mtype=c))
        c <- c+1
    }
}
c <- c-1 # number of mutation types


r <- fread("fitness-stats_summary.txt",data.table=F,sep="\t")
r$rep <- as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",r$V8))
r$countsAnc <- lapply(strsplit(r$V6, ","), function(x) {
    bind_cols(muts,z=as.numeric(x))
})
r$countsNew <- lapply(strsplit(r$V7, ","), function(x) {
    bind_cols(muts,z=as.numeric(x))
})

rl <- r %>%
    select(-V6,-V7,-V8,gen=V1,deme=V2,fit=V3,N=V4,het=V5,countsAnc,countsNew) %>%
    group_by(across(!c(countsAnc,countsNew))) %>%
    summarize(tibble(countsAnc),tibble(countsNew))



rlEveryGenMinimal <- rl %>% group_by(gen,rep) %>% filter(N>10)
rlEveryGenStaticDemes <- rlEveryGenMinimal %>%
    mutate(isEdge=
               ifelse(deme==1,"Core",
               ifelse(deme==25,"Deme 25",
               ifelse(deme==35,"Deme 35",
               ifelse(deme==45,"Deme 45",NA))))) %>%
    filter(!is.na(isEdge)) %>% ungroup()
rlEveryGenEdgeDemes <- rlEveryGenMinimal %>%
    filter(deme==max(deme)) %>%
    mutate(isEdge="Edge") %>% ungroup()
rlEveryGen <- bind_rows(
    rlEveryGenStaticDemes,
    rlEveryGenEdgeDemes
) %>%
    mutate(isEdge=factor(isEdge,
                         levels=c("Core","Deme 25", "Deme 35", "Deme 45", "Edge"),
                         ordered=T)) %>% 
    group_by(rep,gen) %>%
    mutate(demeReached=max(deme)) %>%
    group_by(rep)

## Core at the end of sim
rlf <- rlEveryGen %>% filter(isEdge=="Core",gen==max(gen))
rlfm <- rlf %>%
    select(everything(),counts=countsNew,-countsAnc) %>%
    prepCountData

p2r <- hsbind2dPlot(rlfm)+labs(fill="Count")
p2r

saveplot(p2r,"plots/p2_range_newCore_ms")


## FIG S4
rlfm <- rlf  %>% 
    select(everything(),counts=countsNew,-countsAnc) %>%
    prepCountData
p1r <- hsbind2dPlot(rlfm)+labs(fill="Count")
p1r

saveplot(p1r,"plots/p1_range_newEdge_ms")

