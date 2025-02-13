source("~/code/r/source_me_light.R")

## FIG S1

setwd("../output/20240212_hs_scenarios/")

dSeg <- fread("summary_dom.hs_segr.txt",data.table=F) %>%
    filter(id!="id") %>%
    mutate_at(1:5,as.numeric) %>%
    mutate(par=as.numeric(gsub("(\\d+)(.+$)","\\1",.[[6]])),
           rep=as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",.[[6]])))

dFix <- fread("summary_dom.hs_fixed.txt",data.table=F) %>%
    filter(no!="no") %>%
    select(-type,-pop) %>% 
    mutate_at(1:9,as.numeric) %>%
    mutate(par=as.numeric(gsub("(\\d+)(.+$)","\\1",.[[10]])),
           rep=as.numeric(gsub("(.+rep)(\\d+)(.+$)","\\2",.[[10]])))

parcomb <- fread("parcomb.csv",data.table=F) %>%
    rename(par=V1) %>%
    mutate(N=paste(N),
           h_shape=gsub("[^A-Za-z0-9]", "", h_shape)
           )

dSeg <- inner_join(dSeg,parcomb,by='par') %>%
    filter(!(h_shape=="fixeds" & s!=-0.05)) # these are dummy mutations (no selection)
dFix <- inner_join(dFix,parcomb,by='par') %>%
    filter(!(h_shape=="fixeds" & s!=-0.05)) # these are dummy mutations (no selection)



p61 <- plot_grid(
    dSeg %>%
    filter(h_shape=="inverse"
           ) %>%
    ggplot()+
    geom_density(
        aes(x=h,fill=N, y = after_stat(scaled)),
        alpha=.5)+
    geom_function(
        fun = h_density,
        args = list(int = 0.5, rate = 100, smean = 0.05, scaled = T),
        n = 1000,
        color = "grey50", linewidth = 2)+
    xlim(0,0.5)+
    labs(y="scaled density")+
    background_grid(minor='y'),
    dSeg %>%
    filter(h_shape=="unif") %>%
    ggplot()+
    geom_density(
        aes(x=h,fill=N, y = after_stat(scaled)),
        alpha=.5)+
    xlim(0,0.5)+
    labs(y="scaled density")+
    background_grid(minor='y'),
    ncol=1,
    labels="AUTO"
    )
p61

saveplot(p61,"plots/p61_obs_h_dist_N_input_ms",8,8)


## FIG S2

p81 <- plot_grid(
    dFix %>%
    filter(h_shape=="inverse"
           ) %>%
    ggplot()+
    geom_density(
        aes(x=h,fill=N, y = after_stat(scaled)),
        alpha=.5)+
    geom_function(
        fun = h_density,
        args = list(int = 0.5, rate = 100, smean = 0.05, scaled = T),
        n = 1000,
        color = "grey50", linewidth = 2)+
    xlim(0,0.5)+
    labs(y="scaled density")+
    background_grid(minor='y'),
    dFix %>%
    filter(h_shape=="unif") %>%
    ggplot()+
    geom_density(
        aes(x=h,fill=N, y = after_stat(scaled)),
        alpha=.5)+
    xlim(0,0.5)+
    labs(y="scaled density")+
    background_grid(minor='y'),
    ncol=1,
    labels="AUTO"
    )
p81

saveplot(p81,"plots/p81_obs_h_density_fixed_N_input_ms",8,8)
