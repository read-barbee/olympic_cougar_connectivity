vars <- c("tri + f(id1, tri, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "tpi + f(id2, tpi, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "ndvi + f(id3, ndvi, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "urban + f(id4, urban, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "crop + f(id5, crop, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "grass + f(id6, grass, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "shrub + f(id7, shrub, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "wetland + f(id8, wetland, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "popbig + f(id9, popbig, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "riv_width + f(id10, riv_width, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))",
          "roads_int + f(id11, roads_int, values=1:74,model='iid',hyper=list(theta=list(initial=log(1),fixed=F,prior='pc.prec',param=c(3,0.05))))")

all_comb <- do.call("c", lapply(seq_along(vars), function(i) combn(vars, i, FUN = list)))
response <- "y"
mods <- list()
dic_list <- list()
## Took 32273 seconds (~9hrs) on Legion for 11 covariates (2047 combinations)
system.time(
  for(i in 1:length(all_comb)){
    var_i <- all_comb[[i]]
    form <- as.formula(paste(response, paste("-1", "f(strat, model='iid', hyper=list(theta = list(initial=log(1e-6),fixed=T)))", paste(var_i, collapse="+"), sep="+"), sep="~"))
    
    mods[[i]] <- inla(form, family ="Poisson", data=disp_inla,
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)), control.compute = list(dic = TRUE))
    
    dic_list[[i]] <- mods[[i]]$dic$dic
    
    print(i/length(all_comb))
  }
)