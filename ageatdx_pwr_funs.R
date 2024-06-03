library(tidyverse)
library(faux)
library(tidymodels)


# Set.seed for reproducibility
set.seed(12309)

workhorse_sim = function(n=1400, #Sample size
                         edu=-0.2, #This is the standardised effect of confounder education on age at dx
                         fun=0.4, #This is the standardised effect of confounder functioning on age at dx
                         acc=-0.1, #This is the standardised effect of confounder access to services on age at dx
                         paternal_eff = -0.1, #This is the standardised effect of age_at_dx on paternal wellbeing
                         maternal_eff = -0.2, #This is the standardised effect of age_at_dx on maternal wellbeing
                         sib_eff = -0.05, #This is the standardised effect of age_at_dx on sibling wellbeing
                         ind_eff = -0.2, #This is the standardised effect of age_at_dx on individual wellbeing
                         htest_level="family", #How is H3 specified? Alternatives are "person" and "outcome"
                         n_replicates=100){ #How many times should the scenario be simulated
  
  allres=list()
  for(i in 1:n_replicates){   
    
    if(i %in% n_replicates/10){
    message(paste0("running replicate number ", i, " of ", n_replicates))
    }
    # Simulate putative confounders
    
    simdat = faux::rnorm_multi(n=n,
                               mu = c(educ = 0, funct = 0, access = 0),
                               sd = c(1,1,1),
                               r = c( 1,   -0.2, 0.5, 
                                      -0.2, 1,    0, 
                                      0.5, 0,  1))
    
    # Simulate age at diagnosis with a relationship to confounders, removing low values and rounding
    
    simdat = simdat %>% 
      mutate(age_at_dx = 4 + edu*educ + fun*funct + acc*access + rnorm(nrow(simdat),0, 0.5 )) %>% 
      mutate(age_at_dx = ifelse(age_at_dx<2,2,round(age_at_dx,0)))
    
    # Check relationships
    
    lm(age_at_dx ~ educ + funct + access, data=simdat) %>% broom::tidy()
    
    # Simulate outcomes 
    
    simdat = simdat %>% 
      mutate(p_wellb_1 = rnorm(1,paternal_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             p_wellb_2 = rnorm(1,paternal_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             p_wellb_3 = rnorm(1,paternal_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             m_wellb_1 = rnorm(1,maternal_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             m_wellb_2 = rnorm(1,maternal_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             m_wellb_3 = rnorm(1,maternal_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             s_wellb_1 = rnorm(1,sib_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             s_wellb_2 = rnorm(1,sib_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             s_wellb_3 = rnorm(1,sib_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             i_wellb_1 = rnorm(1,ind_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             i_wellb_2 = rnorm(1,ind_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ),
             i_wellb_3 = rnorm(1,ind_eff,0.05)*age_at_dx + 0.5*educ + -0.2*funct + 0.1*access + rnorm(nrow(simdat),0, 0.5 ))
    
    # Run models
    
    model_func = function(yvar, dat){
      model = paste0(yvar," ~ age_at_dx + educ + funct + access")
      tmp = lm(model, data= dat )%>% broom::tidy()
      return(tmp)
      
    }
    
    linear_mods_res <- 
      names(simdat %>% select(matches("wellb"))) %>%
      as_tibble() %>%
      select("yvar"=value) %>% 
      mutate(results = pmap(., model_func, dat=simdat))
    
    # Perform FDR correction at specified level
    
    if(htest_level=="family"){
      
      res = linear_mods_res$results %>% purrr::reduce(bind_rows) %>% 
        filter(term=="age_at_dx") %>% 
        mutate(outcome= rep(c(1,2,3),4),
               person= rep(c("father","mother","sib","self"), each=3),
               fdrp=p.adjust(p.value, method="fdr"),
               mtc_sig= ifelse(fdrp<0.05, 1,0))
      
    }else if(htest_level=="person"){
      
      res = linear_mods_res$results %>% purrr::reduce(bind_rows) %>% 
        filter(term=="age_at_dx") %>% 
        mutate(outcome= rep(c(1,2,3),4),
               person= rep(c("father","mother","sib","self"), each=3)) %>% 
        group_by(person) %>% 
        mutate(fdrp=p.adjust(p.value, method="fdr"),
               mtc_sig= ifelse(fdrp<0.05, 1,0))%>% 
        ungroup()
      
    }else if(htest_level=="outcome"){
      
      res = linear_mods_res$results %>% purrr::reduce(bind_rows) %>% 
        filter(term=="age_at_dx") %>% 
        mutate(outcome= rep(c(1,2,3),4),
               person= rep(c("father","mother","sib","self"), each=3)) %>% 
        group_by(person,outcome) %>% 
        mutate(fdrp=p.adjust(p.value, method="fdr"),
               mtc_sig= ifelse(fdrp<0.05, 1,0)) %>% 
        ungroup()
      
    }
    
    #Return annotated res 
    allres[[i]]=res %>% 
      mutate(repl=i,
             htest_level=htest_level,
             paternal_eff=paternal_eff,
             maternal_eff= maternal_eff,
             sib_eff=sib_eff,
             ind_eff=ind_eff)
  }
  
  allres= allres %>% 
    purrr::reduce(bind_rows)
  return(allres)
}

pwr_sim = function(n=1400, #Sample size
                   edu=-0.2, #This is the standardised effect of confounder education on age at dx
                   fun=0.4, #This is the standardised effect of confounder functioning on age at dx
                   acc=-0.1, #This is the standardised effect of confounder access to services on age at dx
                   paternal_eff = -0.1, #This is the standardised effect of age_at_dx on paternal wellbeing
                   maternal_eff = -0.2, #This is the standardised effect of age_at_dx on maternal wellbeing
                   sib_eff = -0.05, #This is the standardised effect of age_at_dx on sibling wellbeing
                   ind_eff = -0.2, #This is the standardised effect of age_at_dx on individual wellbeing
                   htest_level="family", #How is H3 specified? Alternatives are "person" and "outcome"
                   n_replicates=100){
  
  # Expand inputs (this allows for concatenated vectors to be supplied for all inputs)
  inputs=expand.grid(n,edu,fun,acc,paternal_eff,maternal_eff,sib_eff,ind_eff,htest_level) %>% 
    `colnames<-`(c("n","edu","fun","acc","paternal_eff","maternal_eff","sib_eff","ind_eff",
                 "htest_level"))
  # Run the simulations
  out= inputs %>% 
    as_tibble %>% 
    mutate(results= pmap(., workhorse_sim, n_replicates=n_replicates, .progress=T)) 

  #Summarise power
  
  out= out %>% 
    mutate(overall_power = purrr::reduce(purrr::map(out$results, function(x){
      x %>% 
        summarise(power=sum(mtc_sig/n())) %>% 
        ungroup()
      
    }), bind_rows), 
    person_power = purrr::map(out$results, function(x){
      x %>% 
        group_by(person) %>% 
        summarise(power=sum(mtc_sig/n()))%>% 
        ungroup()
      
    }), 
    outcome_power = purrr::map(out$results, function(x){
      x %>% 
        group_by(person,outcome) %>% 
        summarise(power=sum(mtc_sig/n()))%>% 
        ungroup()
      
    })) 
    
  
  return(out)
}



  