rm(list = ls())
pacman::p_load(pacman, haven, tidyr, ggplot2, survival, pseudo, readr, matrixStats,dplyr)
devtools::reload(pkgload::inst("dplyr"))
#### Load Data ####
setwd("~/Library/CloudStorage/Dropbox-UniversityofMichigan/Abby Loe/A. Loe-Project 2/CareQoL")
source("adler_rf.R")
source("careQoL_useful_fxns.R")
load("~/Library/CloudStorage/Dropbox-UniversityofMichigan/Abby Loe/data/imputed_data_Sadness_t_score.Rdata",
     verbose = T)
# this is the sad dataset. 
sad = impute_list[[1]]$complete_dat
plot(sad$Sadness_t_score, sad$Sadness_t_score_BL)
# some EDA
sad = sad %>% 
  mutate(Ri = ifelse(Sadness_t_score>=Sadness_t_score_BL +2,
                     0, 1), # if t-score greater than threshold, sad. Not at risk -> Ri = 0.
         cat_Ri = ifelse(Ri == 0, "Sad", "Happy")) %>% 
  relocate(Ri, .after = Sadness_t_score)

ggplot(sad %>% 
         group_by(subjectid) %>% 
         filter(cur_group_id() %in% c(100, 18)))+
  geom_line(aes(x = week, y = Sadness_t_score, 
                group = subjectid, color = subjectid))+
  scale_color_viridis_d()+
  geom_point(aes(x = week, y = Sadness_t_score,
                 color = cat_Ri))

#### CLDF ####

a = 2; tau = 2; max.follow = 12; 
little.t = seq(from = 0, to = max.follow-tau, by = a);little.t
admin.censor = max(little.t)+tau; final.win = max(little.t)

#### Apply CLDF ####
foo = sad %>%
  filter(week>=0) %>% 
  arrange(subjectid, week) %>%
  group_by(subjectid) %>%
  # Mark days where crisis occurs, else NA
  mutate(sad_week = ifelse(Ri == 0, week, NA_integer_),
         next_sad_week = lead(sad_week)) %>%
  # Fill down to all rows until next crisis using tidyr::fill
  fill(next_sad_week, .direction = "up") %>%
  # Calculate days until next crisis for non-crisis days only
  mutate(week_until_crisis = ifelse(Ri == 1, next_sad_week - week, NA_integer_)) %>%
  ungroup() %>% 
  relocate(next_sad_week, .before = Ri)

shifted.ids = sad %>%
  filter(week == 0 & Ri == 0) %>%
  select(subjectid) %>%
  unique() %>%
  unlist() %>%
  as.vector() # 126 subejcts not at risk at the beginning....
# how many subjects are never at risk??
never.at.risk = sad %>%
  group_by(subjectid) %>% 
  filter(all(Ri == 0, na.rm = T)) %>%
  select(subjectid) %>%
  unique() %>%
  unlist() %>%
  as.vector() # 9 subjects are never at risk...
# should we discount them?


# 
min.day = foo %>%
  filter(subjectid %in% shifted.ids) %>%
  filter(Ri==1) %>%
  group_by(subjectid) %>%
  arrange(week) %>%
  filter(row_number()==1) %>%
  rename(first.day.at.risk = week) %>%
  select(subjectid, first.day.at.risk)

shifted = foo %>% 
  left_join(min.day) %>%
  mutate(first.day.at.risk = ifelse(is.na(first.day.at.risk), 0, first.day.at.risk)) %>% 
  group_by(subjectid) %>% 
  mutate(week = week - first.day.at.risk,
         next_sad_week = next_sad_week - first.day.at.risk,
         Ci = max(week) - first.day.at.risk, 
         censor.time = ifelse(Ci>=max.follow, max.follow, Ci) )%>% 
  # relocate(first.day.at.risk,.after = day) %>% 
  group_by(subjectid) %>% 
  filter(week < censor.time) %>% 
  filter(week %in% little.t) %>% 
  filter(!(subjectid%in% never.at.risk))
table(shifted$week, shifted$Ri)
# table(shifted$day, shifted$delta)
table(shifted$Ci)
table(shifted$censor.time)


primary_events = shifted %>% 
  filter(Ri == 1) %>% 
  group_by(subjectid) %>% 
  ungroup() %>% 
  mutate(min_T_tau = pmin(next_sad_week - week, tau, na.rm = T)) %>% 
  mutate(delta = case_when( next_sad_week<=  censor.time~1, # delta is if next observed event time <= censor time.
                            is.na(next_sad_week) ~0, # if there is no observed next sad week, then delta is 0.
                            next_sad_week> censor.time~0),
   ) %>% 
  ungroup() %>% 
  mutate(min_T_tau = ifelse(delta == 0, pmin(censor.time - week, tau), min_T_tau)) %>% 
  filter(week %in% little.t) %>% 
  relocate(censor.time, next_sad_week, min_T_tau, delta, Sadness_t_score_BL, .after = Ri) 
table(primary_events$delta, primary_events$week)
table(primary_events$Ri, primary_events$week)

dat = NULL
start.time <- Sys.time()
for(tp in (little.t)){
  subdata = primary_events %>% 
    filter(week == tp)
  subdata$pseudoval = pseudomean(subdata$min_T_tau, subdata$delta, tmax = tau)
  dat = rbind.data.frame(dat, subdata)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


full.dat = shifted %>% 
  left_join(dat) %>% 
  filter(week %in% little.t)%>% 
  relocate(censor.time, first.day.at.risk, next_sad_week, min_T_tau, delta, pseudoval, .after = Ri) %>% 
  filter(!(subjectid %in% never.at.risk)) %>% 
  group_by(subjectid) %>% 
  fill(censor.time, .direction = "down")
table(full.dat$delta, full.dat$week)
table(full.dat$Ri, full.dat$week)

full.dat %>% 
  group_by(group) %>% 
  summarise(mean(pseudoval, na.rm = T))

##### DONT MESS WTIH THE ABOVE #####

#### History Variables ####

# I need to get at the history variables for the data
# we don't have the history variables for this group.

hist.frame = sad %>%
  arrange(subjectid, week) %>%
  group_by(subjectid) %>%
  # Mark weeks where crisis occurs, else NA
  mutate(sad_week = ifelse(Ri == 0, week, NA_integer_),
         next_sad_week = lead(sad_week)) %>%
  # Fill down to all rows until next crisis using tidyr::fill
  fill(next_sad_week, .direction = "up") %>%
  # Calculate weeks until next crisis for non-crisis weeks only
  mutate(weeks_until_crisis = ifelse(Ri == 1, next_sad_week - week, NA_integer_)) %>%
  ungroup() %>% 
  relocate(next_sad_week, .before = Ri) %>% 
  filter(week<24) %>% 
  select(subjectid, next_sad_week, week) %>%
  left_join(full.dat %>% select(subjectid, censor.time, first.day.at.risk) %>% 
              unique()) %>% 
  mutate( # leave censor time as it is. We want to get the times ebfore the first.week.
    next_sad_week = next_sad_week - first.day.at.risk
  ) %>% 
  select(-first.day.at.risk) %>% 
  pivot_wider(names_from =week,
              values_from = next_sad_week) %>%
  filter(!(subjectid%in% never.at.risk)) %>% 
  as.matrix() %>% 
  apply(., MARGIN = 1, 
        FUN = function(x) 
          c(x[c(1,2)], 
            as.numeric(setdiff(unique(x[3:length(x)]), NA)), 
            rep(NA_integer_, times = length(x) - length(setdiff(unique(x[3:length(x)]), NA)) - 2))
  ) %>% t()

censor.frame = matrix(rep(as.numeric(hist.frame[,"censor.time"]), times = ncol(hist.frame) - 2), ncol = ncol(hist.frame) - 2, nrow = nrow(hist.frame))
events.only = matrix(as.numeric(hist.frame[,3:ncol(hist.frame)]), ncol = ncol(hist.frame) - 2, nrow = nrow(hist.frame))

events.only[events.only>=censor.frame] = NA
events.only = cbind(events.only, censor.frame[,1]) %>% 
  apply(., MARGIN = 1, FUN = shift_left) %>% t()

history.cal = events.only
history.cal[events.only>=0] = NA_integer_
history.cal = history.cal %>% 
  apply(., MARGIN = 1, FUN = shift_left) %>% t()
history.gap.frame = t(apply(history.cal, 1, get_gap_times))
history.gap.frame[,1] = 0 - history.gap.frame[,1]
history.gap.frame[,1 ] = 23 + history.cal[,1]
init.gap = history.gap.frame[, 1]
init.gap[is.na(init.gap)] = 23
history.gap.frame[,1] = init.gap

Rstar_cap = events.only %>% 
  apply(., MARGIN = 1, 
        FUN = function(x) 
          c( 
            as.numeric(setdiff(unique(x), NA)), 
            rep(NA_integer_, times = length(x) - length(setdiff(unique(x), NA))))
  ) %>% t()
Rstar_cap[Rstar_cap<=0] = NA_integer_
Rstar_cap = t(apply(Rstar_cap, 1, shift_left))
R_cap = t(apply(Rstar_cap, 1, get_gap_times))

flipped.hist = flipped_history_vars(
  R_cap = R_cap,
  Rstar_cap = Rstar_cap,
  little.t = little.t,
  history.gap.frame = history.gap.frame,
  history.cal =  history.cal,
  setting.name = "forward.backward.",
  ID = unique(hist.frame[,"subjectid"])
)


final.formatted.dat = full.dat %>% 
  left_join(flipped.hist, by = c("subjectid" = "ID",
                                 "week" = "time.pt")) %>% 
  group_by(subjectid) %>%
  mutate(across(week, list(
    prior_1 = ~ lag(Ri, 1),
    prior_2 = ~ lag(Ri, 2),
    prior_3 = ~ lag(Ri, 3),
    prior_4 = ~ lag(Ri, 4),
    prior_5 = ~ lag(Ri, 5),
    prior_6 = ~ lag(Ri, 6),
    prior_7 = ~ lag(Ri, 7),
    # prior_8 = ~ lag(Ri, 8),
    post_1 = ~ lead(Ri, 1),
    post_2 = ~ lead(Ri, 2),
    post_3 = ~ lead(Ri, 3),
    post_4 = ~ lead(Ri, 4),
    post_5 = ~ lead(Ri, 5),
    post_6 = ~ lead(Ri, 6),
    post_7 = ~ lead(Ri, 7)#,
    # post_8 = ~ lead(Ri, 8)
  ), .names = "Ri_{fn}")) %>%  
  mutate(ID = cur_group_id(),
         week = as.numeric(week),
         dem_race = case_when((dem_race == "american_indian"|
                                dem_race == "asian"|
                                dem_race == "hawaiian"|
                                dem_race == "more")~"Other",
                              dem_race == "white"~"White",
                              dem_race == "black"~ "Black")) %>% 
  ungroup() 



bad.names =  c( "min_T_tau",
                "next_sad_week",
                "pseudoval", "delta", "censor.time",
                "cat_Ri", "subjectid",
                "Ci",
                "forward.backward.direction.hist",
                "week_until_crisis",
                "first.day.at.risk",
                "Sadness_t_score",
                "cas_relationship",
                "waves", "next_T", "sad_week")




# for(numbucket in c(4, 8, 16, 32)){
#   for(numsamp in c(4, 6, 8, 10, 12)){
#     for (numtrees in c(250, 500, 750)){
for(numbucket in c(4, 8, 16, 32)){
  for(numsamp in c(5)){
    for (numtrees in c(250, 500, 750)){      
      cl= parallel::makeCluster(7)
      registerDoParallel(cl)
      adler.dt = setDT(final.formatted.dat %>%
                         ungroup() %>% 
                         dplyr::select(- any_of(bad.names),
                                       - contains("X"),
                                       - contains("TScore"),
                                       -contains("t_score")
                         ) %>%
                         filter(week >0) %>% 
                         mutate(across(c( contains("Ri_prior"), contains("Ri_post")), ~ifelse(is.na(.), -1, .))) %>% 
                         as.data.frame())
      
      names(adler.dt)
      # 30 possible covariates (including the history covs!)
      table(adler.dt$Ri, adler.dt$week)
      
      split.vars = (c(names(adler.dt %>% dplyr::select(contains("forward.backward"))), 
                      "week"))
      start.time.rf = Sys.time()
      adler.res = calibrated.forest(
        n.iter = numtrees,
        data = adler.dt,
        chosen.seed = 121110,
        time.name = "week",
        cov.mat = adler.dt[, !c("Ri", "ID"), with = FALSE],
        ID.vec = adler.dt$ID,
        nsamp = numsamp,
        minbuck = numbucket,
        always.split = split.vars,
        loss.function = "class")
      end.time.rf = Sys.time()
      stopCluster(cl)
      
      rf.time =end.time.rf-start.time.rf
      print(rf.time)
      adler.res = as.data.frame(adler.res)
      
      final.dat = final.formatted.dat %>% 
        left_join(adler.res, by = c("ID" = "ID", 
                                    "week" = "continuous_time")) %>% 
        mutate(avg.rf.pred = ifelse(week == 0, 1, avg.rf.pred))
      # quantile(1/final.dat$avg.rf.pred[which(final.dat$Ri == 1)])
      # sum(1/final.dat$avg.rf.pred[which(final.dat$Ri == 1)])
      # nrow(final.dat)
      # 
      # ggplot(final.dat %>% filter(Ri == 1))+
      #   geom_histogram(aes(x = 1/avg.rf.pred), bins = 50)+
      #   geom_vline(xintercept = 1)
      # 
      # final.dat %>% 
      #   group_by(group) %>% 
      #   filter(Ri == 1) %>% 
      #   summarise(mean(pseudoval,na.rm = T
      #   ))
      
      # dat = read.csv("sad_dat.csv") %>% 
      # relocate(Sadness_t_score_BL, .after = Sadness_t_score)
      
      write.csv(final.dat, file = paste0("RF_Tuning/sad_dat_thresh2_numbucket", numbucket, "_numsamp", numsamp, "_numtrees", numtrees, ".csv"))
    }
  }
}

