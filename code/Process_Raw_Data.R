## Set Up Environment

source("~/epistasis_in_breeding/code/Setup_Environment.R")
source("~/epistasis_in_breeding/code/Simulation_Parameters.R")
source("~/epistasis_in_breeding/code/Functions.R")

focal_qtl = 1

## Small Formats to Raw Data

long_trait_results <- trait_results[-1,] %>% rename(.,Rep=rep) %>% rename(.,Cycle=cycle) %>% rename(.,Trait=trait) %>%
  rename(.,Per_Rep_GMean=GV_Mean)

population_data <- as.tibble(matrix(NA,nrow=(n_qtl*n_chr),ncol = ncol(population_results[[1]])))
names(population_data) <- names(population_results[[1]])

for (rep in 1:n_reps){
  rep_population_results <- population_results[[rep]] %>% as.tibble()
  population_data <- rbind(population_data,rep_population_results)
}

population_data <- population_data[-c(1:(n_qtl*n_chr)),]
rep_id <- rep(c(1:n_reps),each = (n_qtl*n_chr))

popgen_data <- population_data %>% bind_cols(rep_id,.) %>% rename(Rep = ...1) %>% select(-contains("trait"))
names(popgen_data) <- sub("generation", "", colnames(popgen_data))
long_population_data <- gather(popgen_data, key = "Cycle", value = "Estimates",-c(Rep,Primary_QTL)) %>%
  unnest(cols = c(Estimates)) %>%
  group_by(Primary_QTL) %>%
  mutate(QTL_Category = factor(ifelse(Primary_QTL %in% focal_qtl, "Epistasis", "No Epistasis")))

ase_data <- population_data %>% bind_cols(rep_id,.) %>% rename(Rep = ...1) %>% select(-contains("generation"))

saveRDS(long_trait_results,file=paste0("./data/GeneticTrend_Trait",trait,".rds"))
saveRDS(long_population_data,file=paste0("./data/PopData_Trait",trait,".rds"))
saveRDS(ase_data,file=paste0("./data/ASE_Trait",trait,".rds"))

### Process Raw Data for Plotting ----

### Extract Per Trait, Per Selection Cycle, ASE and convert in to a long dataframe format ----
trait_ase_data <- ase_data %>% gather(key = "Trait",value = "Data",-c(Rep,Primary_QTL)) %>%
  unnest(Data) %>%
  select(-contains("effects")) %>%
  gather(key = "Cycle",value = "Estimate",-c(Rep,Primary_QTL,Trait)) %>%
  unnest(Estimate) %>%
  select(Rep,Primary_QTL,Trait,Cycle,ASE,SE_ASE,Rhat_per_locus) %>%
  mutate(Trait = gsub("trait","Trait ",Trait),Cycle = gsub("cycle","",Cycle)) %>%
  group_by(Primary_QTL) %>%
  mutate(QTL_Category = factor(ifelse(Primary_QTL %in% focal_qtl, "Epistasis", "No Epistasis")))

### Extract ASE's in the 1st Selection Cycle ----
initialAseData <- trait_ase_data %>%
  group_by(Trait,Rep) %>%
  filter(Cycle == 1) %>%
  select(Rep,Trait,ASE)

# Extract ASE for Trait 1 in the 1st Selection Cycle
initAseTrait1 <- initialAseData %>%
  filter(Trait=='Trait 1') %>%
  map_df(rev) %>%
  mutate(., QTL_Value_Incr = (0.8*ASE)*2) %>% #Predict Increase in Trait 1 Value at each QTL based on ASE
  mutate(., Cycle = rep(c(2:11),times=5))

# Extract ASE for Trait 2 in the 1st Selection Cycle
initAseTrait2 <- initialAseData %>%
  filter(Trait=='Trait 2') %>%
  map_df(rev) %>%
  mutate(., QTL_Value_Incr = (0.8*ASE)*2) %>% #Predict Increase in Trait 2 Value at each QTL based on ASE
  mutate(., Cycle = rep(c(2:11),times=5))

# Extract ASE for Trait 3 in the 1st Selection Cycle
initAseTrait3 <- initialAseData %>%
  filter(Trait=='Trait 3') %>%
  map_df(rev) %>%
  mutate(., QTL_Value_Incr = (0.8*ASE)*2) %>% #Predict Increase in Trait 3 Value at each QTL based on ASE
  mutate(., Cycle = rep(c(2:11),times=5))

#Recombine into one dataframe, convert Trait column into numbers
initAseData <- bind_rows(initAseTrait1,initAseTrait2,initAseTrait3) %>% mutate(Trait = gsub("Trait ","",Trait),Trait = as.numeric(Trait))

### Calculate Genetic Correlations between Traits in 1st Selection Cycle  ----
# cor.test(initASETrait1$ASE,initASETrait2$ASE) #Trait 1 & 2
cor.test(initAseTrait1$ASE,initAseTrait3$ASE) #Trait 1 & 3
# cor.test(initASETrait2$ASE,initASETrait3$ASE) #Trait 2 & 3

### Calculate Mean and Variance of ASE Over Selection Cycles ----
ase_plot_data <- trait_ase_data %>%
  group_by(Primary_QTL) %>%
  mutate(QTL_Category = factor(ifelse(Primary_QTL %in% focal_qtl, "Epistasis", "No Epistasis"))) %>%
  ungroup() %>%
  group_by(Trait,Cycle,QTL_Category) %>%
  mutate(Avg_ASE = mean(ASE),ASE_Variance = var(ASE)) %>%
  ungroup() %>%
  rename(.,QTL = Primary_QTL) %>%
  rename(.,`QTL Effects` = QTL_Category)

saveRDS(initAseData,file=paste0("./output/1stCycle_ASE_Data.rds"))
saveRDS(ase_plot_data,file=paste0("./output/ASE_Trajectory_Plot_Data.rds"))

# Combine Genetic Mean and Variance Results with ASE Estimated in 1st Selection Cycle

trait_results <- full_join(long_trait_results, initAseData, by=c("Rep","Trait","Cycle")) %>%
  group_by(Rep,Trait,Cycle) %>%
  mutate(QTL_Value_Incr = replace_na(QTL_Value_Incr, Per_Rep_GMean)) %>%
  group_by(Rep,Trait) %>%
  mutate(Projected_GV_Mean = cumsum(QTL_Value_Incr)) %>%
  ungroup() %>%
  melt(., id.vars = c("Rep","Trait","Cycle","Per_Rep_GMean","GV_Variance"), measure.vars = c("Per_Rep_GMean","Projected_GV_Mean"),var = "Selection") %>% #Merge True and Predicted mean genetic value into 1 column, and add "Selection" as a factor column
  rename(Across_Rep_GMean = value) %>%
  unique() %>% #remove duplicate entries
  mutate(Per_Rep_GMean = replace(Per_Rep_GMean,Selection == "Predicted Values",NA)) %>%
  group_by(Trait,Cycle,Selection) %>%
  mutate(Overall_Gv_Mean = mean(Per_Rep_GMean),Overall_Gv_Variance = mean(GV_Variance),Overall_Genetic_Mean = mean(Across_Rep_GMean)) %>%
  ungroup() %>%
  select("Rep","Cycle","Trait","Selection","Per_Rep_GMean","Overall_Gv_Mean","Across_Rep_GMean","Overall_Genetic_Mean")
trait_results$Selection <- trait_results$Selection %>% str_replace_all(c("Projected_GV_Mean" = "Predicted Values","Per_Rep_GMean" = "True Values"))

saveRDS(trait_results,file=paste0("./output/Trait_Trajectory_Plot_Data.rds"))

all_freq_plot_data <- long_population_data %>%
  group_by(Cycle,QTL_Category) %>%
  mutate(avg_all_freq = mean(p),all_freq_variance = var(p)) %>%
  rename(.,QTL = Primary_QTL) %>%
  rename(.,`QTL Effects` = QTL_Category) %>%
  ungroup()

saveRDS(all_freq_plot_data,file=paste0("./output/QTL_Trajectory_Plot_Data.rds"))


