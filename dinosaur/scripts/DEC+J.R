library(BioGeoBEARS)
library(GenSA)
library(FD)
library(parallel)

#DEC
treefile <- np("Langer.newick")
geogfn <- np("phylogeography.txt")
areas_allowed <- np("geography.txt")
timeperiods <- np("timeperiodsR.txt")
dispersal <- np("dispersal.txt")
distances <- np("distances.txt")

tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size <- 5

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- treefile
BioGeoBEARS_run_object$dispersal_multipliers_fn <- dispersal
BioGeoBEARS_run_object$timesfn <- timeperiods 
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$distsfn <- distances
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE
BioGeoBEARS_run_object$areas_allowed_fn <- areas_allowed 
BioGeoBEARS_run_object$use_optimx <- "GenSA"

BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object <- section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table = TRUE, plot_pieces = FALSE, cut_fossils=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object$num_cores_to_use <- 8

res <- bears_optim_run(BioGeoBEARS_run_object)

#DEC+J
treefile <- np("Langer.newick")
geogfn <- np("phylogeography.txt")
areas_allowed <- np("geography.txt")
timeperiods <- np("timeperiodsR.txt")
dispersal <- np("dispersal.txt")
distances <- np("distances.txt")
dstart <- res$outputs@params_table["d","est"]
estart <- res$outputs@params_table["e", "est"]
jstart <- 0.0001

tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)
max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size <- 5

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- treefile
BioGeoBEARS_run_object$dispersal_multipliers_fn <- dispersal
BioGeoBEARS_run_object$timesfn <- timeperiods 
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$distsfn <- distances
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE
BioGeoBEARS_run_object$areas_allowed_fn <- areas_allowed 
BioGeoBEARS_run_object$use_optimx <- "GenSA"

BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object <- section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table = TRUE, plot_pieces = FALSE, cut_fossils=FALSE)

BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

BioGeoBEARS_run_object$num_cores_to_use <- 8

jres <- bears_optim_run(BioGeoBEARS_run_object)
