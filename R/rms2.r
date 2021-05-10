#' Reverse 2-sample MR analysis
#' 
#' @description
#' This analysis is a wrapper around MR and colocalisation functions for systematically searching OpenGWAS for mediating pathways from genotype to trait
#' 
#' @export
rms2 <- R6::R6Class("rms2", list(
	#' @field igd_id The OpenGWAS dataset to analyse as the outcome GWAS
	igd_id = NULL,
	#' @field info Meta data for igd_id
	info = NULL,
	#' @field gwashits Data from of GWAS hits
	gwashits = NULL,
	#' @field rsid_scan List of results of PheWAS search - one item for each GWAS hit
	rsid_scan = list(),
	#' @field coloc_data Data used for colocalisation analysis - list of lists (rsid x candidate trait)
	coloc_data = list(),
	#' @field coloc_result Results from colocalisation analysis - list of data frames, one data frame for each region analysed
	coloc_result = list(),
	#' @field candidate_instruments List of instrument datasets for each candidate
	candidate_instruments = list(),
	#' @field candidate_instruments_outcome Dataframe of instrument-outcome associations
	candidate_instruments_outcome = NULL,
	#' @field mr_scan Results from MR analysis - list of data frames with one data frame for each region
	mr_scan = list(),
	#' @field mv_dat Data and results from mvmr
	mv_dat = list(),


	#' @description
	#' Initialise the rms2 class for a particular outcome GWAS
	#'
	#' @param igd_id The OpenGWAS ID to use
	initialize = function(igd_id)
	{
		gi <- ieugwasr::gwasinfo(igd_id)
		if(nrow(gi) != 1)
		{
			stop(igd_id, " not found in OpenGWAS database")
		}
		message(str(gi))
		self$igd_id <- igd_id
		self$info <- gi
	},


	#' @description
	#' Extract GWAS hits for igd_id
	#'
	#' @param igd_id Default=self$igd_id
	#' @param pval Significance threshold for detecting GWAS hit default=5e-8
	#' @param kb Clumping parameter default=10000
	#' @param r2 Clumping parameter default=0.001
	#' @param pop 1000 genome super population to use for clumping default='EUR'
	extract_gwashits = function(igd_id=self$igd_id, pval=5e-8, kb=10000, r2=0.001, pop='EUR')
	{
		self$gwashits <- ieugwasr::tophits(igd_id, pval=pval, r2=r2, kb=kb, pop=pop)
		message("Found ", nrow(self$gwashits), " GWAS hits for ", igd_id)
	},


	#' @description
	#' Search OpenGWAS for associations with a particular GWAS hit
	#'
	#' @param rsid RSID
	#' @param pval Significance threshold for candidate traits default=5e-8
	#' @param idlist Which traits to search. When NULL (default) will search all available traits
	scan_rsid = function(rsid, pval=5e-8, idlist=NULL)
	{
		if(is.null(idlist))
		{
			res <- ieugwasr::phewas(rsid, pval=pval)
		} else {
			res <- ieugwasr::associations(rsid, idlist) %>%
				dplyr::filter(p < pval)
		}
		res <- subset(res, !id == self$igd_id)
		for(r in rsid)
		{
			message("Found ", nrow(res), " associations for ", r)
			self$rsid_scan[[r]] <- subset(res, rsid == r)
		}
	},

	#' @description
	#' Perform colocalisation between outcome trait and a candidate trait for a particular locus
	#'
	#' @param rsid RSID
	#' @param trait_id trait_id
	#' @param radius Radius around which to perform analysis default=50000
	colocalise_rsid = function(rsid, trait_id, radius=50000)
	{
		if(! rsid %in% names(self$coloc_data))
		{
			self$coloc_data[[rsid]] <- list()
		}
		chrpos <- paste0(self$gwashits$chr[self$gwashits$rsid == rsid], ":", self$gwashits$position[self$gwashits$rsid == rsid] - radius, "-", self$gwashits$position[self$gwashits$rsid == rsid] + radius)
		dat <- gwasglue::ieugwasr_to_coloc(self$igd_id, trait_id, chrpos)
		res <- coloc::coloc.abf(dat[[1]], dat[[2]])
		self$coloc_data[[rsid]][[trait_id]] <- list(data=dat, coloc=res)
	},

	#' @description
	#' Plot region for a colocalisation analysis
	#'
	#' @param rsid rsid
	#' @param trait_id trait_id
	plot_coloc = function(rsid, trait_id)
	{
		temp <- gwasglue::coloc_to_gassocplot(self$coloc_data[[rsid]][[trait_id]][["data"]])
		require(gassocplot)
		gassocplot::stack_assoc_plot(temp$markers, temp$z, temp$corr, traits=temp$traits)
	},

	#' @description
	#' Perform colocalisation systematically for all detected candidate traits for a given region
	#'
	#' @param rsid rsid i.e. one of the outcome GWAS hits
	#' @param radius Radius around which to perform analysis default=50000
	coloc_scan = function(rsid, radius=50000)
	{
		trait_ids <- unique(self$rsid_scan[[rsid]][["id"]])
		for(id in trait_ids)
		{
			message("coloc analysis of ", id)
			self$colocalise_rsid(rsid, id, radius)
		}
		self$coloc_result[[rsid]] <- lapply(names(self$coloc_data[[rsid]]), function(x){
			a <- self$coloc_data[[rsid]][[x]][["coloc"]][["summary"]] %>% 
				as.list() %>% 
				dplyr::as_tibble() %>%
				dplyr::mutate(id=x) %>%
				dplyr::select(id, dplyr::everything())
		}) %>% dplyr::bind_rows()
	},

	#' @description
	#' Perform MR analysis of all candidate traits for a particular rsid
	#'
	#' @param rsid GWAS hit from which candidate traits are obtained
	#' @param steiger_filtering Whether to perform steiger filtering default=TRUE
	#' @param exclude_rsid_region Whether to exclude the known GWAS hit default=TRUE
	#' @param radius Excluded region radius default=250000
	#' @param mrmethod MR method default="mr_ivw"
	mr = function(rsid, steiger_filtering=TRUE, exclude_rsid_region=TRUE, radius=250000, mrmethod="mr_ivw")
	{
		self$mr_scan[[rsid]] <- lapply(self$rsid_scan[[rsid]][["id"]], function(trait_id)
		{
			if(! trait_id %in% names(self$candidate_instruments))
			{
				self$candidate_instruments[[trait_id]] <- TwoSampleMR::extract_instruments(trait_id)
				rsid_needed <- self$candidate_instruments[[trait_id]][["SNP"]][! self$candidate_instruments[[trait_id]][["SNP"]] %in% self$candidate_instruments_outcome[["SNP"]]]
				if(length(rsid_needed) > 0)
				{
					self$candidate_instruments_outcome <- dplyr::bind_rows(self$candidate_instruments_outcome, TwoSampleMR::extract_outcome_data(rsid_needed, self$igd_id))
				}
			}
			dat <- TwoSampleMR::harmonise_data(self$candidate_instruments[[trait_id]], self$candidate_instruments_outcome)
			print(nrow(dat))
			if(steiger_filtering)
			{
				dat <- dat %>% 
					TwoSampleMR::steiger_filtering() %>%
					dplyr::filter(steiger_dir)
			}
			print(nrow(dat))
			if(exclude_rsid_region)
			{
				chr <- self$gwashits$chr[self$gwashits$rsid == rsid]
				position <- self$gwashits$position[self$gwashits$rsid == rsid]
				dat <- subset(dat, ! (chr == chr & pos < (position + radius) & pos > (position - radius)))
			}
			print(nrow(dat))
			if(nrow(dat) > 0)
			{
				return(TwoSampleMR::mr(dat, method_list=c("mr_wald_ratio", mrmethod)))
			} else {
				return(NULL)
			}
		}) %>% dplyr::bind_rows()
	},

	#' @description
	#' Perform multivariable MR for all or some of the candidate traits related to a particular region
	#'
	#' @param rsid RSID from which candidate traits are obtained
	#' @param traitlist List of traits to use. Default (NULL) is to use all, but important to scrutinise and manually exclude those which appear to be synonymous with the outcome
	mvmr = function(rsid, traitlist=NULL)
	{
		if(is.null(traitlist))
		{
			traitlist <- self$rsid_scan[[rsid]][["id"]]
			message("Using all ", length(traitlist), " traits in MVMR model")
		}
		self$mv_dat[[rsid]][["exposure"]] <- TwoSampleMR::mv_extract_exposures(traitlist)
		self$mv_dat[[rsid]][["outcome"]] <- TwoSampleMR::extract_outcome_data(self$mv_dat[[rsid]][["exposure"]][["SNP"]], self$igd_id)
		self$mv_dat[[rsid]][["dat"]] <- TwoSampleMR::mv_harmonise_data(self$mv_dat[[rsid]][["exposure"]], self$mv_dat[[rsid]][["outcome"]])
		self$mv_dat[[rsid]][["result"]] <- TwoSampleMR::mv_multiple(self$mv_dat[[rsid]][["dat"]])
	},

	#' @description
	#' This function will take the output from mvmr, perform feature selection from the analysed traits using MVMR-LASSO, and then re-estimate the MVMR associations with only the retained traits
	#'
	#' @param rsid Region from which MVMR analysis has been performed
	mvmr_lasso = function(rsid)
	{
		res <- TwoSampleMR::mv_lasso_feature_selection(self$mv_dat[[rsid]][["dat"]])
		self$mvmr(rsid, traitlist=res$exposure)

	}

))

