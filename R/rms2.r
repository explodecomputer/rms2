rms2 <- R6::R6Class("rms2", list(
	igd_id = NULL,
	info = NULL,
	gwashits = NULL,
	rsid_scan = list(),
	coloc_data = list(),
	coloc_result = list(),
	candidate_instruments = NULL,

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

	extract_gwashits = function(igd_id=self$igd_id, pval=5e-8, kb=10000, r2=0.001, pop='EUR')
	{
		self$gwashits <- ieugwasr::tophits(igd_id, pval=pval, r2=r2, kb=kb, pop=pop)
		message("Found ", nrow(self$gwashits), " GWAS hits for ", igd_id)
	},

	scan_rsid = function(rsid, pval=5e-8, idlist=NULL)
	{
		if(is.null(idlist))
		{
			res <- ieugwasr::phewas(rsid, pval=pval)
		} else {
			res <- ieugwasr::associations(rsid, idlist) %>%
				dplyr::filter(p < pval)
		}
		for(r in rsid)
		{
			message("Found ", nrow(res), " associations for ", r)
			self$rsid_scan[[r]] <- subset(res, rsid == r)
		}
	},

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

	plot_coloc = function(rsid, trait_id)
	{
		temp <- gwasglue::coloc_to_gassocplot(self$coloc_data[[rsid]][[trait_id]][["data"]])
		require(gassocplot)
		gassocplot::stack_assoc_plot(temp$markers, temp$z, temp$corr, traits=temp$traits)
	},

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

	extract_candidate_instruments = function()
	{
		# something
	}

))