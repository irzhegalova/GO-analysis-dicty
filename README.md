# GO-analysis-dicty

## example file for GO-dicty usage

	source('scripts/functions.R')
	res.for.GO <- read.table('data/res.forGO.example.txt')
res.forGO has three columns: id - GO-term ID
			#     cl - number of group
                        #     GENE = gene name

	GO_list <- lapply(c('BP', 'MF', 'CC'),
		function(term) {
		name <- paste0('data/GO_filtered.', term, '.tsv')
			    # read simplified GO-db
		simp <- read.table(name, header = FALSE)
		modeF <- 'DE'
		go <- GO.analysis.ForGeneLists(res.for.GO, mode = modeF, ONT = term, GOfilter = simp)
		# if needed you can write results to files
		# name <- paste0('results/GO_k4_', modeF, '_forcluster_more10',cluster, '.', term, '.tsv')
		# write.table(go, name, quote = FALSE, row.names = FALSE,col.names = FALSE, sep='\t')
		return(go)
	})
### if needed, you can simplify it
	GO_df <- do.call(rbind.data.frame, GO_list)
