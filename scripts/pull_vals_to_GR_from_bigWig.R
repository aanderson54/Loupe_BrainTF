#functions needed to pull signal values from bigwigs and take median within the union of peaks for a tissue/cell type
#Author: Brian Roberts



require(rtracklayer)

read_and_bin_func_bigWig <- function(bin_GR,bigWig_path,bin_func = function(x){quantile(x,0.5)}, run_unstable=F){
	#the which input to rtracklayer::import cannnot be named; save names and put back later
	names_bin_GR <- names(bin_GR)
	bin_GR <- unname(bin_GR)
	big_wig_data_GR <- import(con= bigWig_path,which= bin_GR)
	if(run_unstable){
		binned_data <- trackGR_to_func_bins_FAST_UNSTABLE(big_wig_data_GR, bin_GR, bin_func= bin_func)
	}
	else{
		binned_data <- trackGR_to_func_bins(big_wig_data_GR, bin_GR, bin_func= bin_func)
	}
	
	bigWig_data_func_to_bin_GR <- granges(bin_GR)
	bigWig_data_func_to_bin_GR$bin_vals <- binned_data
	names(bigWig_data_func_to_bin_GR) <- names_bin_GR
	return(bigWig_data_func_to_bin_GR)
}


#this function is called by read_and_bin_func_bigWig
trackGR_to_func_bins = function(track_GR,bin_GR,bin_func=function(x){quantile(x,0.5)},minoverlap=1,value_name=determine_meta_data_col(mcols(track_GR)),ignore_strand=T,verbose=F){
	#if length(min_overlap) = 1, i.e its not a vector run findOverlaps for speed
	if(length(minoverlap)==1){
		print("minoverlap is a scalar. Running fast overlaps.")
		
		if(minoverlap<0){
			hits_obj <- findOverlaps(bin_GR,track_GR, maxgap= -minoverlap,ignore.strand= ignore_strand)
		}
		else{
			hits_obj <- findOverlaps(bin_GR,track_GR, minoverlap= minoverlap,ignore.strand= ignore_strand)
		}
		
		if(verbose){print("overlaps complete")}
		track_to_bin_hits_with_size_keep_df <- as.data.frame(hits_obj)
		#print("coerced hits to data.frame")
		colnames(track_to_bin_hits_with_size_keep_df) = c("queryHit","subjectHit")
	}
	else{
		track_to_bin_hits_with_size_df = get_overlaps_size(track_GR, bin_GR)
		if(verbose){print("overlaps complete")}
		#print("coerced hits to data.frame")
		mask_keep_hits = which(track_to_bin_hits_with_size_df$overlap_bp >= minoverlap)
		track_to_bin_hits_with_size_keep_df = track_to_bin_hits_with_size_df[mask_keep_hits,]
	}
	out_list_vals <- list()
	out_list_vals[1:length(bin_GR)] = NA
	#out_vals = rep(0,length(bin_GR))
	
	all_val_df = mcols(track_GR)
	
	mask_data_col = which(colnames(all_val_df)== value_name)
	
	all_vals = all_val_df[, mask_data_col]
	
	matched_vals = all_vals[track_to_bin_hits_with_size_keep_df$subjectHit]
	
	quant_vals_list = tapply(matched_vals,INDEX= track_to_bin_hits_with_size_keep_df$queryHit,FUN= bin_func,simplify=F)
	
	uniq_bin_GR_hits <- unique(track_to_bin_hits_with_size_keep_df$queryHit)
	
	out_list_vals[uniq_bin_GR_hits] <- quant_vals_list
	
	#check to see if out_list_vals can be coerced to a numeric vector
	unlist_same_size_check <- length(out_list_vals) == length(unlist(out_list_vals))
	
	all_can_be_numeric_check <- all_can_be_coerced_to_numeric(unlist(out_list_vals))
	
	if(all(c(unlist_same_size_check,all_can_be_numeric_check))){
		out_list_vals = unlist(out_list_vals)
		out_list_vals[which(is.na(out_list_vals))] <-0

	}
	
	return(out_list_vals)
}

determine_meta_data_col = function(extra_data_df,data_set_name=NULL,verbose=F){
	extra_data_colnames = colnames(extra_data_df)
	if(ncol(extra_data_df)==1){
		binned_data_name = extra_data_colnames[1]
		if(verbose){
			print(paste("one metadata column found for", data_set_name,"; using", binned_data_name))
		}
	}
	else{
		mask_found_cnts = which(extra_data_colnames == "cnts")
		if(length(mask_found_cnts)==1){
			binned_data_name = extra_data_colnames[mask_found_cnts]
			if(verbose){
				print(paste("multiple metadata columns found for", data_set_name,"; using", binned_data_name))
			}
		}
		else{
			binned_data_name = extra_data_colnames[1]
			if(verbose){
				print(paste("multiple metadata columns found for", data_set_name,"; using", binned_data_name))
			}
		}
	}
	return(binned_data_name)
}

all_can_be_coerced_to_numeric = function(test_vector){
	no_na_vec = test_vector[which(!is.na(test_vector))]
	op <- options(warn=2)
	try_coerce_to_numeric = try(as.numeric(no_na_vec),silent=T)
	options(op)
	if(is(try_coerce_to_numeric,"try-error")){return(F)}else{return(T)}
}

