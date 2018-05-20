#' @title Preparing ordered counts data for ordered topic model.
#' @description Ordering counts data based on an ordering vector and padding
#' zeros at appropriate locations, so that the total number of columns is a
#' power of 2, a constraint one needs to apply the ordered topic model.
#'
#' @param local_data The abundance counts data for site/species.
#' @param ordering_vec The vector used for ordering the columns of the data
#'                     matrix (species in this scenario).
#' @param cut_off_percentile The percentile cut off used for padding the zeros.
#'
#' @return Returns an ordered counts matrix with the species along the columns
#' ordered by the \item{ordering_vec}.
#'

order_counts <- function(local_data,
                        ordering_vec,
                        threshmax = 5){
if(length(ordering_vec) != dim(local_data)[2]){
  stop("The length of the ordering vector must matches the number of species");
}
 ordered_vec <- ordering_vec[order(ordering_vec, decreasing=FALSE)];
 total_zeros_padded <- 2^{ceiling(log(dim(local_data)[2], base=2))} - dim(local_data)[2];

 ranges <- array(0, length(ordered_vec)-1);
 for(num in 1:(length(ordered_vec)-1)){
   ranges[num] <- ordered_vec[num+1] - ordered_vec[num];
 }
 prop_obs <- ranges/sum(ranges)
 num_zeros_padded_vec <- ceiling(prop_obs*total_zeros_padded);
 num_zeros_padded_vec[num_zeros_padded_vec < 5] = 0;
 num_intervals <- length(which(num_zeros_padded_vec > 0));
 indices_to_be_filled <- order(ranges, decreasing=TRUE)[1:num_intervals];

 new_ordering_vec <- as.numeric();
 for(num in 1:length(indices_to_be_filled)){
   out <- seq(ordered_vec[indices_to_be_filled[num]], ordered_vec[indices_to_be_filled[num]+1], length.out=num_zeros_padded_vec[indices_to_be_filled[num]]+2)
   out <- out[-c(1, length(out))];
   new_ordering_vec <- c(new_ordering_vec, out);
 }

 diff <- (length(new_ordering_vec) +  length(ordered_vec)) - 2^{ceiling(log(dim(local_data)[2], base=2))}
 if(diff  > 0){
   new_ordering_vec <- new_ordering_vec[-(1:diff)];
 }else if(diff < 0){
   k=1;
   while(diff < 0){
     new_ordering_vec <- c(new_ordering_vec, (max(ordered_vec) + k*min(ranges[ranges >0])));
     k=k+1;
     diff = diff + 1;
   }
 }

 padded_zero <- matrix(0, dim(local_data)[1], length(new_ordering_vec));
 pooled_ordering_vec <- c(ordering_vec, new_ordering_vec);
 pooled_data <- cbind(local_data, padded_zero);
 pooled_ordered_data <- pooled_data[,order(pooled_ordering_vec, decreasing=FALSE)];
 pooled_ordered_vec <- pooled_ordering_vec[order(pooled_ordering_vec, decreasing=FALSE)];
 ll <- list("data"=pooled_ordered_data, "ordered_vec"=pooled_ordered_vec)
 return(ll);
}
