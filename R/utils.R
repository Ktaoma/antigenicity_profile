normalize_date <- function(x) {
  # year only -> yyyy-01-01
  x <- ifelse(grepl("^\\d{4}$", x), paste0(x, "-01-01"), x)
  
  # year-month only -> yyyy-mm-01
  x <- ifelse(grepl("^\\d{4}-\\d{2}$", x), paste0(x, "-01"), x)
  
  x_edit <- as.POSIXct(x, format="%Y-%m-%d", tz="UTC")
  
  return(x_edit)
}



translate_AA <- function(X){
  
  ss <- list(translate(subseq(X,1),if.fuzzy.codon = "X"),
             translate(subseq(X,2),if.fuzzy.codon = "X"),
             translate(subseq(X,3),if.fuzzy.codon = "X"),
             translate(subseq(reverseComplement(X),1),if.fuzzy.codon = "X"),
             translate(subseq(reverseComplement(X),2),if.fuzzy.codon = "X"),
             translate(subseq(reverseComplement(X),3),if.fuzzy.codon = "X"))
  
  min_orf <- 0
  final_frame <- list()
  for (t in ss) {
    
    aa_split <- unlist(strsplit(as.character(t), "\\*"))
    aa_longest_idx <- aa_split |> nchar() |> which.max()
    aa_longest <- aa_split[aa_longest_idx] |> AAStringSet()
    print(nchar(aa_longest))
    if (nchar(aa_longest) > min_orf) {
      
      names(aa_longest) <- "sample"
      final_frame[[1]] <- aa_longest
      
      min_orf <- nchar(aa_longest)
    }
  }
  return(final_frame[[1]])
}
