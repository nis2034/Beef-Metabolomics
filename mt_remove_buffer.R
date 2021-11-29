mt_remove_buffer <- function (D)
{
  # Removing 0 sample (buffer)
  
  cat("removing buffer samples\n")
  D = D[,-which(D$Material == "buffer")]
  
}
