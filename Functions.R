`%notin%` <- function(x, table) {
  !(x %in% table)
}
to_excel_numeric <- function(x) {
  if (is.null(x)) return(NA)
  
  # If it's already a Date or POSIXct (real date), convert to Excel Serial
  if (inherits(x, c("Date", "POSIXt"))) {
    return(as.numeric(as.Date(x)) + 25569)
  }
  
  # If it's a character string (e.g., "2020-02-01"), convert to Date then Serial
  if (is.character(x)) {
    # Check if it's a numeric string first (like "46000")
    if (all(grepl("^[0-9.]+$", na.omit(x)))) {
      return(as.numeric(x))
    }
    parsed_date <- as.Date(x)
    return(as.numeric(parsed_date) + 25569)
  }
  
  # If it's already numeric, just return it
  return(as.numeric(x))
}