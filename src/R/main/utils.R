
coalesce <- function(v, default){
  #' Returns the given value if not null; returns #default value otherwise
  ifelse(is.null(v), default, v)
}

lib <- function(p) {
  #' Requires the given package and installs it first if necessary
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}