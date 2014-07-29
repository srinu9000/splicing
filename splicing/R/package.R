
.onLoad <- function(lib, pkg) {
  library.dynam("splicing", pkg, lib)
}

splicing_version <- function() {
  .Call("R_splicing_version", PACKAGE="splicing")
}
