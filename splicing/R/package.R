
.onAttach <- function(library, pkg) {
  invisible()
}

.onLoad <- function(dir, package) {
  library.dynam("splicing", package, dir, local=FALSE)
}

.onUnload <- function(libpath) {
  library.dynam.unload("splicing", libpath)
}

.Last.lib <- function(libpath) {
  .onUnload(libpath)
}

splicing_version <- function() {
  .Call("R_splicing_version", PACKAGE="splicing")
}
