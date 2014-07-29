
splicingtest <- function() {
  do.call(library, list("testthat"))
  tdir <- system.file("tests", package="splicing")
  do.call("test_dir", list(tdir))
}
