
splicingtest <- function() {
  require("testthat")
  tdir <- system.file("tests", package="splicing")
  do.call("test_dir", list(tdir))
}
