.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion("MyoScore")
  packageStartupMessage(
    sprintf("MyoScore v%s", ver),
    "\n  A genetically-informed transcriptomic scoring system",
    "\n  for quantifying human skeletal muscle health.",
    "\n  5 dimensions: Strength, Mass, LeanMuscle, Youth, Resilience",
    "\n  Use myoscore_score() to calculate scores from raw counts."
  )
}
