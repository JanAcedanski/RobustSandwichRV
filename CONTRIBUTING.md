# Contributing

Please open an issue before proposing a change to covariance
conventions, critical values, or certification status semantics. Pull
requests should:

1.  include focused `testthat` coverage;
2.  preserve the frozen Julia reference fixtures;
3.  pass `R CMD check` on Windows, macOS, and Linux;
4.  distinguish verified witnesses from local-search candidates; and
5.  never interpret numerical search failure as global infeasibility.
