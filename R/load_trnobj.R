loadModule("trn_module", TRUE)
RcppParallel::setThreadOptions(numThreads = RcppParallel::defaultNumThreads() - 1)
