# Time-Series Primitives

Status: initial internal univariate time-series primitives added; validation is ongoing.

The layer provides full-sample-mean autocovariance with biased (n) and unbiased (n−lag) denominators, biased autocorrelation, a sampled ACF through a caller-selected maximum lag, and one-step lag differencing. Seasonal differencing is obtained by selecting the seasonal period as the lag; higher orders can be composed.

All input observations must be finite. ACF uses the biased autocovariance sequence and requires positive sample variance. Differencing validates all outputs before writing so invalid input or overflow leaves the caller buffer unchanged. ACF uses a temporary output buffer and commits only after all requested lags succeed.

Native tests cover hand-computable autocovariances/ACF, biased versus unbiased denominators, lag differencing, constant-series behavior, and output preservation.
