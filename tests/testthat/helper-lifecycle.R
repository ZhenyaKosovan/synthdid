# Suppress lifecycle deprecation warnings globally in tests.
# Snapshot tests for deprecation messages should use
# withr::local_options(lifecycle_verbosity = "warning") to override.
options(lifecycle_verbosity = "quiet")
