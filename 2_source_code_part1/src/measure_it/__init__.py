"""measure-it-public: public-data alpha of the "Measure It to Cure It" measurement-validation and deployment engine.

With MEASURE_IT_OFFLINE=1 in the environment (set by `measure-it pipeline --offline`) every process that
imports this package, including joblib/multiprocessing workers, installs the socket-level offline guard
(measure_it.http.install_offline_guard), so no step can reach the network.
"""
import os as _os

if _os.environ.get("MEASURE_IT_OFFLINE", "").strip().lower() in {"1", "true", "yes", "on"}:
    from .http import install_offline_guard as _install_offline_guard

    _install_offline_guard()
