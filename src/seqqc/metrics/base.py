"""Base class and registry for metrics (the plugin interface).

Every metric is a small class with three methods:

* update(batch) - receive one batch of reads (called many times)
* result()      - return a JSON-serialisable dict (called once, at the end)
* plot(ax)      - optional: draw the metric on a matplotlib Axes

The core reads the file once and hands every batch to all registered metrics.
Adding a new metric = one new class + the @register_metric decorator; the core
does not change.
"""

from abc import ABC, abstractmethod

_REGISTRY: dict[str, type["BaseMetric"]] = {}


class BaseMetric(ABC):
    #: unique key used in the JSON report and on the command line
    name: str = ""

    @abstractmethod
    def update(self, batch) -> None:
        """Consume one batch of reads (the batch type is defined in week 2)."""

    @abstractmethod
    def result(self) -> dict:
        """Return the final values as a JSON-serialisable dict."""

    def plot(self, ax) -> None:
        """Optional: draw the metric on ``ax``. The default draws nothing."""
        return None


def register_metric(cls):
    """Class decorator that adds a metric to the registry."""
    if not cls.name:
        raise ValueError(f"{cls.__name__} must define a non-empty 'name'")
    if cls.name in _REGISTRY:
        raise ValueError(f"metric name already registered: {cls.name}")
    _REGISTRY[cls.name] = cls
    return cls


def get_registered_metrics() -> dict[str, type[BaseMetric]]:
    """Return a copy of the registry (name -> class)."""
    return dict(_REGISTRY)
