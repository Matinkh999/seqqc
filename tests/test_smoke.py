import pytest

import seqqc
from seqqc.metrics import base
from seqqc.metrics.base import BaseMetric, get_registered_metrics, register_metric


def test_version():
    assert seqqc.__version__ == "0.1.0"


def test_register_metric():
    @register_metric
    class Dummy(BaseMetric):
        name = "dummy_for_test"

        def update(self, batch):
            pass

        def result(self):
            return {}

    try:
        assert "dummy_for_test" in get_registered_metrics()
    finally:
        base._REGISTRY.pop("dummy_for_test", None)


def test_register_requires_a_name():
    class NoName(BaseMetric):
        def update(self, batch):
            pass

        def result(self):
            return {}

    with pytest.raises(ValueError):
        register_metric(NoName)
