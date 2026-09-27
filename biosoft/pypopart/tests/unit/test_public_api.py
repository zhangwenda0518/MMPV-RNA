"""Smoke tests for the top-level public API."""

import pypopart


class TestPublicAPI:
    """Every advertised name imports and the one-call workflow works."""

    def test_all_names_importable(self):
        """Each name in __all__ resolves to a real attribute."""
        for name in pypopart.__all__:
            assert getattr(pypopart, name) is not None

    def test_version_is_string(self):
        """The package exposes a version string."""
        assert isinstance(pypopart.__version__, str)

    def test_one_call_workflow(self):
        """build_network runs end to end from an Alignment."""
        alignment = pypopart.Alignment(
            [
                pypopart.Sequence('s1', 'AAAA'),
                pypopart.Sequence('s2', 'AAAT'),
                pypopart.Sequence('s3', 'AATT'),
            ]
        )
        network = pypopart.build_network('mst', alignment)
        assert isinstance(network, pypopart.HaplotypeNetwork)
        assert network.is_connected()
        assert network.num_nodes == 3
