import asyncio
import asyncio
import sys
import types

import pytest

if "psutil" not in sys.modules:
    def _cpu_count():
        return 4

    def _virtual_memory():
        return types.SimpleNamespace(total=8 * 1024**3)

    def _disk_usage(_path):
        return types.SimpleNamespace(free=8 * 1024**3)

    sys.modules["psutil"] = types.SimpleNamespace(
        cpu_count=_cpu_count,
        virtual_memory=_virtual_memory,
        disk_usage=_disk_usage,
    )

if "aiohttp" not in sys.modules:
    class _DummySession:  # pragma: no cover - stub for import
        def __init__(self, *args, **kwargs):
            raise RuntimeError("aiohttp is not available in tests")

    sys.modules["aiohttp"] = types.SimpleNamespace(
        ClientSession=_DummySession,
        ClientTimeout=lambda *args, **kwargs: None,
        TCPConnector=lambda *args, **kwargs: None,
    )

if "aiofiles" not in sys.modules:
    async def _raise(*args, **kwargs):  # pragma: no cover - stub for import
        raise RuntimeError("aiofiles is not available in tests")

    sys.modules["aiofiles"] = types.SimpleNamespace(open=_raise)

if "requests" not in sys.modules:
    class _DummyResponse:  # pragma: no cover - stub for import
        status_code = 200

        def raise_for_status(self):
            return None

        def iter_content(self, chunk_size=8192):
            return iter([])

    sys.modules["requests"] = types.SimpleNamespace(get=lambda *args, **kwargs: _DummyResponse())

if "Bio" not in sys.modules:
    seqio_stub = types.SimpleNamespace(parse=lambda *args, **kwargs: iter([]))
    sys.modules["Bio"] = types.SimpleNamespace(SeqIO=seqio_stub)
    sys.modules["Bio.SeqIO"] = seqio_stub

if "networkx" not in sys.modules:
    class _DummyGraph:
        def __init__(self):
            self._nodes = []

        def nodes(self, data=False):
            return []

        def edges(self):
            return []

        def in_degree(self, _node):
            return 0

        def predecessors(self, _node):
            return []

        def remove_nodes_from(self, _nodes):
            return None

    sys.modules["networkx"] = types.SimpleNamespace(
        DiGraph=_DummyGraph,
        is_directed_acyclic_graph=lambda *args, **kwargs: True,
        ancestors=lambda *args, **kwargs: set(),
        descendants=lambda *args, **kwargs: set(),
        single_source_shortest_path_length=lambda *args, **kwargs: {},
    )

if "obonet" not in sys.modules:
    import types as _types

    def _read_obo(*args, **kwargs):
        return sys.modules["networkx"].DiGraph()

    sys.modules["obonet"] = _types.SimpleNamespace(read_obo=_read_obo)

if "pandas" not in sys.modules:
    sys.modules["pandas"] = types.SimpleNamespace(DataFrame=lambda *args, **kwargs: None)

if "scipy" not in sys.modules:
    stats_stub = types.SimpleNamespace(
        fisher_exact=lambda *args, **kwargs: (1.0, 1.0),
        hypergeom=types.SimpleNamespace(sf=lambda *args, **kwargs: 1.0),
    )
    sys.modules["scipy"] = types.SimpleNamespace(stats=stats_stub)
    sys.modules["scipy.stats"] = stats_stub

if "statsmodels" not in sys.modules:
    def _fdrcorrection(p_values, alpha=0.05, method="indep"):
        rejected = [p <= alpha for p in p_values]
        return rejected, [min(alpha, max(0.0, p)) for p in p_values]

    multitest_stub = types.SimpleNamespace(fdrcorrection=_fdrcorrection)
    sys.modules["statsmodels"] = types.SimpleNamespace(stats=types.SimpleNamespace(multitest=multitest_stub))
    sys.modules["statsmodels.stats"] = types.SimpleNamespace(multitest=multitest_stub)
    sys.modules["statsmodels.stats.multitest"] = multitest_stub

from src.main_pipeline import dcGOPipeline
from src.ontology_processor import Annotation
from src.statistical_inference import AssociationResult


class DummyOntologyProcessor:
    def __init__(self):
        self.filter_calls = []
        self.propagate_calls = []

    def apply_optimal_level_filter(
        self,
        significant_associations,
        protein_domain_map,
        protein_go_map,
        min_background_size,
        alpha_threshold,
    ):
        self.filter_calls.append(
            (significant_associations, min_background_size, alpha_threshold)
        )
        # Return the first association to simulate filtering
        return list(significant_associations[:1])

    def propagate_annotations(self, direct_associations):
        self.propagate_calls.append(direct_associations)
        if not direct_associations:
            return []

        source = direct_associations[0]
        score = max(1.0, min(100.0, getattr(source, "hyper_score", 1.0)))
        q_value = getattr(source, "q_value", 0.01) or 0.01

        return [
            Annotation(
                domain=source.domain,
                go_term=source.go_term,
                q_value=q_value,
                association_score=score,
                annotation_type="direct",
                direct_source_term=source.go_term,
            ),
            Annotation(
                domain=source.domain,
                go_term="GO:ANCESTOR",
                q_value=q_value,
                association_score=score,
                annotation_type="propagated",
                direct_source_term=source.go_term,
            ),
        ]


def test_apply_true_path_rule_disabled():
    pipeline = dcGOPipeline(enable_async=False)
    pipeline.ontology_processor = DummyOntologyProcessor()

    result = asyncio.run(
        pipeline._apply_true_path_rule(
            significant_associations=[],
            protein_domain_map={"P1": ["PF00001"]},
            protein_go_map={"P1": {"GO:0001"}},
            parameters={"apply_true_path": False},
        )
    )

    assert result["enabled"] is False
    assert result["filtered_associations"] == []
    assert result["propagated_annotations"] == []

def test_apply_true_path_rule_enabled_propagates():
    pipeline = dcGOPipeline(enable_async=False)
    stub_processor = DummyOntologyProcessor()
    pipeline.ontology_processor = stub_processor

    associations = [
        AssociationResult(
            domain="PF00001",
            go_term="GO:0001",
            a=5,
            b=3,
            c=2,
            d=10,
            p_value=1e-5,
            odds_ratio=5.0,
            hyper_score=40.0,
            q_value=1e-3,
        ),
        AssociationResult(
            domain="PF00002",
            go_term="GO:0002",
            a=4,
            b=2,
            c=1,
            d=12,
            p_value=2e-4,
            odds_ratio=3.5,
            hyper_score=30.0,
            q_value=5e-3,
        ),
    ]

    result = asyncio.run(
        pipeline._apply_true_path_rule(
            significant_associations=associations,
            protein_domain_map={"P1": ["PF00001", "PF00002"]},
            protein_go_map={"P1": {"GO:0001", "GO:0002"}},
            parameters={
                "apply_true_path": True,
                "true_path_min_background": 5,
                "true_path_alpha_threshold": 0.05,
            },
        )
    )

    assert result["enabled"] is True
    assert len(result["filtered_associations"]) == 1
    assert len(result["propagated_annotations"]) == 2
    assert result["statistics"]["direct_annotations"] == 1
    assert result["statistics"]["propagated_annotations"] == 1

    assert stub_processor.filter_calls, "Filter should be invoked"
    call = stub_processor.filter_calls[0]
    assert call[1] == 5
    assert call[2] == 0.05

    assert stub_processor.propagate_calls, "Propagation should be invoked"
    assert len(stub_processor.propagate_calls[0]) == 1
