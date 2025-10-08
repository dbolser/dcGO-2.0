import asyncio
import importlib
import sys
import types

import pytest


@pytest.fixture
def pipeline_components(monkeypatch):
    """Provide core pipeline classes with external dependencies safely stubbed."""

    def ensure_module(name: str, stub: types.SimpleNamespace) -> None:
        if name not in sys.modules:
            monkeypatch.setitem(sys.modules, name, stub)

    ensure_module(
        "psutil",
        types.SimpleNamespace(
            cpu_count=lambda: 4,
            virtual_memory=lambda: types.SimpleNamespace(total=8 * 1024**3),
            disk_usage=lambda _path: types.SimpleNamespace(free=8 * 1024**3),
        ),
    )

    ensure_module(
        "aiohttp",
        types.SimpleNamespace(
            ClientSession=lambda *args, **kwargs: (_ for _ in ()).throw(
                RuntimeError("aiohttp is not available in tests")
            ),
            ClientTimeout=lambda *args, **kwargs: None,
            TCPConnector=lambda *args, **kwargs: None,
        ),
    )

    ensure_module(
        "aiofiles",
        types.SimpleNamespace(
            open=lambda *args, **kwargs: (_ for _ in ()).throw(
                RuntimeError("aiofiles is not available in tests")
            )
        ),
    )

    ensure_module(
        "requests",
        types.SimpleNamespace(
            get=lambda *args, **kwargs: types.SimpleNamespace(
                status_code=200,
                raise_for_status=lambda: None,
                iter_content=lambda chunk_size=8192: iter([]),
            )
        ),
    )

    if "Bio" not in sys.modules:
        seqio_stub = types.SimpleNamespace(parse=lambda *args, **kwargs: iter([]))
        ensure_module("Bio", types.SimpleNamespace(SeqIO=seqio_stub))
        ensure_module("Bio.SeqIO", seqio_stub)

    class _DummyGraph:
        def __init__(self):
            self._nodes = []

        def nodes(self, data=False):  # pragma: no cover - simple stub
            return []

        def edges(self):  # pragma: no cover - simple stub
            return []

        def in_degree(self, _node):  # pragma: no cover - simple stub
            return 0

        def predecessors(self, _node):  # pragma: no cover - simple stub
            return []

        def remove_nodes_from(self, _nodes):  # pragma: no cover - simple stub
            return None

    ensure_module(
        "networkx",
        types.SimpleNamespace(
            DiGraph=_DummyGraph,
            is_directed_acyclic_graph=lambda *args, **kwargs: True,
            ancestors=lambda *args, **kwargs: set(),
            descendants=lambda *args, **kwargs: set(),
            single_source_shortest_path_length=lambda *args, **kwargs: {},
        ),
    )

    ensure_module(
        "obonet",
        types.SimpleNamespace(
            read_obo=lambda *args, **kwargs: sys.modules["networkx"].DiGraph()
        ),
    )

    ensure_module(
        "pandas",
        types.SimpleNamespace(DataFrame=lambda *args, **kwargs: None),
    )

    stats_stub = types.SimpleNamespace(
        fisher_exact=lambda *args, **kwargs: (1.0, 1.0),
        hypergeom=types.SimpleNamespace(sf=lambda *args, **kwargs: 1.0),
    )
    ensure_module("scipy", types.SimpleNamespace(stats=stats_stub))
    ensure_module("scipy.stats", stats_stub)

    multitest_stub = types.SimpleNamespace(
        fdrcorrection=lambda p_values, alpha=0.05, method="indep": (
            [p <= alpha for p in p_values],
            [min(alpha, max(0.0, p)) for p in p_values],
        )
    )
    ensure_module(
        "statsmodels",
        types.SimpleNamespace(stats=types.SimpleNamespace(multitest=multitest_stub)),
    )
    ensure_module("statsmodels.stats", types.SimpleNamespace(multitest=multitest_stub))
    ensure_module("statsmodels.stats.multitest", multitest_stub)

    main_pipeline = importlib.import_module("src.main_pipeline")
    ontology_module = importlib.import_module("src.ontology_processor")
    stats_module = importlib.import_module("src.statistical_inference")

    return types.SimpleNamespace(
        dcGOPipeline=main_pipeline.dcGOPipeline,
        Annotation=ontology_module.Annotation,
        AssociationResult=stats_module.AssociationResult,
    )


@pytest.fixture
def dummy_ontology_processor(pipeline_components):
    Annotation = pipeline_components.Annotation

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
            return list(significant_associations[:1])

        def propagate_annotations(self, direct_associations):
            self.propagate_calls.append(direct_associations)
            if not direct_associations:
                return []

            source = direct_associations[0]
            score = max(1.0, min(100.0, getattr(source, "hyper_score", 1.0)))
            q_value = getattr(source, "q_value", None)
            if q_value is None:
                q_value = 0.01

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

    return DummyOntologyProcessor()


def test_apply_true_path_rule_disabled(pipeline_components, dummy_ontology_processor):
    pipeline = pipeline_components.dcGOPipeline(enable_async=False)
    pipeline.ontology_processor = dummy_ontology_processor

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


def test_apply_true_path_rule_enabled_propagates(
    pipeline_components, dummy_ontology_processor
):
    pipeline = pipeline_components.dcGOPipeline(enable_async=False)
    pipeline.ontology_processor = dummy_ontology_processor
    AssociationResult = pipeline_components.AssociationResult

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

    assert dummy_ontology_processor.filter_calls, "Filter should be invoked"
    call = dummy_ontology_processor.filter_calls[0]
    assert call[1] == 5
    assert call[2] == 0.05

    assert dummy_ontology_processor.propagate_calls, "Propagation should be invoked"
    assert len(dummy_ontology_processor.propagate_calls[0]) == 1
