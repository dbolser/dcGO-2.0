# dcGO Pipeline - TODO List

## 🚨 Critical Issues (Must Fix First)

### 1. Dependency Configuration
- [ ] Fix `sqlite3` dependency in pyproject.toml (sqlite3 is built-in to Python, should be removed)
- [ ] Remove conflicting dependencies causing uv sync failures
- [ ] Test clean `uv sync` installation
- [ ] Verify all scientific computing dependencies install correctly

### 2. InterProScan Strategy Decision
- [ ] Research InterProScan download vs local execution trade-offs
- [ ] Document recommended approach for domain annotation acquisition
- [ ] Update data acquisition module based on chosen strategy
- [ ] Create fallback options for different compute environments

### 3. Taxonomic Scope Decision
- [ ] Review papers in docs/ to understand taxonomic background requirements
- [ ] Determine if algorithm needs "all of biology" or can work with human-only
- [ ] Document taxonomic scope recommendations
- [ ] Update configuration and documentation accordingly

## 🧪 Testing Infrastructure (High Priority)

### 4. Unit Test Suite
- [ ] Create comprehensive unit tests for all 5 core modules:
  - [ ] `tests/unit/test_data_acquisition.py` - Download, validation, InterProScan setup
  - [ ] `tests/unit/test_domain_scanning.py` - Domain parsing, supra-domain generation
  - [ ] `tests/unit/test_statistical_inference.py` - Fisher's tests, FDR correction
  - [ ] `tests/unit/test_ontology_processor.py` - True Path Rule, propagation
  - [ ] `tests/unit/test_database_manager.py` - SQLite operations, exports
- [ ] Mock external dependencies (InterProScan, file downloads)
- [ ] Achieve >80% test coverage

### 5. Integration Tests
- [ ] `tests/integration/test_pipeline_stages.py` - Multi-module workflows
- [ ] `tests/integration/test_data_flow.py` - Data passing between modules
- [ ] Test checkpoint/resume functionality
- [ ] Test error recovery and validation

### 6. End-to-End Tests
- [ ] `tests/e2e/test_small_dataset.py` - Complete pipeline with toy data
- [ ] Create sample datasets for testing (100-1000 proteins)
- [ ] Test CLI interface and command-line arguments
- [ ] Validate output database structure and content

## 🔧 Technical Improvements (Medium Priority)

### 7. Configuration Enhancements
- [ ] Add environment-specific configs (dev/test/prod)
- [ ] Implement config validation with pydantic
- [ ] Add CLI config override options
- [ ] Document all configuration parameters

### 8. Error Handling & Logging
- [ ] Standardize error handling across all modules
- [ ] Improve progress reporting for long-running operations
- [ ] Add structured logging with correlation IDs
- [ ] Create error recovery documentation

### 9. Performance Optimization
- [ ] Profile memory usage during large dataset processing
- [ ] Optimize correspondence matrix construction
- [ ] Add parallel processing for statistical tests
- [ ] Implement streaming for large file operations

## 📚 Documentation (Medium Priority)

### 10. User Documentation
- [ ] Create step-by-step installation guide
- [ ] Write beginner tutorial with example data
- [ ] Document hardware requirements and scaling guidance
- [ ] Create troubleshooting guide for common issues

### 11. Developer Documentation
- [ ] API reference documentation with sphinx
- [ ] Architecture decision records (ADRs)
- [ ] Contributing guidelines
- [ ] Code style and review guidelines

## 🚀 Production Readiness (Lower Priority)

### 12. Deployment & Distribution
- [ ] Create Docker images for different environments
- [ ] Set up CI/CD pipeline with GitHub Actions
- [ ] Create conda package for easy installation
- [ ] Add version management and release process

### 13. Advanced Features
- [ ] Web interface for pipeline results visualization
- [ ] REST API for programmatic access
- [ ] Integration with Galaxy workflow system
- [ ] Support for custom ontologies beyond GO

## 📊 Data & Validation

### 14. Dataset Management
- [ ] Create curated test datasets
- [ ] Implement data versioning for reproducibility
- [ ] Add dataset integrity validation
- [ ] Create benchmark datasets for performance testing

### 15. Scientific Validation
- [ ] Compare results with original dcGO implementation
- [ ] Validate statistical methodology implementation
- [ ] Cross-validate with other function prediction tools
- [ ] Document known limitations and assumptions

## 🎯 Immediate Next Steps (Today's Work)

1. **Fix pyproject.toml** - Remove sqlite3, test uv sync
2. **Review docs PDFs** - Understand taxonomic requirements
3. **Create basic test suite** - Get testing infrastructure working
4. **Test InterProScan strategy** - Research download options
5. **Validate end-to-end flow** - Run with minimal test data

## 📝 Progress Tracking

### Completed ✅
- [x] Core pipeline architecture (5 modules, 5000+ lines)
- [x] Configuration system
- [x] Database schema and management
- [x] Statistical inference engine
- [x] Ontology processing with True Path Rule
- [x] Comprehensive documentation (README, IMPLEMENTATION_GUIDE, CLAUDE.md)
- [x] Project structure and build configuration

### In Progress 🔄
- [ ] Dependency resolution and environment setup
- [ ] Test suite development
- [ ] InterProScan integration strategy
- [ ] Taxonomic scope determination

### Blocked ⚠️
- Testing requires dependency resolution
- Full pipeline testing requires InterProScan decision
- Performance validation requires taxonomic scope decision

---

**Priority Order for Completion:**
1. Fix dependencies → 2. Basic tests → 3. Strategic decisions → 4. Full validation → 5. Documentation → 6. Production features