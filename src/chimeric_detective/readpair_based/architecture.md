# Read-Pair Based Chimera Detection Architecture

## Core Components

### 1. Data Layer
- **BamParser**: Extract read-pair information from BAM files
- **ReadPairExtractor**: Process individual read pairs
- **InsertSizeCalculator**: Compute insert size distributions

### 2. Analysis Layer
- **BaselineEstimator**: Establish baseline metrics from proper pairs
- **SlidingWindowAnalyzer**: Window-based anomaly detection
- **StatisticalDetector**: Pluggable statistical methods
- **BreakpointScorer**: Calculate confidence scores

### 3. Configuration Layer
- **ConfigManager**: Handle all parameters and thresholds
- **ThresholdManager**: Dynamic threshold adjustment
- **WindowManager**: Flexible window sizing strategies

### 4. Output Layer
- **OutputFormatter**: Multiple format support (JSON, TSV, etc.)
- **VisualizationEngine**: Generate plots and reports
- **ResultAggregator**: Combine and rank predictions

### 5. Utilities
- **Logger**: Comprehensive logging with debug modes
- **Validator**: Input validation and error handling
- **MetricsCalculator**: Common statistical functions

## Key Design Principles

1. **Separation of Concerns**: Each component has a single responsibility
2. **Dependency Injection**: Components receive dependencies, not create them
3. **Configuration-Driven**: All parameters externally configurable
4. **Format Agnostic**: Support multiple input/output formats
5. **Extensible**: Easy to add new statistical methods or metrics
6. **Testable**: Each component independently testable
7. **Performance-Aware**: Efficient processing of large BAM files

## Data Flow

1. BAM file → BamParser → ReadPairExtractor
2. Read pairs → BaselineEstimator → Baseline metrics
3. Read pairs + Baseline → SlidingWindowAnalyzer → Anomaly signals
4. Anomaly signals → StatisticalDetector → Breakpoint candidates
5. Candidates → BreakpointScorer → Scored predictions
6. Predictions → OutputFormatter → Results (JSON/TSV/etc.)

## Configuration Schema

```yaml
# Example configuration structure
read_quality:
  min_mapping_quality: 20
  min_base_quality: 20
  require_proper_pairs: true

insert_size:
  outlier_method: "mad"  # or "iqr", "zscore", "kde"
  outlier_threshold: 3.0
  min_pairs_for_baseline: 1000

sliding_window:
  window_size: 1000
  step_size: 500
  min_pairs_per_window: 50

detection:
  proper_pair_drop_threshold: 0.5
  insert_size_shift_threshold: 2.0
  discordant_pair_threshold: 0.3

output:
  formats: ["json", "tsv"]
  include_debug_info: false
  visualization: true
```