"""
Command-line interface for Chimeric Detective.
"""

import os
import sys
import tempfile
import logging
from pathlib import Path
from typing import Optional, List, Tuple, Dict
import click
from tqdm import tqdm

from . import __version__
from .detector import ChimeraDetector
from .analyzer import ChimeraAnalyzer
from .resolver import ChimeraResolver
from .visualizer import ChimeraVisualizer
from .utils import (
    validate_file_exists, create_output_directory, check_external_tools,
    parse_reads_pattern, setup_logging
)
from .multi_sample import MultiSampleProcessor
from .config import ConfigManager


@click.command()
@click.option('--config', '-c', type=click.Path(exists=True),
              help='Configuration file (YAML or JSON)')
@click.option('--preset', type=click.Choice(['small', 'large', 'sensitive', 'conservative', 'hpc', 'development']),
              help='Use predefined configuration preset')
@click.option('--generate-config', type=click.Path(),
              help='Generate a configuration file with current settings and exit')
@click.option('--list-presets', is_flag=True,
              help='List available configuration presets and exit')
@click.option('--assembly', '-a', required=True, type=click.Path(exists=True),
              help='Input assembly file in FASTA format')
@click.option('--bam', '-b', type=click.Path(exists=True),
              help='Input BAM file with aligned reads')
@click.option('--reads1', '-1', type=click.Path(exists=True),
              help='Forward reads file (FASTQ/FASTA, optionally gzipped)')
@click.option('--reads2', '-2', type=click.Path(exists=True),
              help='Reverse reads file (FASTQ/FASTA, optionally gzipped)')
@click.option('--reads', '-r', type=click.Path(exists=True),
              help='Single-end reads file (FASTQ/FASTA, optionally gzipped)')
@click.option('--reads-dir', type=click.Path(exists=True, file_okay=False),
              help='Directory containing multiple read files')
@click.option('--reads-pattern', default='*_R{1,2}.fastq.gz',
              help='Pattern for finding read files in reads-dir [default: *_R{1,2}.fastq.gz]')
@click.option('--multi-sample-mode', type=click.Choice(['separate', 'merged', 'batch']), 
              default='separate', help='How to process multiple samples [default: separate]')
@click.option('--max-workers', default=4, type=int,
              help='Maximum parallel workers for multi-sample processing [default: 4]')
@click.option('--batch-size', default=5, type=int,
              help='Batch size for batch processing mode [default: 5]')
@click.option('--parallel/--no-parallel', default=True,
              help='Enable parallel processing of samples [default: enabled]')
@click.option('--out', '-o', required=True, type=click.Path(),
              help='Output directory')
@click.option('--reference', type=click.Path(exists=True),
              help='Reference database for taxonomic classification (optional)')
@click.option('--min-contig-length', default=1000, type=int,
              help='Minimum contig length to analyze [default: 1000]')
@click.option('--min-coverage', default=5.0, type=float,
              help='Minimum coverage threshold [default: 5.0]')
@click.option('--coverage-fold-change', default=2.0, type=float,
              help='Minimum coverage fold change to consider [default: 2.0]')
@click.option('--gc-content-threshold', default=0.1, type=float,
              help='Minimum GC content difference threshold [default: 0.1]')
@click.option('--kmer-distance-threshold', default=0.3, type=float,
              help='Minimum k-mer distance threshold [default: 0.3]')
@click.option('--confidence-threshold', default=0.5, type=float,
              help='Minimum confidence threshold for splitting [default: 0.5]')
@click.option('--min-split-length', default=500, type=int,
              help='Minimum length for split contigs [default: 500]')
@click.option('--threads', '-t', default=1, type=int,
              help='Number of threads to use [default: 1]')
@click.option('--sensitivity', type=click.Choice(['low', 'medium', 'high']), default='medium',
              help='Detection sensitivity level [default: medium]')
@click.option('--split-technical/--no-split-technical', default=True,
              help='Split technical artifacts [default: enabled]')
@click.option('--split-pcr/--no-split-pcr', default=True,
              help='Split PCR chimeras [default: enabled]')
@click.option('--preserve-biological/--no-preserve-biological', default=True,
              help='Preserve biological recombination [default: enabled]')
@click.option('--generate-report/--no-generate-report', default=True,
              help='Generate interactive HTML report [default: enabled]')
@click.option('--keep-intermediates/--no-keep-intermediates', default=False,
              help='Keep intermediate files [default: disabled]')
@click.option('--log-level', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']), 
              default='INFO', help='Logging level [default: INFO]')
@click.version_option(version=__version__)
def main(**kwargs):
    """
    Chimeric Detective: Detect and resolve chimeric contigs in viral metagenomic assemblies.
    
    BASIC USAGE:
    
    With BAM file:
        chimeric_detective -a assembly.fasta -b aligned.bam -o results/
    
    With paired-end reads:
        chimeric_detective -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/
    
    With single-end reads:
        chimeric_detective -a assembly.fasta -r reads.fastq.gz -o results/
    
    CONFIGURATION USAGE:
    
    With configuration file:
        chimeric_detective -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -c config.yaml -o results/
    
    With preset configuration:
        chimeric_detective -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz --preset large -o results/
    
    Generate configuration file:
        chimeric_detective --generate-config myconfig.yaml
    
    List available presets:
        chimeric_detective --list-presets
    
    MULTI-SAMPLE USAGE:
    
    With multiple samples (separate analysis):
        chimeric_detective -a assembly.fasta --reads-dir /path/to/reads/ --reads-pattern "*_R{1,2}.fastq.gz" -o results/
    
    With multiple samples (merged analysis):
        chimeric_detective -a assembly.fasta --reads-dir /path/to/reads/ --multi-sample-mode merged -o results/
    
    With multiple samples (parallel processing):
        chimeric_detective -a assembly.fasta --reads-dir /path/to/reads/ --max-workers 8 --parallel -o results/
    
    ADVANCED USAGE:
    
    With reference database and custom parameters:
        chimeric_detective -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz \\
                          --reference viral_db.fasta --sensitivity high \\
                          --min-coverage 10 --threads 8 -o results/
    """
    try:
        # Handle special commands first
        if kwargs.get('list_presets'):
            _list_presets()
            return
        
        if kwargs.get('generate_config'):
            _generate_config(kwargs)
            return
        
        # Initialize configuration manager
        config_manager = ConfigManager()
        
        # Load configuration from various sources (in order of precedence)
        # 1. Auto-load from default locations
        auto_config_path = config_manager.auto_load_config()
        if auto_config_path:
            click.echo(f"📁 Auto-loaded configuration from: {auto_config_path}")
        
        # 2. Load preset if specified
        if kwargs.get('preset'):
            config_manager.load_preset(kwargs['preset'])
            click.echo(f"⚙️  Applied preset: {kwargs['preset']}")
        
        # 3. Load config file if specified (overrides preset)
        if kwargs.get('config'):
            config_manager.load_config(kwargs['config'])
            click.echo(f"📁 Loaded configuration from: {kwargs['config']}")
        
        # 4. Apply environment variable overrides
        config_manager.apply_env_overrides()
        
        # 5. Apply CLI argument overrides (highest precedence)
        _apply_cli_overrides(config_manager, kwargs)
        
        # Validate final configuration
        config_manager.validate_config()
        
        # Merge configuration with CLI arguments
        merged_kwargs = config_manager.to_cli_args()
        # CLI arguments override config (keep original CLI values if provided)
        for key, value in kwargs.items():
            if value is not None and key not in ['config', 'preset', 'generate_config', 'list_presets']:
                merged_kwargs[key] = value
        # Create output directory first (needed for log file)
        create_output_directory(merged_kwargs['out'])
        
        # Setup logging
        log_file = Path(merged_kwargs['out']) / 'chimeric_detective.log' if merged_kwargs['out'] else None
        logger = setup_logging(merged_kwargs['log_level'], str(log_file) if log_file else None)
        
        logger.info(f"Starting Chimeric Detective v{__version__}")
        logger.info(f"Command: {' '.join(sys.argv)}")
        
        # Validate inputs and setup
        _validate_inputs(**merged_kwargs)
        _check_dependencies(logger)
        
        # Adjust parameters based on sensitivity
        detector_params, analyzer_params = _adjust_sensitivity_parameters(merged_kwargs)
        
        # Apply sensitivity adjustments to merged config
        merged_kwargs.update(detector_params)
        merged_kwargs.update(analyzer_params)
        
        # Check if this is multi-sample processing
        if merged_kwargs.get('reads_dir'):
            _run_multi_sample_pipeline(logger, **merged_kwargs)
        else:
            _run_pipeline(logger, **merged_kwargs)
        
        logger.info("Analysis completed successfully!")
        click.echo(f"✅ Results written to: {merged_kwargs['out']}")
        
    except Exception as e:
        if 'logger' in locals():
            logger.error(f"Pipeline failed: {str(e)}")
        click.echo(f"❌ Error: {str(e)}", err=True)
        sys.exit(1)


def _validate_inputs(**kwargs):
    """Validate input parameters."""
    
    # Check assembly file
    validate_file_exists(kwargs['assembly'], "Assembly file")
    
    # Check read inputs
    has_bam = kwargs.get('bam') is not None
    has_reads1 = kwargs.get('reads1') is not None
    has_reads2 = kwargs.get('reads2') is not None
    has_single_reads = kwargs.get('reads') is not None
    has_reads_dir = kwargs.get('reads_dir') is not None
    
    read_input_count = sum([has_bam, has_reads1 or has_reads2, has_single_reads, has_reads_dir])
    
    if read_input_count == 0:
        raise click.BadParameter("Must provide either --bam, --reads1/--reads2, --reads, or --reads-dir")
    
    if read_input_count > 1:
        raise click.BadParameter("Can only provide one type of read input")
    
    # Validate read files if provided
    if has_bam:
        validate_file_exists(kwargs['bam'], "BAM file")
    
    if has_reads1:
        validate_file_exists(kwargs['reads1'], "Forward reads file")
        if has_reads2:
            validate_file_exists(kwargs['reads2'], "Reverse reads file")
    
    if has_single_reads:
        validate_file_exists(kwargs['reads'], "Single reads file")
    
    if has_reads_dir:
        if not os.path.isdir(kwargs['reads_dir']):
            raise click.BadParameter(f"Reads directory does not exist: {kwargs['reads_dir']}")
    
    # Validate reference database if provided
    if kwargs.get('reference'):
        validate_file_exists(kwargs['reference'], "Reference database")
    
    # Validate parameter ranges
    if not 0 < kwargs['confidence_threshold'] <= 1:
        raise click.BadParameter("Confidence threshold must be between 0 and 1")
    
    if kwargs['min_contig_length'] < 100:
        raise click.BadParameter("Minimum contig length must be at least 100")
    
    if kwargs['min_coverage'] < 0:
        raise click.BadParameter("Minimum coverage must be non-negative")


def _check_dependencies(logger):
    """Check for required external dependencies."""
    
    tools = check_external_tools()
    
    missing_tools = [tool for tool, available in tools.items() if not available]
    
    if missing_tools:
        logger.warning(f"Missing optional tools: {', '.join(missing_tools)}")
        logger.warning("Some functionality may be limited")
    else:
        logger.info("All external tools found")


def _adjust_sensitivity_parameters(kwargs):
    """Adjust detection parameters based on sensitivity level."""
    
    sensitivity = kwargs['sensitivity']
    
    if sensitivity == 'low':
        detector_params = {
            'coverage_fold_change': kwargs['coverage_fold_change'] * 1.5,
            'gc_content_threshold': kwargs['gc_content_threshold'] * 1.5,
            'kmer_distance_threshold': kwargs['kmer_distance_threshold'] * 1.5,
        }
        analyzer_params = {
            'min_blast_identity': 85.0,
            'min_blast_coverage': 60.0,
        }
    elif sensitivity == 'high':
        detector_params = {
            'coverage_fold_change': kwargs['coverage_fold_change'] * 0.7,
            'gc_content_threshold': kwargs['gc_content_threshold'] * 0.7,
            'kmer_distance_threshold': kwargs['kmer_distance_threshold'] * 0.7,
        }
        analyzer_params = {
            'min_blast_identity': 75.0,
            'min_blast_coverage': 40.0,
        }
    else:  # medium
        detector_params = {
            'coverage_fold_change': kwargs['coverage_fold_change'],
            'gc_content_threshold': kwargs['gc_content_threshold'],
            'kmer_distance_threshold': kwargs['kmer_distance_threshold'],
        }
        analyzer_params = {
            'min_blast_identity': 80.0,
            'min_blast_coverage': 50.0,
        }
    
    return detector_params, analyzer_params


def _run_multi_sample_pipeline(logger, **kwargs):
    """Run the multi-sample chimera detection and resolution pipeline."""
    
    logger.info("🔬 Starting multi-sample analysis")
    logger.info(f"Processing mode: {kwargs['multi_sample_mode']}")
    logger.info(f"Max workers: {kwargs['max_workers']}")
    
    # Initialize multi-sample processor
    processor = MultiSampleProcessor(
        processing_mode=kwargs['multi_sample_mode'],
        max_workers=kwargs['max_workers'],
        log_level=kwargs['log_level']
    )
    
    # Process all samples - remove conflicting params from kwargs
    processing_kwargs = kwargs.copy()
    processing_kwargs.pop('assembly', None)
    processing_kwargs.pop('reads_dir', None) 
    processing_kwargs.pop('reads_pattern', None)
    processing_kwargs.pop('out', None)
    processing_kwargs.pop('reads1', None)
    processing_kwargs.pop('reads2', None)
    processing_kwargs.pop('reads', None)
    processing_kwargs.pop('bam', None)
    
    results = processor.process_samples_directory(
        assembly_file=kwargs['assembly'],
        reads_dir=kwargs['reads_dir'],
        reads_pattern=kwargs['reads_pattern'],
        output_dir=kwargs['out'],
        **processing_kwargs
    )
    
    # Print summary
    _print_multi_sample_summary(results, kwargs['multi_sample_mode'])


def _run_pipeline(logger, **kwargs):
    """Run the complete chimera detection and resolution pipeline."""
    
    temp_dir = None
    if not kwargs['keep_intermediates']:
        temp_dir = tempfile.mkdtemp(prefix="chimeric_detective_")
    
    try:
        output_dir = Path(kwargs['out'])
        
        # Step 1: Chimera Detection
        logger.info("🔍 Step 1: Detecting chimeric contigs...")
        
        detector = ChimeraDetector(
            min_contig_length=kwargs['min_contig_length'],
            min_coverage=kwargs['min_coverage'],
            coverage_fold_change=kwargs['coverage_fold_change'],
            gc_content_threshold=kwargs['gc_content_threshold'],
            kmer_distance_threshold=kwargs['kmer_distance_threshold'],
            log_level=kwargs['log_level']
        )
        
        # Prepare reads input
        bam_file, reads1, reads2 = _prepare_reads_input(kwargs, temp_dir)
        
        candidates = detector.detect_chimeras(
            assembly_file=kwargs['assembly'],
            bam_file=bam_file,
            reads1=reads1,
            reads2=reads2,
            temp_dir=temp_dir
        )
        
        if not candidates:
            logger.info("No chimeric contigs detected")
            click.echo("✅ No chimeric contigs detected. Assembly appears clean.")
            return
        
        logger.info(f"Detected {len(candidates)} chimera candidates")
        
        # Step 2: Chimera Analysis
        logger.info("🧬 Step 2: Analyzing and classifying chimeras...")
        
        analyzer = ChimeraAnalyzer(
            reference_db=kwargs['reference'],
            log_level=kwargs['log_level']
        )
        
        analyses = analyzer.analyze_chimeras(
            candidates=candidates,
            assembly_file=kwargs['assembly'],
            bam_file=bam_file
        )
        
        # Step 3: Chimera Resolution
        logger.info("✂️ Step 3: Resolving chimeras and creating cleaned assembly...")
        
        resolver = ChimeraResolver(
            split_technical=kwargs['split_technical'],
            split_pcr=kwargs['split_pcr'],
            preserve_biological=kwargs['preserve_biological'],
            min_split_length=kwargs['min_split_length'],
            confidence_threshold=kwargs['confidence_threshold'],
            log_level=kwargs['log_level']
        )
        
        output_files = resolver.resolve_chimeras(
            analyses=analyses,
            assembly_file=kwargs['assembly'],
            output_dir=str(output_dir)
        )
        
        # Step 4: Generate Report
        if kwargs['generate_report']:
            logger.info("📊 Step 4: Generating interactive report...")
            
            visualizer = ChimeraVisualizer(log_level=kwargs['log_level'])
            
            report_path = visualizer.create_report(
                analyses=analyses,
                decisions=resolver.splitting_decisions,
                output_dir=str(output_dir),
                assembly_file=kwargs['assembly']
            )
            
            logger.info(f"Interactive report generated: {report_path}")
        
        # Print summary
        _print_summary(analyses, resolver.splitting_decisions, output_files)
        
    finally:
        # Cleanup temporary files if not keeping intermediates
        if temp_dir and not kwargs['keep_intermediates']:
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)


def _prepare_reads_input(kwargs, temp_dir) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Prepare reads input for the pipeline."""
    
    if kwargs.get('bam'):
        return kwargs['bam'], None, None
    
    elif kwargs.get('reads1'):
        return None, kwargs['reads1'], kwargs.get('reads2')
    
    elif kwargs.get('reads'):
        return None, kwargs['reads'], None
    
    elif kwargs.get('reads_dir'):
        # Multiple samples - handle via multi-sample processor
        return None, None, None  # Will be handled by multi-sample processor
    
    else:
        raise click.BadParameter("No valid read input provided")


def _print_summary(analyses, decisions, output_files):
    """Print a summary of the analysis results."""
    
    click.echo("\n" + "="*60)
    click.echo("📋 SUMMARY")
    click.echo("="*60)
    
    # Chimera detection summary
    chimera_types = {}
    for analysis in analyses:
        chimera_type = analysis.chimera_type
        chimera_types[chimera_type] = chimera_types.get(chimera_type, 0) + 1
    
    click.echo(f"🔍 Detected {len(analyses)} chimeric contigs:")
    for chimera_type, count in chimera_types.items():
        click.echo(f"   - {chimera_type}: {count}")
    
    # Decision summary
    decision_types = {}
    for decision in decisions:
        action = decision.action
        decision_types[action] = decision_types.get(action, 0) + 1
    
    click.echo(f"\n✂️ Resolution decisions:")
    for action, count in decision_types.items():
        action_icon = "🔪" if action == "split" else "💾" if action == "preserve" else "🏷️"
        click.echo(f"   {action_icon} {action}: {count} contigs")
    
    # Output files
    click.echo(f"\n📁 Output files:")
    click.echo(f"   - Cleaned assembly: {output_files.get('cleaned_assembly', 'N/A')}")
    click.echo(f"   - Splitting decisions: {output_files.get('splitting_decisions', 'N/A')}")
    click.echo(f"   - Detailed results: {output_files.get('results_json', 'N/A')}")
    
    if 'chimeric_detective_report.html' in str(output_files.get('results_json', '')):
        report_path = str(output_files['results_json']).replace('results.json', 'report.html')
        click.echo(f"   - Interactive report: {report_path}")
    
    click.echo("\n🎉 Analysis complete! Check the output directory for detailed results.")


def _print_multi_sample_summary(results: Dict[str, str], processing_mode: str):
    """Print summary for multi-sample analysis."""
    
    click.echo("\n" + "="*60)
    click.echo("📋 MULTI-SAMPLE SUMMARY")
    click.echo("="*60)
    
    click.echo(f"🔬 Processing mode: {processing_mode}")
    click.echo(f"📁 Total samples processed: {len(results)}")
    
    if processing_mode == "separate":
        click.echo(f"\n📊 Individual sample results:")
        for sample_name, output_dir in results.items():
            click.echo(f"   - {sample_name}: {output_dir}")
            
        click.echo(f"\n📈 View combined results:")
        # Assuming first result gives us the parent directory
        if results:
            parent_dir = Path(list(results.values())[0]).parent
            multi_report = parent_dir / "multi_sample_report.html"
            if multi_report.exists():
                click.echo(f"   🌐 Multi-sample report: {multi_report}")
            
            summary_table = parent_dir / "multi_sample_summary.tsv"
            if summary_table.exists():
                click.echo(f"   📋 Summary table: {summary_table}")
    
    elif processing_mode == "merged":
        click.echo(f"\n📊 Merged analysis results:")
        for sample_name, output_dir in results.items():
            click.echo(f"   - Combined results: {output_dir}")
    
    click.echo(f"\n💡 Tips for multi-sample analysis:")
    click.echo(f"   - Use 'separate' mode to analyze each sample independently")
    click.echo(f"   - Use 'merged' mode to combine all reads for higher coverage")
    click.echo(f"   - Use 'batch' mode for memory-efficient processing of many samples")
    click.echo(f"   - Adjust --max-workers based on your system resources")
    
    click.echo("\n🎉 Multi-sample analysis complete!")


def _list_presets():
    """List available configuration presets."""
    click.echo("🎛️  Available Configuration Presets:\n")
    
    presets = ConfigManager.list_presets()
    for name, description in presets.items():
        click.echo(f"  {name:12} - {description}")
    
    click.echo(f"\nUsage: chimeric_detective --preset <name> [other options]")
    click.echo(f"Example: chimeric_detective --preset large -a assembly.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o results/")


def _generate_config(kwargs):
    """Generate a configuration file with current settings."""
    config_path = kwargs['generate_config']
    
    # Initialize config manager and apply any settings
    config_manager = ConfigManager()
    
    # Apply preset if specified
    if kwargs.get('preset'):
        config_manager.load_preset(kwargs['preset'])
    
    # Apply CLI overrides
    _apply_cli_overrides(config_manager, kwargs)
    
    # Determine format from file extension
    config_path = Path(config_path)
    if config_path.suffix.lower() in ['.yaml', '.yml']:
        format_type = 'yaml'
    elif config_path.suffix.lower() == '.json':
        format_type = 'json'
    else:
        # Default to YAML if no extension
        format_type = 'yaml'
        if not config_path.suffix:
            config_path = config_path.with_suffix('.yaml')
    
    try:
        config_manager.save_config(config_path, format_type)
        click.echo(f"✅ Configuration file generated: {config_path}")
        click.echo(f"📝 Edit the file to customize settings, then use: --config {config_path}")
    except Exception as e:
        click.echo(f"❌ Failed to generate configuration file: {e}", err=True)
        sys.exit(1)


def _apply_cli_overrides(config_manager: ConfigManager, kwargs):
    """Apply CLI argument overrides to configuration."""
    # Detection parameters
    if kwargs.get('min_contig_length') is not None:
        config_manager.config.detection.min_contig_length = kwargs['min_contig_length']
    if kwargs.get('min_coverage') is not None:
        config_manager.config.detection.min_coverage = kwargs['min_coverage']
    if kwargs.get('coverage_fold_change') is not None:
        config_manager.config.detection.coverage_fold_change = kwargs['coverage_fold_change']
    if kwargs.get('gc_content_threshold') is not None:
        config_manager.config.detection.gc_content_threshold = kwargs['gc_content_threshold']
    if kwargs.get('kmer_distance_threshold') is not None:
        config_manager.config.detection.kmer_distance_threshold = kwargs['kmer_distance_threshold']
    if kwargs.get('confidence_threshold') is not None:
        config_manager.config.detection.confidence_threshold = kwargs['confidence_threshold']
    if kwargs.get('min_split_length') is not None:
        config_manager.config.detection.min_split_length = kwargs['min_split_length']
    
    # Processing parameters
    if kwargs.get('threads') is not None:
        config_manager.config.processing.threads = kwargs['threads']
    if kwargs.get('max_workers') is not None:
        config_manager.config.processing.max_workers = kwargs['max_workers']
    if kwargs.get('parallel') is not None:
        config_manager.config.processing.parallel = kwargs['parallel']
    if kwargs.get('keep_intermediates') is not None:
        config_manager.config.processing.keep_intermediates = kwargs['keep_intermediates']
    if kwargs.get('batch_size') is not None:
        config_manager.config.processing.batch_size = kwargs['batch_size']
    
    # Output parameters
    if kwargs.get('generate_report') is not None:
        config_manager.config.output.generate_report = kwargs['generate_report']
    if kwargs.get('log_level') is not None:
        config_manager.config.output.log_level = kwargs['log_level']
    
    # Behavior parameters
    if kwargs.get('split_technical') is not None:
        config_manager.config.behavior.split_technical = kwargs['split_technical']
    if kwargs.get('split_pcr') is not None:
        config_manager.config.behavior.split_pcr = kwargs['split_pcr']
    if kwargs.get('preserve_biological') is not None:
        config_manager.config.behavior.preserve_biological = kwargs['preserve_biological']
    if kwargs.get('sensitivity') is not None:
        config_manager.config.behavior.sensitivity = kwargs['sensitivity']
    
    # Multi-sample parameters
    if kwargs.get('multi_sample_mode') is not None:
        config_manager.config.multi_sample.multi_sample_mode = kwargs['multi_sample_mode']
    if kwargs.get('reads_pattern') is not None:
        config_manager.config.multi_sample.reads_pattern = kwargs['reads_pattern']


if __name__ == '__main__':
    main()