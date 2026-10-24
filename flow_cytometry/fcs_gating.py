"""Gating and per-well statistics for flow cytometry FCS files (Fig 5G).

Pruned from the lab's fcs_analysis.py to the steps used for Fig 5G: load FCS files with a metadata
table, compensate with the acquisition spillover matrix, gate cells (FSC/SSC rectangle), singlets
(FSC width percentiles) and fluorophore-positive events, combine gates, and export per-well channel
statistics. The functions below are unchanged from the original module.

    from fcs_gating import FCAnalyzer
    analyzer = FCAnalyzer(output_dir='out')
    analyzer.load_files(fcs_paths, metadata_csv='metadata.csv')
    analyzer.apply_compensation()
    analyzer.gate_cells(fsc='FSC-H', ssc='SSC-H', fsc_min=..., ...)
    analyzer.gate_singlets(area='FSC-A', height='FSC-H', width='FSC-W', method='width')
    analyzer.gate_fluorophore(channel='Comp-BL1-A', threshold=..., positive_name='GFP+')
    stats = analyzer.get_channel_statistics(channels=[...], gate=['GFP+'], stats=['median'])
"""

import warnings
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Union

import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import inv
import fcsparser


def parse_spillover_matrix(spillover_str: str) -> tuple[list[str], np.ndarray]:
    """
    Parse the $SPILLOVER or $SPILL string from FCS metadata.
    
    Format is typically: n,channel1,channel2,...,channelN,val1,val2,...,valN*N
    
    Returns:
        Tuple of (channel_names, spillover_matrix)
    """
    parts = spillover_str.split(',')
    n = int(parts[0])
    channels = parts[1:n+1]
    values = [float(v) for v in parts[n+1:]]
    
    if len(values) != n * n:
        raise ValueError(f"Expected {n*n} matrix values, got {len(values)}")
    
    matrix = np.array(values).reshape(n, n)
    return channels, matrix


@dataclass
class GateResult:
    """Stores results from a single gating operation."""
    name: str
    parent_count: int
    gated_count: int
    percent_of_parent: float
    percent_of_total: float
    mask: np.ndarray
    gate_params: dict = field(default_factory=dict)


@dataclass 
class FileAnalysis:
    """Stores complete analysis results for a single FCS file."""
    filepath: str
    filename: str
    metadata: dict
    raw_data: pd.DataFrame
    channel_map: dict  # original -> renamed
    total_events: int
    gates: dict = field(default_factory=dict)  # gate_name -> GateResult
    current_mask: np.ndarray = None  # cumulative mask after all gates
    
    def get_gated_data(self, gate_name: str = None) -> pd.DataFrame:
        """Return data filtered by specified gate or current cumulative gate."""
        if gate_name and gate_name in self.gates:
            return self.raw_data[self.gates[gate_name].mask]
        elif self.current_mask is not None:
            return self.raw_data[self.current_mask]
        return self.raw_data


class FCAnalyzer:
    """
    Batch gating of FCS files.

    Workflow:
        1. Load FCS files (with a metadata CSV)
        2. Compensate with the embedded spillover matrix
        3. Apply sequential gates (cells, singlets) and fluorophore gates
        4. Export per-file channel statistics
    """

    def __init__(self, output_dir: str = './fcs_results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.files: dict[str, FileAnalysis] = {}  # filename -> FileAnalysis
        self.channel_rename_map: dict[str, str] = {}  # global rename mapping
        self.gate_sequence: list[str] = []  # ordered list of applied gates
        self.metadata_df: Optional[pd.DataFrame] = None
        
        # Figure settings
        self.fig_dpi = 150
        self.fig_format = 'png'
        self.grid_rows = 3
        self.grid_cols = 4
        
    def load_files(
        self, 
        path: Union[str, list[str]], 
        metadata_csv: str = None,
        file_column: str = 'filename',
        recursive: bool = False
    ) -> 'FCAnalyzer':
        """
        Load FCS files from directory or list of paths.
        
        Args:
            path: Directory path or list of file paths
            metadata_csv: Optional CSV with file metadata. Must have a column 
                         matching filenames.
            file_column: Column name in metadata CSV containing filenames
            recursive: If True, search subdirectories
            
        Returns:
            self for method chaining
        """
        # Load metadata if provided
        if metadata_csv:
            self.metadata_df = pd.read_csv(metadata_csv)
            if file_column not in self.metadata_df.columns:
                raise ValueError(f"Column '{file_column}' not found in metadata CSV")
            self.metadata_df = self.metadata_df.set_index(file_column)
        
        # Gather file paths
        fcs_paths = []
        if isinstance(path, list):
            fcs_paths = [Path(p) for p in path]
        else:
            path = Path(path)
            if path.is_file():
                fcs_paths = [path]
            else:
                pattern = '**/*.fcs' if recursive else '*.fcs'
                fcs_paths = list(path.glob(pattern))
                # Also check for .FCS extension
                fcs_paths.extend(path.glob(pattern.replace('.fcs', '.FCS')))
        
        print(f"Found {len(fcs_paths)} FCS files")
        
        # Load each file
        for fcs_path in fcs_paths:
            try:
                meta, data = fcsparser.parse(str(fcs_path), reformat_meta=True)
                filename = fcs_path.name
                
                # Get file-specific metadata from CSV if available
                file_metadata = {}
                if self.metadata_df is not None and filename in self.metadata_df.index:
                    file_metadata = self.metadata_df.loc[filename].to_dict()
                
                # Store FCS header metadata
                file_metadata['_fcs_meta'] = {
                    k: v for k, v in meta.items() 
                    if isinstance(v, (str, int, float))
                }
                
                self.files[filename] = FileAnalysis(
                    filepath=str(fcs_path),
                    filename=filename,
                    metadata=file_metadata,
                    raw_data=data,
                    channel_map={c: c for c in data.columns},
                    total_events=len(data),
                    current_mask=np.ones(len(data), dtype=bool)
                )
                print(f"  Loaded: {filename} ({len(data):,} events)")
                
            except Exception as e:
                warnings.warn(f"Failed to load {fcs_path}: {e}")
        
        return self
    
    def get_spillover_matrix(self, filename: str = None) -> tuple[list[str], np.ndarray]:
        """
        Extract spillover matrix from FCS file metadata.
        
        Args:
            filename: Specific file to get matrix from. If None, uses first file.
            
        Returns:
            Tuple of (channel_names, spillover_matrix)
        """
        if filename is None:
            filename = list(self.files.keys())[0]
        
        fa = self.files[filename]
        meta = fa.metadata.get('_fcs_meta', {})
        
        # Look for spillover matrix in common locations
        spillover_str = None
        for key in ['$SPILLOVER', 'SPILLOVER', '$SPILL', 'SPILL', '$COMP', 'COMP']:
            if key in meta:
                spillover_str = meta[key]
                break
        
        if spillover_str is None:
            raise ValueError(
                f"No spillover matrix found in {filename}. "
                f"Available metadata keys: {list(meta.keys())}"
            )
        
        return parse_spillover_matrix(spillover_str)
    
    def apply_compensation(
        self,
        matrix: np.ndarray = None,
        channels: list[str] = None,
        use_file: str = None,
        prefix: str = 'Comp-',
        replace: bool = False
    ) -> 'FCAnalyzer':
        """
        Apply compensation to fluorescence channels.
        
        Compensation corrects for spectral overlap between fluorophores.
        Uses the spillover matrix embedded in FCS files (from acquisition)
        or a custom matrix you provide.
        
        Args:
            matrix: Optional custom compensation matrix (NxN numpy array).
                   If None, extracts from FCS file metadata.
            channels: Channel names corresponding to matrix rows/columns.
                     Required if providing custom matrix.
            use_file: Which file to extract spillover matrix from (default: first file).
                     Ignored if custom matrix provided.
            prefix: Prefix for compensated channel names (default: 'Comp-')
            replace: If True, replace original channels. If False, add new columns.
            
        Returns:
            self for method chaining
            
        Example:
            # Use embedded matrix from files
            analyzer.apply_compensation()
            
            # Use custom matrix
            analyzer.apply_compensation(
                matrix=my_matrix,
                channels=['FITC-A', 'PE-A', 'APC-A']
            )
        """
        # Get spillover matrix
        if matrix is None:
            channels, spillover = self.get_spillover_matrix(use_file)
            # Invert spillover to get compensation matrix
            comp_matrix = inv(spillover)
            print(f"Extracted spillover matrix for channels: {channels}")
        else:
            if channels is None:
                raise ValueError("Must provide channel names when using custom matrix")
            comp_matrix = matrix
            spillover = None
        
        # Store for reference
        self._compensation_matrix = comp_matrix
        self._compensation_channels = channels
        
        # Apply to each file
        for filename, fa in self.files.items():
            # Check which channels are present in this file
            present_channels = [ch for ch in channels if ch in fa.raw_data.columns]
            missing = set(channels) - set(present_channels)
            
            if missing:
                warnings.warn(
                    f"{filename}: Missing compensation channels {missing}, skipping"
                )
                continue
            
            if len(present_channels) != len(channels):
                warnings.warn(
                    f"{filename}: Only {len(present_channels)}/{len(channels)} "
                    f"compensation channels present"
                )
                continue
            
            # Extract fluorescence data matrix (events x channels)
            fluor_data = fa.raw_data[channels].values
            
            # Apply compensation: compensated = raw @ compensation_matrix
            compensated = fluor_data @ comp_matrix
            
            # Add or replace channels
            for i, ch in enumerate(channels):
                new_name = ch if replace else f"{prefix}{ch}"
                fa.raw_data[new_name] = compensated[:, i]
                fa.channel_map[f"comp_{ch}"] = new_name
        
        n_channels = len(channels)
        action = "Replaced" if replace else "Added"
        print(f"{action} {n_channels} compensated channels across {len(self.files)} files")
        
        return self
    
    def get_common_channels(self) -> list[str]:
        """Return channels common to all loaded files."""
        if not self.files:
            return []
        channel_sets = [set(fa.raw_data.columns) for fa in self.files.values()]
        return list(set.intersection(*channel_sets))
    
    def _validate_channels(self, channels: list[str]) -> None:
        """Check that channels exist in all files."""
        common = set(self.get_common_channels())
        missing = set(channels) - common
        if missing:
            raise ValueError(
                f"Channels not found in all files: {missing}. "
                f"Available common channels: {common}"
            )
    
    def gate_cells(
        self,
        fsc: str = 'FSC-A',
        ssc: str = 'SSC-A',
        fsc_min: float = None,
        fsc_max: float = None,
        ssc_min: float = None,
        ssc_max: float = None,
        method: str = 'rectangle',
        percentile_low: float = 1,
        percentile_high: float = 99,
        name: str = 'Cells'
    ) -> 'FCAnalyzer':
        """
        Gate cell-like objects based on forward and side scatter.
        
        Args:
            fsc: Forward scatter channel name
            ssc: Side scatter channel name
            fsc_min/max: Manual FSC thresholds
            ssc_min/max: Manual SSC thresholds
            method: 'rectangle' for manual bounds, 'auto' for percentile-based
            percentile_low/high: For auto method, exclude events outside these percentiles
            name: Name for this gate
        """
        self._validate_channels([fsc, ssc])
        self.gate_sequence.append(name)
        
        for filename, fa in self.files.items():
            parent_mask = fa.current_mask.copy()
            parent_count = parent_mask.sum()
            
            data = fa.raw_data
            
            if method == 'auto':
                # Use percentiles from parent population
                parent_data = data.loc[parent_mask]
                if parent_count == 0:
                    # No parent events - use data range as fallback
                    fsc_lo = data[fsc].min()
                    fsc_hi = data[fsc].max()
                    ssc_lo = data[ssc].min()
                    ssc_hi = data[ssc].max()
                else:
                    fsc_lo = np.percentile(parent_data[fsc], percentile_low)
                    fsc_hi = np.percentile(parent_data[fsc], percentile_high)
                    ssc_lo = np.percentile(parent_data[ssc], percentile_low)
                    ssc_hi = np.percentile(parent_data[ssc], percentile_high)
            else:
                fsc_lo = fsc_min if fsc_min is not None else data[fsc].min()
                fsc_hi = fsc_max if fsc_max is not None else data[fsc].max()
                ssc_lo = ssc_min if ssc_min is not None else data[ssc].min()
                ssc_hi = ssc_max if ssc_max is not None else data[ssc].max()
            
            gate_mask = (
                parent_mask &
                (data[fsc] >= fsc_lo) & (data[fsc] <= fsc_hi) &
                (data[ssc] >= ssc_lo) & (data[ssc] <= ssc_hi)
            )
            
            gated_count = gate_mask.sum()
            
            fa.gates[name] = GateResult(
                name=name,
                parent_count=parent_count,
                gated_count=gated_count,
                percent_of_parent=100 * gated_count / parent_count if parent_count > 0 else 0,
                percent_of_total=100 * gated_count / fa.total_events,
                mask=gate_mask,
                gate_params={
                    'fsc': fsc, 'ssc': ssc,
                    'fsc_min': fsc_lo, 'fsc_max': fsc_hi,
                    'ssc_min': ssc_lo, 'ssc_max': ssc_hi,
                    'method': method
                }
            )
            fa.current_mask = gate_mask
        
        self._generate_gate_figures(name, 'scatter', x_channel=fsc, y_channel=ssc)
        return self
    
    def gate_singlets(
        self,
        area: str = 'FSC-A',
        height: str = 'FSC-H',
        width: str = None,
        method: str = 'diagonal',
        slope_tolerance: float = 0.15,
        fit_percentile: float = 40,
        width_percentile_low: float = 5,
        width_percentile_high: float = 95,
        name: str = 'Singlets'
    ) -> 'FCAnalyzer':
        """
        Gate singlets using area/height ratio or width gating.
        
        Args:
            area: Area channel (e.g., FSC-A)
            height: Height channel (e.g., FSC-H)
            width: Width channel (optional, e.g., FSC-W)
            method: 'diagonal' uses A vs H correlation, 'width' uses W gating
            slope_tolerance: For diagonal method, tolerance around ideal slope
            fit_percentile: For diagonal method, fit only events below this 
                           percentile of Area. Lower values = less doublet 
                           contamination in fit. Default 40 (bottom 40%).
            width_percentile_low/high: For width method, percentile bounds
            name: Name for this gate
        """
        channels = [area, height]
        if width:
            channels.append(width)
        self._validate_channels(channels)
        self.gate_sequence.append(name)
        
        for filename, fa in self.files.items():
            parent_mask = fa.current_mask.copy()
            parent_count = parent_mask.sum()
            data = fa.raw_data
            parent_data = data.loc[parent_mask]
            
            # Handle empty parent population
            if parent_count == 0:
                gate_mask = parent_mask.copy()  # All False
                gate_params = {
                    'area': area, 'height': height, 'method': method,
                    'note': 'empty_parent_population'
                }
                if width:
                    gate_params['width'] = width
                    
            elif method == 'diagonal':
                # Fit linear relationship using only LOWER portion of data
                # This avoids doublets/multiplets skewing the fit
                a_vals = parent_data[area].values
                h_vals = parent_data[height].values
                
                # Only fit events in the lower percentile of Area
                # (small cells are almost always singlets)
                area_threshold = np.percentile(a_vals[a_vals > 0], fit_percentile)
                fit_mask = (a_vals > 0) & (h_vals > 0) & (a_vals <= area_threshold)
                
                if fit_mask.sum() > 100:
                    slope, intercept, _, _, _ = stats.linregress(
                        a_vals[fit_mask], h_vals[fit_mask]
                    )
                elif (a_vals > 0).sum() > 100:
                    # Fallback to all data if not enough in lower range
                    valid = (a_vals > 0) & (h_vals > 0)
                    slope, intercept, _, _, _ = stats.linregress(
                        a_vals[valid], h_vals[valid]
                    )
                else:
                    slope = 1.0
                    intercept = 0
                
                # Calculate deviation from diagonal for ALL events
                expected_h = slope * data[area] + intercept
                deviation = np.abs(data[height] - expected_h) / (expected_h + 1)
                
                gate_mask = parent_mask & (deviation < slope_tolerance)
                
                gate_params = {
                    'area': area, 'height': height, 'method': method,
                    'slope': slope, 'intercept': intercept, 'tolerance': slope_tolerance,
                    'fit_percentile': fit_percentile, 'fit_threshold': area_threshold
                }
                
            elif method == 'width' and width:
                # Gate based on width distribution
                w_lo = np.percentile(parent_data[width], width_percentile_low)
                w_hi = np.percentile(parent_data[width], width_percentile_high)
                
                gate_mask = parent_mask & (data[width] >= w_lo) & (data[width] <= w_hi)
                
                gate_params = {
                    'area': area, 'height': height, 'width': width,
                    'method': method, 'width_min': w_lo, 'width_max': w_hi
                }
            else:
                raise ValueError(f"Invalid method '{method}' or missing width channel")
            
            gated_count = gate_mask.sum()
            
            fa.gates[name] = GateResult(
                name=name,
                parent_count=parent_count,
                gated_count=gated_count,
                percent_of_parent=100 * gated_count / parent_count if parent_count > 0 else 0,
                percent_of_total=100 * gated_count / fa.total_events,
                mask=gate_mask,
                gate_params=gate_params
            )
            fa.current_mask = gate_mask
        
        self._generate_gate_figures(name, 'singlet', x_channel=area, y_channel=height)
        return self
    
    def gate_fluorophore(
        self,
        channel: str,
        threshold: float = None,
        method: str = 'manual',
        percentile: float = 95,
        positive_name: str = None,
        negative_name: str = None,
        log_scale: bool = True
    ) -> 'FCAnalyzer':
        """
        Gate positive/negative populations based on fluorophore intensity.
        
        Args:
            channel: Fluorescence channel name
            threshold: Manual threshold value
            method: 'manual' for fixed threshold, 'percentile' for automatic
            percentile: For percentile method, events above this are positive
            positive_name: Name for positive gate (default: '{channel}+')
            negative_name: Name for negative gate (default: '{channel}-')
            log_scale: Whether to use log scale for visualization
        """
        self._validate_channels([channel])
        
        pos_name = positive_name or f'{channel}+'
        neg_name = negative_name or f'{channel}-'
        
        self.gate_sequence.append(pos_name)
        self.gate_sequence.append(neg_name)
        
        for filename, fa in self.files.items():
            parent_mask = fa.current_mask.copy()
            parent_count = parent_mask.sum()
            data = fa.raw_data
            parent_data = data.loc[parent_mask]
            
            if method == 'percentile':
                if parent_count == 0:
                    # No parent events - use a default threshold
                    thresh = 0
                else:
                    thresh = np.percentile(parent_data[channel], percentile)
            else:
                thresh = threshold
            
            if thresh is None:
                raise ValueError("Must provide threshold or use method='percentile'")
            
            pos_mask = parent_mask & (data[channel] >= thresh)
            neg_mask = parent_mask & (data[channel] < thresh)
            
            pos_count = pos_mask.sum()
            neg_count = neg_mask.sum()
            
            fa.gates[pos_name] = GateResult(
                name=pos_name,
                parent_count=parent_count,
                gated_count=pos_count,
                percent_of_parent=100 * pos_count / parent_count if parent_count > 0 else 0,
                percent_of_total=100 * pos_count / fa.total_events,
                mask=pos_mask,
                gate_params={
                    'channel': channel, 'threshold': thresh, 
                    'direction': 'positive', 'log_scale': log_scale
                }
            )
            
            fa.gates[neg_name] = GateResult(
                name=neg_name,
                parent_count=parent_count,
                gated_count=neg_count,
                percent_of_parent=100 * neg_count / parent_count if parent_count > 0 else 0,
                percent_of_total=100 * neg_count / fa.total_events,
                mask=neg_mask,
                gate_params={
                    'channel': channel, 'threshold': thresh,
                    'direction': 'negative', 'log_scale': log_scale
                }
            )
        
        # Generate histogram figures
        self._generate_gate_figures(pos_name, 'histogram', x_channel=channel, 
                                   threshold=threshold or thresh, log_scale=log_scale)
        return self
    
    def combine_gates(
        self,
        gates: list[str],
        operation: str = 'AND',
        name: str = None
    ) -> 'FCAnalyzer':
        """
        Combine existing gates with AND/OR logic.
        
        Args:
            gates: List of gate names to combine
            operation: 'AND' (intersection) or 'OR' (union)
            name: Name for new combined gate (default: auto-generated)
            
        Returns:
            self for method chaining
            
        Examples:
            # Double positive (must pass both gates)
            analyzer.combine_gates(['GFP+', 'mCherry+'], operation='AND', name='DoublePositive')
            
            # Either positive (passes at least one gate)
            analyzer.combine_gates(['GFP+', 'mCherry+'], operation='OR', name='EitherPositive')
            
            # Triple positive
            analyzer.combine_gates(['GFP+', 'mCherry+', 'BFP+'], operation='AND', name='TriplePos')
        """
        if len(gates) < 2:
            raise ValueError("Must provide at least 2 gates to combine")
        
        operation = operation.upper()
        if operation not in ('AND', 'OR'):
            raise ValueError("operation must be 'AND' or 'OR'")
        
        # Auto-generate name if not provided
        if name is None:
            joiner = '&' if operation == 'AND' else '|'
            name = f"({joiner.join(gates)})"
        
        self.gate_sequence.append(name)
        
        for filename, fa in self.files.items():
            # Check all gates exist for this file
            missing = [g for g in gates if g not in fa.gates]
            if missing:
                warnings.warn(f"{filename}: Missing gates {missing}, skipping combine")
                continue
            
            # Get the parent (use first gate's parent as reference)
            first_gate = fa.gates[gates[0]]
            parent_count = first_gate.parent_count
            
            # Combine masks
            if operation == 'AND':
                combined_mask = fa.gates[gates[0]].mask.copy()
                for gate_name in gates[1:]:
                    combined_mask = combined_mask & fa.gates[gate_name].mask
            else:  # OR
                combined_mask = fa.gates[gates[0]].mask.copy()
                for gate_name in gates[1:]:
                    combined_mask = combined_mask | fa.gates[gate_name].mask
            
            gated_count = combined_mask.sum()
            
            fa.gates[name] = GateResult(
                name=name,
                parent_count=parent_count,
                gated_count=gated_count,
                percent_of_parent=100 * gated_count / parent_count if parent_count > 0 else 0,
                percent_of_total=100 * gated_count / fa.total_events,
                mask=combined_mask,
                gate_params={
                    'combined_gates': gates,
                    'operation': operation
                }
            )
        
        print(f"Created combined gate '{name}' from {gates} ({operation})")
        return self
    
    def _get_combined_mask(
        self, 
        fa: 'FileAnalysis', 
        gates: Union[str, list[str]], 
        operation: str = 'AND'
    ) -> np.ndarray:
        """
        Get a mask combining multiple gates (internal helper).
        
        Args:
            fa: FileAnalysis object
            gates: Single gate name or list of gate names
            operation: 'AND' or 'OR' for combining multiple gates
            
        Returns:
            Combined boolean mask
        """
        if isinstance(gates, str):
            # Single gate
            if gates in fa.gates:
                return fa.gates[gates].mask
            else:
                return fa.current_mask
        
        # Multiple gates
        if len(gates) == 0:
            return fa.current_mask
        
        masks = []
        for gate_name in gates:
            if gate_name in fa.gates:
                masks.append(fa.gates[gate_name].mask)
        
        if len(masks) == 0:
            return fa.current_mask
        
        if operation.upper() == 'AND':
            result = masks[0].copy()
            for m in masks[1:]:
                result = result & m
        else:  # OR
            result = masks[0].copy()
            for m in masks[1:]:
                result = result | m
        
        return result
    
    def _generate_gate_figures(self, *args, **kwargs) -> None:
        """Gate QC figures (per-well gate plots) are not part of the paper's processing and were
        removed with the plotting code; the gate methods still call this hook."""
        return None
    
    def get_channel_statistics(
        self, 
        channels: list[str] = None, 
        gate: Union[str, list[str]] = None,
        gate_operation: str = 'AND',
        include_metadata: list[str] = None,
        include_counts: bool = False,
        stats: list[str] = ['mean', 'median', 'std', 'cv']
    ) -> pd.DataFrame:
        """
        Get intensity statistics for specified channels.
        
        Args:
            channels: List of channels to analyze (default: all fluorescence channels)
            gate: Gate name OR list of gate names to combine
            gate_operation: 'AND' (intersection) or 'OR' (union) when multiple gates provided
            include_metadata: Metadata columns to include (None=none, 'all' or list of names)
            include_counts: If True, include count columns for all gates
            stats: Which statistics to compute: 'mean', 'median', 'std', 'cv', 'gmean', 'min', 'max'
            
        Example:
            # Basic usage
            df = analyzer.get_channel_statistics(['FITC-A', 'PE-A'], gate='Singlets')
            
            # With metadata and counts
            df = analyzer.get_channel_statistics(
                channels=['Comp-FITC-A', 'Comp-PE-A'],
                gate='Singlets',
                include_metadata=['condition', 'replicate', 'timepoint'],
                include_counts=True,
                stats=['median']  # just medians
            )
            
            # Double positive cells (AND logic)
            df = analyzer.get_channel_statistics(
                channels=['Comp-FITC-A'],
                gate=['GFP+', 'mCherry+'],  # must be in BOTH gates
                gate_operation='AND'
            )
            
            # Either positive (OR logic)
            df = analyzer.get_channel_statistics(
                channels=['Comp-FITC-A'],
                gate=['GFP+', 'mCherry+'],
                gate_operation='OR'
            )
        """
        rows = []
        
        for filename, fa in self.files.items():
            # Get mask - handles single gate, multiple gates, or None
            mask = self._get_combined_mask(fa, gate, gate_operation)
            
            data = fa.raw_data.loc[mask]
            chans = channels or [c for c in data.columns if 'FSC' not in c and 'SSC' not in c]
            
            row = {'filename': filename}
            
            # Add population count for the selected gate(s)
            row['population_count'] = mask.sum()
            
            # Add metadata if requested
            if include_metadata == 'all':
                for key, val in fa.metadata.items():
                    if not key.startswith('_'):
                        row[key] = val
            elif include_metadata:
                for key in include_metadata:
                    row[key] = fa.metadata.get(key, None)
            
            # Add gate counts if requested
            if include_counts:
                row['total_events'] = fa.total_events
                for gate_name in self.gate_sequence:
                    if gate_name in fa.gates:
                        row[f'{gate_name}_count'] = fa.gates[gate_name].gated_count
                        row[f'{gate_name}_pct'] = fa.gates[gate_name].percent_of_parent
            
            # Add channel statistics
            for ch in chans:
                if ch in data.columns:
                    vals = data[ch]
                    if 'mean' in stats:
                        row[f'{ch}_mean'] = vals.mean()
                    if 'median' in stats:
                        row[f'{ch}_median'] = vals.median()
                    if 'std' in stats:
                        row[f'{ch}_std'] = vals.std()
                    if 'cv' in stats:
                        row[f'{ch}_cv'] = vals.std() / vals.mean() if vals.mean() != 0 else 0
                    if 'gmean' in stats:
                        positive = vals[vals > 0]
                        row[f'{ch}_gmean'] = np.exp(np.log(positive).mean()) if len(positive) > 0 else np.nan
                    if 'min' in stats:
                        row[f'{ch}_min'] = vals.min()
                    if 'max' in stats:
                        row[f'{ch}_max'] = vals.max()
            
            rows.append(row)
        
        return pd.DataFrame(rows)
