# HDF5 Attribute Size Limit Fix

## Problem

When processing large gene expression datasets, the pipeline encountered this error:

```
OSError: Unable to synchronously create attribute (object header message is too large)
```

### Root Cause

HDF5 attributes have a strict **64 KB size limit** per attribute. In the original implementation, gene and sample names were stored as attributes:

```python
ds.attrs["gene_names"] = df.index.astype(str).tolist()  # ❌ Can exceed 64 KB
ds.attrs["sample_names"] = df.columns.astype(str).tolist()
```

For datasets with many genes (e.g., 20,000+ genes with long gene IDs like "FBgn0123456"), the serialized list of gene names can easily exceed 64 KB.

## Solution

**Store gene and sample names as separate HDF5 datasets instead of attributes.**

Datasets in HDF5 have **no size limit** and are designed for large arrays, while attributes are meant for small metadata.

### Updated HDF5 Schema

**Before (problematic):**
```
expression.h5
  └── /expr
      └── .attrs
          ├── n_genes
          ├── n_samples
          ├── gene_names      ❌ Size limited to 64 KB
          └── sample_names    ❌ Size limited to 64 KB
```

**After (fixed):**
```
expression.h5
  ├── /expr
  │   └── .attrs
  │       ├── n_genes         ✅ Scalar metadata
  │       └── n_samples       ✅ Scalar metadata
  ├── /gene_names             ✅ Dataset (unlimited size)
  └── /sample_names           ✅ Dataset (unlimited size)
```

## Changes Made

### Modified File: `src/scripts/00preprocess/00convert_expr_to_hdf5.py`

#### 1. Function `convert_tsv_to_hdf5()` (lines 76-81)

**Before:**
```python
# Store metadata as attributes
ds.attrs["n_genes"] = expr.shape[0]
ds.attrs["n_samples"] = expr.shape[1]
ds.attrs["gene_names"] = df.index.astype(str).tolist()      # ❌ REMOVED
ds.attrs["sample_names"] = df.columns.astype(str).tolist()  # ❌ REMOVED

# Optional: store gene/sample names as separate datasets
f.create_dataset("gene_names", data=df.index.astype(str).tolist(), dtype=h5py.string_dtype())
f.create_dataset("sample_names", data=df.columns.astype(str).tolist(), dtype=h5py.string_dtype())
```

**After:**
```python
# Store metadata as attributes
ds.attrs["n_genes"] = expr.shape[0]
ds.attrs["n_samples"] = expr.shape[1]

# Store gene/sample names as separate datasets (avoids 64 KB attribute size limit)
f.create_dataset("gene_names", data=df.index.astype(str).tolist(), dtype=h5py.string_dtype())
f.create_dataset("sample_names", data=df.columns.astype(str).tolist(), dtype=h5py.string_dtype())
```

#### 2. Function `generate_toy_data()` (lines 186-191)

**Before:**
```python
with h5py.File(h5_path, "w") as f:
    ds = f.create_dataset("expr", data=expr, dtype=np.float32)
    ds.attrs["n_genes"] = n_genes
    ds.attrs["n_samples"] = n_samples
    ds.attrs["gene_names"] = gene_names        # ❌ REMOVED
    ds.attrs["sample_names"] = sample_names    # ❌ REMOVED
    f.create_dataset("gene_names", data=gene_names, dtype=h5py.string_dtype())
    f.create_dataset("sample_names", data=sample_names, dtype=h5py.string_dtype())
```

**After:**
```python
with h5py.File(h5_path, "w") as f:
    ds = f.create_dataset("expr", data=expr, dtype=np.float32)
    ds.attrs["n_genes"] = n_genes
    ds.attrs["n_samples"] = n_samples
    # Store names as datasets (avoids 64 KB attribute size limit)
    f.create_dataset("gene_names", data=gene_names, dtype=h5py.string_dtype())
    f.create_dataset("sample_names", data=sample_names, dtype=h5py.string_dtype())
```

#### 3. Updated Documentation (lines 15-20)

Updated the output schema documentation in the module docstring to reflect the new structure.

## Downstream Compatibility

✅ **No changes required** for downstream scripts! All scripts were already reading from datasets:

| Script | Line | Code |
|--------|------|------|
| `02a_calc_base_correlations.py` | 368 | `gene_names = [...] f["gene_names"][:]` |
| `02b_bootstrap_significant_edges.py` | 165, 389 | `gene_names = [...] h5["gene_names"][:]` |
| `03_reconstruct_diff_network.py` | 734, 1064 | `gene_names = [...] h5["gene_names"][:]` |
| `04_collect_focus_gene_topology.py` | 1064 | `all_gene_names = [...] h5["gene_names"][:]` |
| `05_prepare_visualization_data.py` | 62, 102 | `gene_names = [...] h5["gene_names"][:]` |

The downstream scripts use the pattern:
```python
if "gene_names" in h5:
    gene_names = [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
```

This works correctly with datasets and includes proper handling of byte/string decoding.

## Benefits

1. **No Size Limit**: HDF5 datasets can store arbitrarily large arrays of strings
2. **Better Performance**: Datasets are optimized for array data with chunking and compression
3. **Backwards Compatible**: All existing downstream scripts already use the correct access pattern
4. **Cleaner Design**: Attributes for scalar metadata, datasets for array data

## Technical Details

### Why Attributes Have Size Limits

HDF5 stores attributes in the object header, which has a maximum size. This design is optimized for small metadata (numbers, short strings, small arrays). When attributes grow too large, they fragment the header and cause the error.

### Why Datasets Don't Have Size Limits

Datasets store data separately from the object header and use a chunked storage layout. They can grow to petabyte scale without issues.

### String Encoding

The fix uses `h5py.string_dtype()` which creates variable-length UTF-8 strings:
- Efficient storage (no padding)
- Proper Unicode support
- Automatic handling by h5py

### Reading Strings

When reading, strings may be returned as bytes (Python 2 compatibility). The downstream scripts handle this:
```python
gene_names = [x.decode() if isinstance(x, bytes) else x for x in h5["gene_names"][:]]
```

## Testing

After applying this fix, you can process large datasets without hitting the 64 KB limit:

```bash
# Example: 20,000 genes with long gene IDs
python src/scripts/00preprocess/00convert_expr_to_hdf5.py \
    --expr-tsv data/large_dataset.txt \
    --out-h5 data/expression.h5

# Verify the structure
python -c "
import h5py
with h5py.File('data/expression.h5', 'r') as f:
    print('Datasets:', list(f.keys()))
    print('Expression attrs:', list(f['expr'].attrs.keys()))
    print('Genes:', len(f['gene_names']))
    print('First gene:', f['gene_names'][0])
"
```

Expected output:
```
Datasets: ['expr', 'gene_names', 'sample_names']
Expression attrs: ['n_genes', 'n_samples']
Genes: 20000
First gene: FBgn0000001
```

## Migration

If you have existing HDF5 files with gene names stored as attributes, you can migrate them:

```python
import h5py

def migrate_h5_file(input_path, output_path):
    """Convert old format (attrs) to new format (datasets)."""
    with h5py.File(input_path, 'r') as f_in, h5py.File(output_path, 'w') as f_out:
        # Copy expression data
        f_out.create_dataset('expr', data=f_in['expr'][:],
                            compression='gzip', compression_opts=4)

        # Copy scalar attributes
        f_out['expr'].attrs['n_genes'] = f_in['expr'].attrs['n_genes']
        f_out['expr'].attrs['n_samples'] = f_in['expr'].attrs['n_samples']

        # Migrate gene_names from attribute to dataset
        if 'gene_names' in f_in['expr'].attrs:
            gene_names = f_in['expr'].attrs['gene_names']
            f_out.create_dataset('gene_names', data=gene_names,
                                dtype=h5py.string_dtype())
        elif 'gene_names' in f_in:
            # Already in correct format
            f_out.create_dataset('gene_names', data=f_in['gene_names'][:],
                                dtype=h5py.string_dtype())

        # Same for sample_names
        if 'sample_names' in f_in['expr'].attrs:
            sample_names = f_in['expr'].attrs['sample_names']
            f_out.create_dataset('sample_names', data=sample_names,
                                dtype=h5py.string_dtype())
        elif 'sample_names' in f_in:
            f_out.create_dataset('sample_names', data=f_in['sample_names'][:],
                                dtype=h5py.string_dtype())

# Usage
migrate_h5_file('data/old_expression.h5', 'data/expression.h5')
```

## Related Issues

- HDF5 attribute size limit: https://github.com/h5py/h5py/issues/1269
- HDF5 documentation on attributes vs datasets: https://docs.h5py.org/en/stable/high/attr.html

## Status

✅ **Fixed and tested** - All scripts updated and verified compatible.

**Last updated:** 2026-02-13
