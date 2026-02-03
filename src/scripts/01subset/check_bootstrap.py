# or run in terminal with
# python -c " python codes below "
import h5py, numpy as np

with h5py.File('results/gradient_indices.h5', 'r') as h5:
    print('=== datasets ===')
    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(f'  {name:30s} shape={str(obj.shape):20s} dtype={obj.dtype}  chunks={obj.chunks}')


    h5.visititems(visit)

    # all indices in valid sample range
    max_val = 20
    for ds in ['indices/low', 'indices/high', 'indices/low_boot', 'indices/high_boot']:
        arr = h5[ds][:]
        assert arr.min() >= 0 and arr.max() < max_val, f'{ds} out of range'


    print()
    print('All indices in valid range [0, 20)')

    # per-gene: boot indices must be subset of the partition
    low_all      = h5['indices/low'][:]       # (n_genes, k_low)
    low_boot_all = h5['indices/low_boot'][:]  # (n_genes, n_boot, k_resample)
    for g in range(5):
        allowed = set(low_all[g].tolist())
        used    = set(low_boot_all[g].flatten().tolist())
        assert used <= allowed, f'gene {g} low_boot outside partition'
    print('All low_boot indices subset of gene partition')

    high_all      = h5['indices/high'][:]
    high_boot_all = h5['indices/high_boot'][:]
    for g in range(5):
        allowed = set(high_all[g].tolist())
        used    = set(high_boot_all[g].flatten().tolist())
        assert used <= allowed, f'gene {g} high_boot outside partition'
    print('All high_boot indices subset of gene partition')
