import flowkit as fk

def inspect_fcs_file(filepath, verbose=False):
    """Load an FCS file. Returns the FlowKit sample and its raw data as a DataFrame."""
    sample = fk.Sample(filepath)
    df_sample = sample.as_dataframe(source='raw', col_multi_index=False)
    if verbose:
        metadata = sample.metadata
        print(f"Total events: {sample.event_count}")
        print(f"PNN labels: {sample.pnn_labels}")
        print(f"PNS labels: {sample.pns_labels}")
        print(f"Instrument (cyt): {metadata.get('cyt')}")
        print(f"Date: {metadata.get('date')}")
        print(f"System (sys): {metadata.get('sys')}")
    return sample, df_sample


def check_spillover_matrix(sample):
    """Return any metadata keys that look like a compensation (spillover) matrix."""
    spill_keys = [k for k in sample.metadata if "spill" in k.lower() or "comp" in k.lower()]
    if spill_keys:
        print("Spillover matrix found")
    else:
        print("No spillover matrix found. Compensation not applied.")
    return spill_keys


def asinh_transformation(sample):
    """Arcsinh transform for linear-scale fluorescence (used for the Attune file in Week 1).
    param_t is the top of that file's measurement range. Not used for the HSCT files,
    which are already log-scaled (0-1023)."""
    asinh_xform = fk.transforms.AsinhTransform(param_t=16777216.0, param_m=4.0, param_a=0.0)
    sample.apply_transform(asinh_xform)
    return sample.as_dataframe(source='xform', col_multi_index=False)