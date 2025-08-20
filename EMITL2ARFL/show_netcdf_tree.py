import netCDF4

def show_netcdf_tree(filename, indent=0, group=None):
    """
    Recursively returns the structure of a NetCDF file as a tree string.
    """
    lines = []
    if group is None:
        ds = netCDF4.Dataset(filename, 'r')
        group = ds
        lines.append(f"File: {filename}")
        close_ds = True
    else:
        ds = group
        close_ds = False

    prefix = "  " * indent
    # Dimensions
    if ds.dimensions:
        lines.append(f"{prefix}Dimensions:")
        for dim_name, dim in ds.dimensions.items():
            lines.append(f"{prefix}  {dim_name}: size={len(dim)}")

    # Variables
    if ds.variables:
        lines.append(f"{prefix}Variables:")
        for var_name, var in ds.variables.items():
            lines.append(f"{prefix}  {var_name}: shape={var.shape}, dtype={var.dtype}")

    # Attributes
    if ds.ncattrs():
        lines.append(f"{prefix}Attributes:")
        for attr in ds.ncattrs():
            lines.append(f"{prefix}  {attr}: {getattr(ds, attr)}")

    # Groups
    if hasattr(ds, 'groups') and ds.groups:
        lines.append(f"{prefix}Groups:")
        for group_name, subgroup in ds.groups.items():
            lines.append(f"{prefix}  Group: {group_name}")
            lines.append(show_netcdf_tree(filename, indent + 2, subgroup))

    if close_ds:
        ds.close()
    return "\n".join(lines)
