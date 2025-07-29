from .EMITNetCDF import EMITNetCDF

class EMITL2AMASKNetCDF(EMITNetCDF):
    """
    This class encapsulates access to the EMIT L2A Reflectance Mask NetCDF file.
    """
    def __repr__(self) -> str:
        return f"EMITL2AMASKNetCDF(\"{self.filename}\")"
