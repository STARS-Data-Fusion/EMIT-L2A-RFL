from .EMITNetCDF import EMITNetCDF

class EMITL2ARFLNetCDF(EMITNetCDF):
    """
    This class encapsulates access to the EMIT L2A Reflectance Geolocation NetCDF file.
    """
    def __repr__(self) -> str:
        return f"EMITL2ARFLNetCDF(\"{self.filename}\")"
