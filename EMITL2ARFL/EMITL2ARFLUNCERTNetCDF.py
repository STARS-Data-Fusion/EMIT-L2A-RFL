from .EMITNetCDF import EMITNetCDF

class EMITL2ARFLUNCERTNetCDF(EMITNetCDF):
    """
    This class encapsulates access to the EMIT L2A Reflectance Uncertainty NetCDF file.
    """
    def __repr__(self) -> str:
        return f"EMITL2ARFLUNCERTNetCDF(\"{self.filename}\")"
