
from rasters import Raster, RasterGeolocation
import numpy as np

class GLT(np.ndarray):
    """
    Numpy ndarray subclass for the GLT array, with geolocation attribute and rows/cols properties.
    """
    def __new__(cls, GLT_array: np.ndarray, geolocation: RasterGeolocation = None):
        # Cast input array to GLT subclass
        if not isinstance(geolocation, RasterGeolocation):
            raise TypeError(f"geolocation must be a RasterGeolocation, got {type(geolocation)}")
        obj = np.asarray(GLT_array).view(cls)
        obj.geolocation = geolocation
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.geolocation = getattr(obj, 'geolocation', None)

    @property
    def rows(self) -> Raster:
        """Return the row indices (glt_y) as a Raster object with geometry from geolocation."""
        return Raster(self[..., 0], geometry=self.geolocation)

    @property
    def cols(self) -> Raster:
        """Return the column indices (glt_x) as a Raster object with geometry from geolocation."""
        return Raster(self[..., 1], geometry=self.geolocation)
