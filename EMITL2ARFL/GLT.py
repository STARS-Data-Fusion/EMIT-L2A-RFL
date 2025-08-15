
from rasters import Raster, RasterGrid
import numpy as np
from rasterio.windows import Window

from .constants import GLT_NODATA_VALUE

class GeometryLookupTable(np.ndarray):
    """
    Geometry Lookup Table (GLT):
    A gridded, georeferenced array storing 1-based row and column indices that map each pixel in the grid
    to its corresponding location in a swath or geolocated array. This enables geospatial referencing and reprojection.

    This class subclasses numpy.ndarray, attaches geolocation metadata via a RasterGrid, and provides convenient
    access to row/col indices as Raster objects with geometry.
    """
    def __new__(
            cls, 
            GLT_array: np.ndarray, 
            GLT_nodata_value: int = GLT_NODATA_VALUE,
            geometry: RasterGrid = None):
        # Cast input array to GLT subclass
        """
        Create a GeometryLookupTable from a GLT array and geolocation metadata.
        Parameters:
            GLT_array (np.ndarray): A (rows, cols, 2) ndarray of 1-based indices.
            geometry (RasterGrid): Georeferencing metadata describing the grid.
        Returns:
            GeometryLookupTable: GLT array with attached geometry.
        """
        if not isinstance(geometry, RasterGrid):
            raise TypeError(f"geometry must be a RasterGrid, got {type(geometry)}")
        obj = np.asarray(GLT_array).view(cls)
        obj.geometry = geometry
        obj.GLT_nodata_value = GLT_nodata_value
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.geometry = getattr(obj, 'geometry', None)
        if hasattr(obj, 'GLT_nodata_value'):
            self.GLT_nodata_value = obj.GLT_nodata_value
        else:
            self.GLT_nodata_value = GLT_NODATA_VALUE

    @property
    def rows(self) -> Raster:
        """
        Returns the zero-based row indices (glt_y) as a Raster object with geometry.
        Each value is a row index into the swath/geolocated array.
        """
        return Raster(self[..., 0] - 1, geometry=self.geometry)

    @property
    def min_row(self) -> int:
        """
        Returns the minimum row index (glt_y) for the GLT, ignoring nodata values.
        """
        valid = self[..., 0][self[..., 0] != self.GLT_nodata_value]
        if valid.size == 0:
            raise ValueError(f"No valid row indices found in GLT (all are nodata value {self.GLT_nodata_value})")
        return int(np.nanmin(valid)) - 1

    @property
    def max_row(self) -> int:
        """
        Returns the maximum row index (glt_y) for the GLT, ignoring nodata values.
        """
        valid = self[..., 0][self[..., 0] != self.GLT_nodata_value]
        if valid.size == 0:
            raise ValueError(f"No valid row indices found in GLT (all are nodata value {self.GLT_nodata_value})")
        return int(np.nanmax(valid)) - 1

    @property
    def cols(self) -> Raster:
        """
        Returns the zero-based column indices (glt_x) as a Raster object with geometry.
        Each value is a column index into the swath/geolocated array.
        """
        return Raster(self[..., 1] - 1, geometry=self.geometry)

    @property
    def min_col(self) -> int:
        """
        Returns the minimum column index (glt_x) for the GLT, ignoring nodata values.
        """
        valid = self[..., 1][self[..., 1] != self.GLT_nodata_value]
        if valid.size == 0:
            raise ValueError(f"No valid column indices found in GLT (all are nodata value {self.GLT_nodata_value})")
        return int(np.nanmin(valid)) - 1

    @property
    def max_col(self) -> int:
        """
        Returns the maximum column index (glt_x) for the GLT, ignoring nodata values.
        """
        valid = self[..., 1][self[..., 1] != self.GLT_nodata_value]
        if valid.size == 0:
            raise ValueError(f"No valid column indices found in GLT (all are nodata value {self.GLT_nodata_value})")
        return int(np.nanmax(valid)) - 1

    @property
    def swath_window(self) -> Window:
        """
        Returns the swath window for the GLT.
        """
        return Window(
            col_off=self.min_col,
            row_off=self.min_row,
            width=self.max_col - self.min_col + 1,
            height=self.max_row - self.min_row + 1
        )