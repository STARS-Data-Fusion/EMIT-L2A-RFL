
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

    def __repr__(self):
        return f"GeometryLookupTable(shape={self.shape}, min_row={self.min_row}, max_row={self.max_row}, min_col={self.min_col}, max_col={self.max_col})"

    def __str__(self):
        return self.__repr__()

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
    
    def adjust_indices(self, window: Window) -> "GeometryLookupTable":
        """
        Adjusts the GLT indices based on a provided window.
        This is useful when subsetting the GLT to a specific region.
        
        Parameters:
            window (Window): The window to adjust indices by.
        
        Returns:
            GeometryLookupTable: A new GLT with adjusted indices.
        """
        if window is None:
            return self

        adjusted_array = self.copy()
        mask = (adjusted_array[..., 0] != self.GLT_nodata_value) & (adjusted_array[..., 1] != self.GLT_nodata_value)
        adjusted_array[..., 0][mask] -= int(window.row_off)
        adjusted_array[..., 1][mask] -= int(window.col_off)

        # Calculate bounds for the swath window
        swath_height = int(window.height)
        swath_width = int(window.width)
        max_row = int(np.nanmax(adjusted_array[..., 0][mask])) - 1 if np.any(mask) else -1
        max_col = int(np.nanmax(adjusted_array[..., 1][mask])) - 1 if np.any(mask) else -1

        if max_row >= swath_height or max_col >= swath_width:
            raise ValueError(
                f"Adjusted GLT indices are out of bounds for the swath window: "
                f"max_row={max_row}, max_col={max_col}, "
                f"swath_height={swath_height}, swath_width={swath_width}, "
                f"window(row_off={window.row_off}, col_off={window.col_off}, height={window.height}, width={window.width})"
            )

        return GeometryLookupTable(GLT_array=adjusted_array, geometry=self.geometry, GLT_nodata_value=self.GLT_nodata_value)
    