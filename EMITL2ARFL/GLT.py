
from rasters import Raster, RasterGrid
import numpy as np

class GeometryLookupTable(np.ndarray):
    """
    Geometry Lookup Table (GLT):
    A gridded, georeferenced array storing 1-based row and column indices that map each pixel in the grid
    to its corresponding location in a swath or geolocated array. This enables geospatial referencing and reprojection.

    This class subclasses numpy.ndarray, attaches geolocation metadata via a RasterGrid, and provides convenient
    access to row/col indices as Raster objects with geometry.
    """
    def __new__(cls, GLT_array: np.ndarray, geometry: RasterGrid = None):
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
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.geometry = getattr(obj, 'geometry', None)

    @property
    def rows(self) -> Raster:
        """
        Returns the 1-based row indices (glt_y) as a Raster object with geometry.
        Each value is a row index into the swath/geolocated array.
        """
        return Raster(self[..., 0], geometry=self.geometry)

    @property
    def cols(self) -> Raster:
        """
        Returns the 1-based column indices (glt_x) as a Raster object with geometry.
        Each value is a column index into the swath/geolocated array.
        """
        return Raster(self[..., 1], geometry=self.geometry)
    
    
