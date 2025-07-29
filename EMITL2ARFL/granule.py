from glob import glob
from os.path import join, abspath, dirname, expanduser
from typing import List

import numpy as np

from rasters import Raster, RasterGeometry, RasterGeolocation

from .read_netcdf_raster import read_netcdf_raster
from .read_netcdf_array import read_netcdf_array
from .read_latitude_array import read_latitude_array
from .read_longitude_array import read_longitude_array
from .read_geolocation import read_geolocation

from .constants import *
from .emit_ortho_raster import emit_ortho_raster
from .quality_mask import quality_mask
from .extract_GLT_array import extract_GLT_array
from .extract_GLT import extract_GLT
from .GLT import GLT
from .show_netcdf_tree import show_netcdf_tree

class EMITNetCDF:
    def __init__(self, filename: str) -> None:
        """
        Initialize an EMITNetCDF object for a given NetCDF file.

        Args:
            filename (str): Path to the NetCDF file.
        """
        self.filename: str = abspath(expanduser(filename))

    def __repr__(self) -> str:
        """
        Return a string representation of the EMITNetCDF object.
        """
        return f"EMITNetCDF(filename=\"{self.filename}\")"

    def show_tree(self, indent: int = 0) -> str:
        """
        Show the NetCDF file structure as a tree.

        Args:
            indent (int, optional): Indentation level for the tree. Defaults to 0.

        Returns:
            str: String representation of the NetCDF file tree.
        """
        return show_netcdf_tree(self.filename, indent)
    
    tree = property(show_tree)

    @property
    def lat(self) -> np.ndarray:
        """
        Read the latitude array from the NetCDF file.

        Returns:
            np.ndarray: Latitude values.
        """
        return read_latitude_array(self.filename)

    @property
    def lon(self) -> np.ndarray:
        """
        Read the longitude array from the NetCDF file.

        Returns:
            np.ndarray: Longitude values.
        """
        return read_longitude_array(self.filename)

    @property
    def geolocation(self) -> RasterGeolocation:
        """
        Read the geolocation information from the NetCDF file.

        Returns:
            RasterGeolocation: Geolocation object containing spatial info.
        """
        return read_geolocation(self.filename)
    
    def extract_GLT(self) -> GLT:
        """
        Extract the Geometry Lookup Table (GLT) from the NetCDF file.

        Returns:
            GLT: Geometry Lookup Table object.
        """
        return extract_GLT(self.filename)

    GLT = property(extract_GLT)

    def read(self, group: str, variable: str):
        """
        Read a variable as a Raster object from a specified group in the NetCDF file.

        Args:
            group (str): Name of the group in the NetCDF file.
            variable (str): Name of the variable to read.

        Returns:
            Raster: Raster object of the requested variable with geolocation.
        """
        return read_netcdf_raster(
            filename=self.filename,
            group=group,
            variable=variable
        )

    def read_elevation(self) -> Raster:
        """
        Read the elevation raster from the NetCDF file.

        Returns:
            Raster: Elevation raster object.
        """
        return read_netcdf_raster(
            filename=self.filename,
            group="location",
            variable="elev"
        )
    
    elevation = property(read_elevation)

    def read_array(self, group: str, variable: str) -> np.ndarray:
        """
        Read a variable array from a specified group in the NetCDF file.

        Args:
            group (str): Name of the group in the NetCDF file.
            variable (str): Name of the variable to read.

        Returns:
            np.ndarray: Array of the requested variable.
        """
        return read_netcdf_array(self.filename, group, variable)

class EMITL2ARFLNetCDF(EMITNetCDF):
    """
    This class encapsulates access to the EMIT L2A Reflectance Geolocation NetCDF file.

    Example EMIT L2A Reflectance Geolocation NetCDF file structure:

    File: /Users/halverso/data/EMIT_L2A_RFL/EMIT_L2A_RFL_001_20241006T165200_2428011_004/EMIT_L2A_RFL_001_20241006T165200_2428011_004.nc
    Dimensions:
    downtrack: size=1280
    crosstrack: size=1242
    bands: size=285
    ortho_y: size=1886
    ortho_x: size=2298
    Variables:
    reflectance: shape=(1280, 1242, 285), dtype=float32
    Attributes:
    ncei_template_version: NCEI_NetCDF_Swath_Template_v2.0
    summary: The Earth Surface Mineral Dust Source Investigation (EMIT) is an Earth Ventures-Instrument (EVI-4) Mission that maps the surface mineralogy of arid dust source regions via imaging spectroscopy in the visible and short-wave infrared (VSWIR). Installed on the International Space Station (ISS), the EMIT instrument is a Dyson imaging spectrometer that uses contiguous spectroscopic measurements from 410 to 2450 nm to resolve absoprtion features of iron oxides, clays, sulfates, carbonates, and other dust-forming minerals. During its one-year mission, EMIT will observe the sunlit Earth's dust source regions that occur within +/-52° latitude and produce maps of the source regions that can be used to improve forecasts of the role of mineral dust in the radiative forcing (warming or cooling) of the atmosphere.\n\nThis file contains L2A estimated surface reflectances and geolocation data. Reflectance estimates are created using an Optimal Estimation technique - see ATBD for details. Reflectance values are reported as fractions (relative to 1). Geolocation data (latitude, longitude, height) and a lookup table to project the data are also included.
    keywords: Imaging Spectroscopy, minerals, EMIT, dust, radiative forcing
    Conventions: CF-1.63
    sensor: EMIT (Earth Surface Mineral Dust Source Investigation)
    instrument: EMIT
    platform: ISS
    institution: NASA Jet Propulsion Laboratory/California Institute of Technology
    license: https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/
    naming_authority: LPDAAC
    date_created: 2024-10-09T02:59:07Z
    keywords_vocabulary: NASA Global Change Master Directory (GCMD) Science Keywords
    stdname_vocabulary: NetCDF Climate and Forecast (CF) Metadata Convention
    creator_name: Jet Propulsion Laboratory/California Institute of Technology
    creator_url: https://earth.jpl.nasa.gov/emit/
    project: Earth Surface Mineral Dust Source Investigation
    project_url: https://earth.jpl.nasa.gov/emit/
    publisher_name: NASA LPDAAC
    publisher_url: https://lpdaac.usgs.gov
    publisher_email: lpdaac@usgs.gov
    identifier_product_doi_authority: https://doi.org
    flight_line: emit20241006t165200_o28011_s000
    time_coverage_start: 2024-10-06T16:52:00+0000
    time_coverage_end: 2024-10-06T16:52:12+0000
    software_build_version: 010619
    software_delivery_version: 010619
    product_version: V001
    history: PGE Run Command: {python /beegfs/store/emit/ops/repos/emit-sds-l2a/spectrum_quality.py /tmp/emit/ops/emit20241006t165200_emit.L2AReflectance_20241008t183054/output/emit20241006t165200_rfl /tmp/emit/ops/emit20241006t165200_emit.L2AReflectance_20241008t183054/output/emit20241006t165200_rfl_quality.txt}, PGE Input Files: {radiance_file=/beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_rdn_b0106_v01.img, pixel_locations_file=/beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_loc_b0106_v01.img, observation_parameters_file=/beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_obs_b0106_v01.img, surface_model_config=/beegfs/store/emit/ops/repos/emit-sds-l2a/surface/surface_20221020.json}
    crosstrack_orientation: as seen on ground
    easternmost_longitude: -114.70034790991087
    northernmost_latitude: 33.3865592156276
    westernmost_longitude: -115.94639824146
    southernmost_latitude: 32.36390868242409
    spatialResolution: 0.000542232520256367
    spatial_ref: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
    geotransform: [-1.15946398e+02  5.42232520e-04 -0.00000000e+00  3.33865592e+01
    -0.00000000e+00 -5.42232520e-04]
    day_night_flag: Day
    title: EMIT L2A Estimated Surface Reflectance 60 m V001
    Groups:
    Group: sensor_band_parameters
        Variables:
        wavelengths: shape=(285,), dtype=float32
        fwhm: shape=(285,), dtype=float32
        good_wavelengths: shape=(285,), dtype=uint8
    Group: location
        Variables:
        lon: shape=(1280, 1242), dtype=float64
        lat: shape=(1280, 1242), dtype=float64
        elev: shape=(1280, 1242), dtype=float64
        glt_x: shape=(1886, 2298), dtype=int32
        glt_y: shape=(1886, 2298), dtype=int32
    """
    def __repr__(self) -> str:
        """
        Return a string representation of the EMITL2ARFLNetCDF object.
        """
        return f"EMITL2ARFLNetCDF(\"{self.filename}\")"

class EMITL2AMASKNetCDF(EMITNetCDF):
    """
    This class encapsulates access to the EMIT L2A Reflectance Mask NetCDF file.

    Example EMIT L2A Reflectance Mask NetCDF file structure:

    File: /Users/halverso/data/EMIT_L2A_RFL/EMIT_L2A_RFL_001_20241006T165200_2428011_004/EMIT_L2A_MASK_001_20241006T165200_2428011_004.nc
    Dimensions:
    downtrack: size=1280
    crosstrack: size=1242
    bands: size=8
    ortho_y: size=1886
    ortho_x: size=2298
    packed_wavelength_bands: size=36
    Variables:
    mask: shape=(1280, 1242, 8), dtype=float32
    band_mask: shape=(1280, 1242, 36), dtype=uint8
    Attributes:
    ncei_template_version: NCEI_NetCDF_Swath_Template_v2.0
    summary: The Earth Surface Mineral Dust Source Investigation (EMIT) is an Earth Ventures-Instrument (EVI-4) Mission that maps the surface mineralogy of arid dust source regions via imaging spectroscopy in the visible and short-wave infrared (VSWIR). Installed on the International Space Station (ISS), the EMIT instrument is a Dyson imaging spectrometer that uses contiguous spectroscopic measurements from 410 to 2450 nm to resolve absoprtion features of iron oxides, clays, sulfates, carbonates, and other dust-forming minerals. During its one-year mission, EMIT will observe the sunlit Earth's dust source regions that occur within +/-52° latitude and produce maps of the source regions that can be used to improve forecasts of the role of mineral dust in the radiative forcing (warming or cooling) of the atmosphere.\n\nThis file contains masks for L2A estimated surface reflectances and geolocation data. Masks account for clouds, cloud shadows (via buffering), spacecraft interference, and poor atmospheric conditions. Geolocation data (latitude, longitude, height) and a lookup table to project the data are also included.
    keywords: Imaging Spectroscopy, minerals, EMIT, dust, radiative forcing
    Conventions: CF-1.63
    sensor: EMIT (Earth Surface Mineral Dust Source Investigation)
    instrument: EMIT
    platform: ISS
    institution: NASA Jet Propulsion Laboratory/California Institute of Technology
    license: https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/
    naming_authority: LPDAAC
    date_created: 2024-10-09T02:59:49Z
    keywords_vocabulary: NASA Global Change Master Directory (GCMD) Science Keywords
    stdname_vocabulary: NetCDF Climate and Forecast (CF) Metadata Convention
    creator_name: Jet Propulsion Laboratory/California Institute of Technology
    creator_url: https://earth.jpl.nasa.gov/emit/
    project: Earth Surface Mineral Dust Source Investigation
    project_url: https://earth.jpl.nasa.gov/emit/
    publisher_name: NASA LPDAAC
    publisher_url: https://lpdaac.usgs.gov
    publisher_email: lpdaac@usgs.gov
    identifier_product_doi_authority: https://doi.org
    flight_line: emit20241006t165200_o28011_s000
    time_coverage_start: 2024-10-06T16:52:00+0000
    time_coverage_end: 2024-10-06T16:52:12+0000
    software_build_version: 010619
    software_delivery_version: 010619
    product_version: V001
    history: PGE Run Command: {python /beegfs/store/emit/ops/repos/emit-sds-l2a/make_emit_masks.py /beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_rdn_b0106_v01.img /beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_loc_b0106_v01.img /beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l2a/emit20241006t165200_o28011_s000_l2a_atm_b0106_v01.img /beegfs/store/emit/ops/repos/emit-sds-l2a/data/kurudz_0.1nm.dat /beegfs/scratch/emit/ops/tmp/emit20241006t165200_emit.L2AMask_20241008t193231/output/emit20241006t165200_o28011_s000_l2a_mask_b0106_v01.img --n_cores 40}, PGE Input Files: {radiance_file=/beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_rdn_b0106_v01.img, pixel_locations_file=/beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_loc_b0106_v01.img, atmosphere_file=/beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l2a/emit20241006t165200_o28011_s000_l2a_atm_b0106_v01.img, solar_irradiance_file=/beegfs/store/emit/ops/repos/emit-sds-l2a/data/kurudz_0.1nm.dat}
    crosstrack_orientation: as seen on ground
    easternmost_longitude: -114.70034790991087
    northernmost_latitude: 33.3865592156276
    westernmost_longitude: -115.94639824146
    southernmost_latitude: 32.36390868242409
    spatialResolution: 0.000542232520256367
    spatial_ref: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
    geotransform: [-1.15946398e+02  5.42232520e-04 -0.00000000e+00  3.33865592e+01
    -0.00000000e+00 -5.42232520e-04]
    day_night_flag: Day
    title: EMIT L2A Masks 60 m V001
    Groups:
    Group: sensor_band_parameters
        Variables:
        mask_bands: shape=(8,), dtype=<class 'str'>
    Group: location
        Variables:
        lon: shape=(1280, 1242), dtype=float64
        lat: shape=(1280, 1242), dtype=float64
        elev: shape=(1280, 1242), dtype=float64
        glt_x: shape=(1886, 2298), dtype=int32
        glt_y: shape=(1886, 2298), dtype=int32
    """
    def __repr__(self) -> str:
        """
        Return a string representation of the EMITL2AMASKNetCDF object.
        """
        return f"EMITL2AMASKNetCDF(\"{self.filename}\")"

class EMITL2ARFLUNCERTNetCDF(EMITNetCDF):
    """
    This class encapsulates access to the EMIT L2A Reflectance Uncertainty NetCDF file.
    
    Example EMIT L2A Reflectance Uncertainty NetCDF file structure:

    File: /Users/halverso/data/EMIT_L2A_RFL/EMIT_L2A_RFL_001_20241006T165200_2428011_004/EMIT_L2A_RFLUNCERT_001_20241006T165200_2428011_004.nc
    Dimensions:
    downtrack: size=1280
    crosstrack: size=1242
    bands: size=285
    ortho_y: size=1886
    ortho_x: size=2298
    Variables:
    reflectance_uncertainty: shape=(1280, 1242, 285), dtype=float32
    Attributes:
    ncei_template_version: NCEI_NetCDF_Swath_Template_v2.0
    summary: The Earth Surface Mineral Dust Source Investigation (EMIT) is an Earth Ventures-Instrument (EVI-4) Mission that maps the surface mineralogy of arid dust source regions via imaging spectroscopy in the visible and short-wave infrared (VSWIR). Installed on the International Space Station (ISS), the EMIT instrument is a Dyson imaging spectrometer that uses contiguous spectroscopic measurements from 410 to 2450 nm to resolve absoprtion features of iron oxides, clays, sulfates, carbonates, and other dust-forming minerals. During its one-year mission, EMIT will observe the sunlit Earth's dust source regions that occur within +/-52° latitude and produce maps of the source regions that can be used to improve forecasts of the role of mineral dust in the radiative forcing (warming or cooling) of the atmosphere.\n\nThis file contains L2A estimated surface reflectance uncertainties and geolocation data. Reflectance uncertainty estimates are created using an Optimal Estimation technique - see ATBD for details. Reflectance uncertainty values are reported as fractions (relative to 1). Geolocation data (latitude, longitude, height) and a lookup table to project the data are also included.
    keywords: Imaging Spectroscopy, minerals, EMIT, dust, radiative forcing
    Conventions: CF-1.63
    sensor: EMIT (Earth Surface Mineral Dust Source Investigation)
    instrument: EMIT
    platform: ISS
    institution: NASA Jet Propulsion Laboratory/California Institute of Technology
    license: https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/
    naming_authority: LPDAAC
    date_created: 2024-10-09T02:59:28Z
    keywords_vocabulary: NASA Global Change Master Directory (GCMD) Science Keywords
    stdname_vocabulary: NetCDF Climate and Forecast (CF) Metadata Convention
    creator_name: Jet Propulsion Laboratory/California Institute of Technology
    creator_url: https://earth.jpl.nasa.gov/emit/
    project: Earth Surface Mineral Dust Source Investigation
    project_url: https://earth.jpl.nasa.gov/emit/
    publisher_name: NASA LPDAAC
    publisher_url: https://lpdaac.usgs.gov
    publisher_email: lpdaac@usgs.gov
    identifier_product_doi_authority: https://doi.org
    flight_line: emit20241006t165200_o28011_s000
    time_coverage_start: 2024-10-06T16:52:00+0000
    time_coverage_end: 2024-10-06T16:52:12+0000
    software_build_version: 010619
    software_delivery_version: 010619
    product_version: V001
    history: PGE Run Command: {python /beegfs/store/emit/ops/repos/emit-sds-l2a/spectrum_quality.py /tmp/emit/ops/emit20241006t165200_emit.L2AReflectance_20241008t183054/output/emit20241006t165200_rfl /tmp/emit/ops/emit20241006t165200_emit.L2AReflectance_20241008t183054/output/emit20241006t165200_rfl_quality.txt}, PGE Input Files: {radiance_file=/beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_rdn_b0106_v01.img, pixel_locations_file=/beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_loc_b0106_v01.img, observation_parameters_file=/beegfs/store/emit/ops/data/acquisitions/20241006/emit20241006t165200/l1b/emit20241006t165200_o28011_s000_l1b_obs_b0106_v01.img, surface_model_config=/beegfs/store/emit/ops/repos/emit-sds-l2a/surface/surface_20221020.json}
    crosstrack_orientation: as seen on ground
    easternmost_longitude: -114.70034790991087
    northernmost_latitude: 33.3865592156276
    westernmost_longitude: -115.94639824146
    southernmost_latitude: 32.36390868242409
    spatialResolution: 0.000542232520256367
    spatial_ref: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
    geotransform: [-1.15946398e+02  5.42232520e-04 -0.00000000e+00  3.33865592e+01
    -0.00000000e+00 -5.42232520e-04]
    day_night_flag: Day
    title: EMIT L2A Estimated Surface Reflectance Uncertainty 60 m V001
    Groups:
    Group: sensor_band_parameters
        Variables:
        wavelengths: shape=(285,), dtype=float32
        fwhm: shape=(285,), dtype=float32
        good_wavelengths: shape=(285,), dtype=uint8
    Group: location
        Variables:
        lon: shape=(1280, 1242), dtype=float64
        lat: shape=(1280, 1242), dtype=float64
        elev: shape=(1280, 1242), dtype=float64
        glt_x: shape=(1886, 2298), dtype=int32
        glt_y: shape=(1886, 2298), dtype=int32
    """
    def __repr__(self) -> str:
        """
        Return a string representation of the EMITL2ARFLUNCERTNetCDF object.
        """
        return f"EMITL2ARFLUNCERTNetCDF(\"{self.filename}\")"

class EMITL2ARFLGranule:
    def __init__(self, reflectance_filename: str, mask_filename: str, uncertainty_filename: str) -> None:
        """
        Initialize an EMITL2ARFL object for a set of EMIT L2A NetCDF files.

        Args:
            reflectance_filename (str): Path to the reflectance NetCDF file.
            mask_filename (str): Path to the mask NetCDF file.
            uncertainty_filename (str): Path to the uncertainty NetCDF file.
        """
        self.reflectance_filename: str = abspath(expanduser(reflectance_filename))
        self.mask_filename: str = abspath(expanduser(mask_filename))
        self.uncertainty_filename: str = abspath(expanduser(uncertainty_filename))

    def __repr__(self) -> str:
        """
        Return a string representation of the EMITL2ARFL object.
        """
        return (f"EMITL2ARFL(reflectance_filename=\"{self.reflectance_filename}\", "
                f"mask_filename=\"{self.mask_filename}\", "
                f"uncertainty_filename=\"{self.uncertainty_filename}\")")

    @property
    def reflectance_netcdf(self) -> EMITL2ARFLNetCDF:
        """
        Get the EMITL2ARFLNetCDF object for the reflectance file.

        Returns:
            EMITL2ARFLNetCDF: Object for reflectance NetCDF file.
        """
        return EMITL2ARFLNetCDF(self.reflectance_filename)

    @property
    def mask_netcdf(self) -> EMITL2AMASKNetCDF:
        """
        Get the EMITL2AMASKNetCDF object for the mask file.

        Returns:
            EMITL2AMASKNetCDF: Object for mask NetCDF file.
        """
        return EMITL2AMASKNetCDF(self.mask_filename)

    @property
    def uncertainty_netcdf(self) -> EMITL2ARFLUNCERTNetCDF:
        """
        Get the EMITL2ARFLUNCERTNetCDF object for the uncertainty file.

        Returns:
            EMITL2ARFLUNCERTNetCDF: Object for uncertainty NetCDF file.
        """
        return EMITL2ARFLUNCERTNetCDF(self.uncertainty_filename)

    @property
    def lat(self) -> np.ndarray:
        """
        Get the latitude array from the reflectance NetCDF file.

        Returns:
            np.ndarray: Latitude values.
        """
        return self.reflectance_netcdf.lat

    @property
    def lon(self) -> np.ndarray:
        """
        Get the longitude array from the reflectance NetCDF file.

        Returns:
            np.ndarray: Longitude values.
        """
        return self.reflectance_netcdf.lon

    @property
    def geolocation(self) -> RasterGeolocation:
        """
        Get the geolocation object from the reflectance NetCDF file.

        Returns:
            RasterGeolocation: Geolocation object.
        """
        return self.reflectance_netcdf.geolocation

    def extract_GLT(self) -> GLT:
        """
        Extract the Geometry Lookup Table (GLT) from the reflectance NetCDF file.

        Returns:
            GLT: Geometry Lookup Table object.
        """
        return self.reflectance_netcdf.extract_GLT()

    GLT = property(extract_GLT)

    def quality_mask(self, quality_bands: List[int] = QUALITY_BANDS) -> np.ndarray:
        """
        Generate a quality mask from the mask NetCDF file for the specified bands.

        Args:
            quality_bands (List[int], optional): List of quality band indices. Defaults to QUALITY_BANDS.

        Returns:
            np.ndarray: Quality mask array.
        """
        qmask: np.ndarray = quality_mask(
            filepath=self.mask_filename,
            quality_bands=quality_bands
        )
        return qmask

    def reflectance(self, geometry: RasterGeometry = None) -> Raster:
        """
        Get the reflectance raster, optionally reprojected to a given geometry.

        Args:
            geometry (RasterGeometry, optional): Target geometry for reprojection. Defaults to None.

        Returns:
            Raster: Reflectance raster object.
        """
        # Generate a quality mask for filtering bad pixels
        qmask: np.ndarray = self.quality_mask()
        # Create the reflectance raster, applying the quality mask
        raster: Raster = emit_ortho_raster(
            filepath=self.reflectance_filename,
            layer_name="reflectance",
            qmask=qmask
        )
        # Optionally reproject to the specified geometry
        if geometry is not None:
            raster = raster.to_geometry(geometry)
        return raster
