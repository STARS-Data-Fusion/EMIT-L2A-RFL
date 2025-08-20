from .EMITNetCDF import EMITNetCDF

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
