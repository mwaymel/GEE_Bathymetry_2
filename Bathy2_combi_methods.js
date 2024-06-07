//This is a Javascript code for running in Google Earth Engine

//Authors: Mathilde Waymel for Kartverket
//Inspired by Jiwei Li et al. (2021) code: https://github.com/CoralMapping/GEE_Sentinel2_Bathymetry_Paper/tree/main
//from the article Automated Global Shallow Water Bathymetry Mapping Using Google Earth Engine (https://doi.org/10.3390/rs13081469)

//___________________________________________________________________________________________________
//___________________________________________________________________________________________________

//              This code generates depth map in shallow waters using Sentinel2 and Landsat8 data
//                Map goes up to 50 meters from the coast in the region of interest


//_________________________________________PARAMETERS TO SET_________________________________________

//Region of interest coordinates (format: list of x, y - Example: [x1, y1, x2, y2, x3, y3])
//Example Fjoloy area: var ROI_coords = [302400, 6558000,307200, 6558000,307200, 6553800,302400, 6553800]
var ROI_coords = [302400, 6558000,
                  307200, 6558000,
                  307200, 6553800,
                  302400, 6553800];

//Name of the files and folder for Google Drive exporting  
var name_export = 'Bathymetry_Fjoloy';
var export_folder = 'GEE_export';
var scale_export = 10; //In meters, 10 = native RGB Sentinel 2 resolution, 30 = native RGB Landsat 8

//Smoothing kernel radius (default value 7)
var radius_kernel = 7;

//_______Ground truth available?______ 
//If so, import it in Assets and replace the path below
var name_ground_truth = 'users/majuwa/bathy_norge/vector_Fjoloy_truth_0_10';
//If not, comment lines 214 to 275, using /* at the beginning and */ at the end, 
//and from line 275, comment all lines using ground_truth and difference

var name_Mean_sea_correction = 'users/majuwa/bathy_norge/MeanSeaLevel1996-2014_above_NN2000_v2023b';

//Comment or decomment on lines 275-295 to choose the layers displayed


//___________________________________________________________________________________________________
//____________________________________COLLECTING AND PROCESSING DATA_________________________________


var roi = ee.Algorithms.GeometryConstructors.Polygon(
      ROI_coords,
      'EPSG:25832', false
  ); 
  
//set the filter input data to Sentinel-2 depth data
var sentinel = ee.ImageCollection('COPERNICUS/S2_SR').filter(ee.Filter.bounds(roi));

//Set up the date range and filter, this example uses two years window to build the clean water mosaic
sentinel = sentinel.filter(ee.Filter.date(ee.Date.fromYMD(2021,1,1),ee.Date.fromYMD(2021,12,31)));

//Import landsat 8 data
var landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filter(ee.Filter.bounds(roi));
landsat = landsat.filter(ee.Filter.date(ee.Date.fromYMD(2021,1,1),ee.Date.fromYMD(2021,12,31)));


//_______________MASKING FUNCTION_________________
//building the clean mosiac image based on different filters
var cloudBitMask = ee.Number(2).pow(10).int();
var cirrusBitMask = ee.Number(2).pow(11).int();

//this function is used to build clean water mosaic in the Google Earth Engine
//the threshold value could be revised, the current value is suggested for a common clean coral reefs waters
function mask_sentinel(img){
  var qa = img.select('QA60');
  var ma = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
  ma = ma.and(img.select(['SCL']).neq(3));
  ma = ma.and(img.select(['SCL']).neq(4));
  ma = ma.and(img.select(['SCL']).neq(5));
  ma = ma.and(img.select(['SCL']).neq(8));
  ma = ma.and(img.select(['SCL']).neq(9));
  ma = ma.and(img.select(['SCL']).neq(10));
  //ma = ma.and(img.select(['B9']).lt(300));
  //ma = ma.and(img.select(['B9']).gt(50));
  //ma = ma.and(img.select(['B3']).gt(100));//.focal_min({kernel: ee.Kernel.circle({radius: 5}), iterations: 1}));
  //ma = ma.focal_min({kernel: ee.Kernel.circle({radius: 1}), iterations: 1});
  img = img.mask(ma);
 
  //adjust for mask bad data
    img = img.updateMask(img.select([4]).lt(1000));
    img = img.updateMask(img.select([7]).lt(300));
  
  var ndwi_revise = (img.select([2]).subtract(img.select([7]))).divide(img.select([2]).add(img.select([7])));
  img = img.updateMask(ndwi_revise.gt(0));

  return img;
}

function mask_landsat(image) {
    // Develop masks for unwanted pixels (fill, cloud, cloud shadow).
    var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0); //px without cloud, cloud shadow and snow
    var saturationMask = image.select('QA_RADSAT').eq(0); //only pixels QA flagged as without saturation
  
    // Apply the scaling factors to the appropriate bands.
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  
    // Replace the original bands with the scaled ones and apply the masks.
    return image.addBands(opticalBands, null, true)//;
        .updateMask(qaMask)
        .updateMask(saturationMask)
        .select('SR_B.')
        .rename('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7');
  }

//run the mask function
sentinel = sentinel.map(mask_sentinel).select('B.');
landsat = landsat.map(mask_landsat);

var proj = roi.projection().getInfo();
var crs = proj['crs'];

function resample_landsat(img) {
  return img.resample('bilinear').reproject({'crs': crs, 'scale': 10.0});
}
var landsat_10m = landsat.map(resample_landsat);
//print(landsat_10m.first());
//print(sentinel.first());

//_____________2ND MASKING OF LAND________________

function mask2(img){
  var qa = img.select('QA60');
  var ma = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
  img = img.mask(ma);
  return img;
}
var sentinel2 = ee.ImageCollection('COPERNICUS/S2_SR').filter(ee.Filter.bounds(roi));
sentinel2 = sentinel2.filter(ee.Filter.date(ee.Date.fromYMD(2021,1,1),ee.Date.fromYMD(2021,12,31)));
sentinel2 = sentinel2.map(mask2);
var median2 = sentinel2.reduce(ee.Reducer.median());
var ndwi_land_mask = (median2.select([2]).subtract(median2.select([7]))).divide(median2.select([2]).add(median2.select([7])));

function mask_land_2(img){
  img = img.updateMask(ndwi_land_mask.gt(0));
  return img;
}
//run the 2nd land mask function
sentinel = sentinel.map(mask_land_2);
landsat_10m = landsat_10m.map(mask_land_2);


//_______DEFINING THE COASTAL AREA (50m to the coast and holes filled)________


//Building of a land polygon to be able to use buffer
var land = ee.Image(1).mask(ndwi_land_mask.lte(0));
var land_vector = land.addBands(land).reduceToVectors({
  geometry: roi,
  crs: roi.projection(),
  scale: scale_export,
  geometryType: 'polygon',
  eightConnected: true,
  labelProperty: 'mask',
  reducer: ee.Reducer.min()
});
land_vector = land_vector.select(['1']);
var geom_land = ee.FeatureCollection(land_vector).geometry();
geom_land = geom_land.geometries();

//Removing land areas corresponding to 1 px = 100m2 (artefacts)
geom_land = ee.FeatureCollection(ee.List(geom_land).map(function(geom){
    var geom_area = ee.Geometry(geom).area(1);
    return ee.Feature(ee.Geometry(geom)).set('area', geom_area);
  }));
//Filter above 110m2 to be robust to vectorisation imperfections
var geom_land = geom_land.filter(ee.Filter.gt('area', 110));     

//Applying 2 buffers (200m and -150m) to fill holes and sea corridors up to 400m wide
var buffer_land = geom_land.geometry().buffer(200).buffer(-150);
var coastal_area = buffer_land.difference({'right': geom_land, 'maxError': 1});


Export.table.toDrive({
  collection: ee.FeatureCollection(ee.Feature(ee.Geometry(coastal_area))),
  description: 'coastal_area_Fjoloy',
  folder: 'Fjoloy',
  fileFormat: 'SHP'
});


//Import of area if already exist
//var coastal_area = ee.FeatureCollection('users/majuwa/bathy_norge/coastal_area_Fjoloy');

//THE AREA OF STUDY IS coastal_area

function clip_coastal_area(img){
  img = img.clip(coastal_area);
  return img;
}

//_____________________________________BUILDING GROUD TRUTH________________________________________

var ground_truth = ee.FeatureCollection(name_ground_truth);

var ground_truth_raster = ground_truth.reduceToImage({
    properties: ['Depth'],
    reducer: ee.Reducer.first()
});

//___________________________________________________________________________________________________
//___________________________________________COMPOSITE METHOD________________________________________

//get the median value of it
var s2_composite_median = sentinel.reduce(ee.Reducer.median());
var l8_composite_median = landsat_10m.reduce(ee.Reducer.median());

//______SENTINEL_______
//calculate the big Rrs, rrs,and rrs*1000
var s2_comp_bigrrs = s2_composite_median.divide(ee.Number(31415.926));
var s2_comp_rrsvec = s2_comp_bigrrs.divide((s2_comp_bigrrs.multiply(ee.Number(1.7))).add(ee.Number(0.52)));
var s2_comp_rrsvec1k = s2_comp_rrsvec.multiply(ee.Number(1000));
//calculate rrs vec
var s2_comp_lnrrsvec = s2_comp_rrsvec1k.log();
//set the chla value for depth processing
//chla median value for our area: 0.3012
var s2_comp_chla = -0.38;
var s2_scale_factor_m1 = 2.1;

var s2_comp_m0 = ee.Number(52.083 * Math.exp(0.957*s2_comp_chla));    //m0 = 69.40359032247949 while mean of our tests except 0: 69,306
var s2_comp_m1 = ee.Number(50.156 * Math.exp(0.957*s2_comp_chla)).subtract(s2_scale_factor_m1);    //m1 = 66.83575209212759 while mean of our tests except 0: 66,758

var s2_composite_depth = ((s2_comp_lnrrsvec.select([1]).divide(s2_comp_lnrrsvec.select([2]))).multiply(s2_comp_m0)).subtract(s2_comp_m1);

//Correction of Chart Datum
var MeanSeaLevel2NN2000 = ee.Image(name_Mean_sea_correction).select(['b1']);
s2_composite_depth = s2_composite_depth.subtract(MeanSeaLevel2NN2000);


//Setting values between 0 and 20m deep (setting values below 0 to 0, and values above 20 to 20):
var s2_composite_depth = s2_composite_depth.where(s2_composite_depth.lt(0), ee.Number(0));
var s2_composite_depth = s2_composite_depth.where(s2_composite_depth.gt(10), ee.Number(10));

var Kernel = ee.Kernel.gaussian({radius: radius_kernel, sigma:3});
var s2_comp_depth_smooth = s2_composite_depth.convolve(Kernel);
s2_comp_depth_smooth = s2_comp_depth_smooth.clip(coastal_area);

var s2_comp_diff_gt = ground_truth_raster.add(s2_comp_depth_smooth);


//______LANDSAT_______
//calculate the big Rrs, rrs,and rrs*1000
var l8_comp_bigrrs = l8_composite_median.divide(ee.Number(31415.926));
var l8_comp_rrsvec = l8_comp_bigrrs.divide((l8_comp_bigrrs.multiply(ee.Number(1.7))).add(ee.Number(0.52)));
var l8_comp_rrsvec1k = l8_comp_rrsvec.multiply(ee.Number(1000));
//calculate rrs vec
var l8_comp_lnrrsvec = l8_comp_rrsvec1k.log();
//set the chla value for depth processing
//chla median value for our area: 0.3012
var l8_comp_chla = 0.48;
var l8_scale_factor_m1 = -1.75;

var l8_comp_m0 = ee.Number(52.083 * Math.exp(0.957*l8_comp_chla));    //m0 = 69.40359032247949 while mean of our tests except 0: 69,306
var l8_comp_m1 = ee.Number(50.156 * Math.exp(0.957*l8_comp_chla)).subtract(l8_scale_factor_m1);    //m1 = 66.83575209212759 while mean of our tests except 0: 66,758

var l8_composite_depth = ((l8_comp_lnrrsvec.select([1]).divide(l8_comp_lnrrsvec.select([2]))).multiply(l8_comp_m0)).subtract(l8_comp_m1).subtract(10).multiply(-1);

//Correction of Chart Datum
var MeanSeaLevel2NN2000 = ee.Image(name_Mean_sea_correction).select(['b1']);
l8_composite_depth = l8_composite_depth.subtract(MeanSeaLevel2NN2000);


//Setting values between 0 and 20m deep (setting values below 0 to 0, and values above 20 to 20):
var l8_composite_depth = l8_composite_depth.where(l8_composite_depth.lt(0), ee.Number(0));
var l8_composite_depth = l8_composite_depth.where(l8_composite_depth.gt(10), ee.Number(10));

var Kernel = ee.Kernel.gaussian({radius: radius_kernel, sigma:3});
var l8_comp_depth_smooth = l8_composite_depth.convolve(Kernel);
l8_comp_depth_smooth = l8_comp_depth_smooth.clip(coastal_area);

var l8_comp_diff_gt = ground_truth_raster.add(l8_comp_depth_smooth);


//___________________________________________________________________________________________________
//_________________________________________LOW TURBIDITY METHOD______________________________________

//_______________SELECTION OF IMAGES WITH LOWEST TURBIDITY________________


var landsat_coastal_area = landsat_10m.map(clip_coastal_area);
/*
//Smoothing landsat resampled
var Kernel =ee.Kernel.gaussian({radius: radius_kernel, sigma:1});

var landsat_10m_smooth = landsat_coastal_area.map(function convolve_imgcoll(img) {
  return img.convolve(Kernel).clip(coastal_area);
});
*/
var landsat_10m_smooth = landsat_coastal_area;


var visualization = {bands: ['B4', 'B3', 'B2'],min: 0.0,max: 0.3,};
//Map.addLayer(landsat_10m_smooth.first(), visualization);

//calculating the turbidity index NDTI

var nb_total_pixel = ee.Number(sentinel.first()
    .unmask().mask(ndwi_land_mask.gt(0))
    .select('B3').reduceRegion({
      reducer: ee.Reducer.count(), geometry: coastal_area, scale: 10}).get('B3'));

function add_ndti(img){
  var ndti_band = img.normalizedDifference(['B4', 'B3']).rename('NDTI');
  var mean_ndti = ndti_band.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: coastal_area, scale: 10});
  mean_ndti = ee.Number(mean_ndti.get('NDTI'));
  var nb_px = ndti_band.reduceRegion({
      reducer: ee.Reducer.count(), geometry: coastal_area, scale: 10});
  var perc_px = ee.Number(nb_px.get('NDTI')).divide(nb_total_pixel);
  return img.addBands(ndti_band).set('mean_NDTI', mean_ndti, 'percent_unmasked_px', perc_px);
}

var sentinel_ndti = sentinel.map(add_ndti);
var landsat_ndti = landsat_10m_smooth.map(add_ndti);


sentinel_ndti = sentinel_ndti.filter(ee.Filter.notNull(['mean_NDTI']))
                           .filter(ee.Filter.gte('percent_unmasked_px', 0.9));
sentinel_ndti = sentinel_ndti.sort('mean_NDTI', true);
var sentinel_selection = sentinel_ndti.limit(23).sort('mean_NDTI', false).limit(20).sort('mean_NDTI', true);
var sentinel_coastal_area = sentinel_selection.map(clip_coastal_area);

landsat_ndti = landsat_ndti.filter(ee.Filter.notNull(['mean_NDTI']))
                           .filter(ee.Filter.gte('percent_unmasked_px', 0.75));
landsat_ndti = landsat_ndti.sort('mean_NDTI', true);
var landsat_selection = landsat_ndti.limit(20);
var landsat_coastal_area = landsat_selection.map(clip_coastal_area);

print('sentinel', sentinel_coastal_area);
print(sentinel_coastal_area.aggregate_array('system:index'));
print('landsat', landsat_coastal_area);
print(landsat_coastal_area.aggregate_array('system:index'));

var vizParams = {bands: ['B4', 'B3', 'B2'], min: 0, max: 3000};
//Map.addLayer(sentinel_coastal_area.first(), vizParams, 'RGB - Lowest NDTI image');
//Map.addLayer(sentinel_coastal_area.first().select(['NDTI']), {min: -1, max: 0, palette: ['0000FF', '00FFFF', '00FF00', 'FFFC00','FFCD00','FF0000']}, 'NDTI - Lowest NDTI image');
//Map.addLayer(sentinel_coastal_area.limit(2).sort('mean_NDTI', false).first().select(['NDTI']), {min: -1, max: 0, palette: ['0000FF', '00FFFF', '00FF00', 'FFFC00','FFCD00','FF0000']}, '2 NDTI - Lowest NDTI image');
//Map.addLayer(sentinel_coastal_area.limit(3).sort('mean_NDTI', false).first().select(['NDTI']), {min: -1, max: 0, palette: ['0000FF', '00FFFF', '00FF00', 'FFFC00','FFCD00','FF0000']}, '3 NDTI - Lowest NDTI image');
//Map.addLayer(sentinel_coastal_area.limit(4).sort('mean_NDTI', false).first().select(['NDTI']), {min: -1, max: 0, palette: ['0000FF', '00FFFF', '00FF00', 'FFFC00','FFCD00','FF0000']}, '4 NDTI - Lowest NDTI image');
//Map.addLayer(sentinel_coastal_area.sort('mean_NDTI', false).first().select(['NDTI']), {min: -1, max: 0, palette: ['0000FF', '00FFFF', '00FF00', 'FFFC00','FFCD00','FF0000']}, '10 NDTI - Lowest NDTI image');
//Map.addLayer(landsat_coastal_area.first().select(['NDTI']), {min: -1, max: 0, palette: ['0000FF', '00FFFF', '00FF00', 'FFFC00','FFCD00','FF0000']}, 'NDTI - Lowest NDTI image');
//Map.addLayer(landsat_coastal_area.limit(2).sort('mean_NDTI', false).first().select(['NDTI']), {min: -1, max: 0, palette: ['0000FF', '00FFFF', '00FF00', 'FFFC00','FFCD00','FF0000']}, '2 NDTI - Lowest NDTI image');



//_______________CALCULATING DEPTH________________

//get the median value of it
//var median = sentinel.reduce(ee.Reducer.median());

function compute_depth_sentinel_img(img) {
  //calculate the big Rrs, rrs,and rrs*1000
  var bigrrs = img.divide(ee.Number(31415.926));//(img.subtract(img.select([4]))).divide(ee.Number(3.1415926));
  var rrsvec = bigrrs.divide((bigrrs.multiply(ee.Number(1.7))).add(ee.Number(0.52)));
  var rrsvec1k = rrsvec.multiply(ee.Number(1000));
  
  //calculate rrs vec
  var lnrrsvec = rrsvec1k.log();
  
  
  //Tests to set the Chla value in the area: 
  var w = bigrrs.select([2]).subtract(bigrrs.select([3]).multiply(0.46)).subtract(bigrrs.select([1]).multiply(0.54));
  var Chla_map = ee.Image(10).pow(w.multiply(191.659).subtract(0.4909));
  var m0_map =  Chla_map.multiply(0.9).exp().multiply(15);//Chla_map.multiply(0.957).exp().multiply(52.083);
  var m1_map =  Chla_map.multiply(0.9).exp().multiply(20);//Chla_map.multiply(0.957).exp().multiply(50.156);
  
  //set the chla value for depth processing
  //chla median value for our area: 0.3012
  var s2_comp_chla = -0.38;
  var s2_scale_factor_m1 = 2.65;
  
  var m0 = ee.Number(52.083 * Math.exp(0.957*s2_comp_chla));    //m0 = 69.40359032247949 while mean of our tests except 0: 69,306
  var m1 = ee.Number(50.156 * Math.exp(0.957*s2_comp_chla)).subtract(s2_scale_factor_m1);    //m1 = 66.83575209212759 while mean of our tests except 0: 66,758
  
  var depth = ((lnrrsvec.select([1]).divide(lnrrsvec.select([2]))).multiply(m0)).subtract(m1);
  depth = depth.rename('depth');
  
  var depth_coastal_blue = ((lnrrsvec.select([0]).divide(lnrrsvec.select([2]))).multiply(m0)).subtract(m1);
  depth_coastal_blue = depth_coastal_blue.rename('depth_coastal_blue');
  //Correction of Chart Datum
//  var MeanSeaLevel2NN2000 = ee.Image(name_Mean_sea_correction).select(['b1']);
//  depth = depth.subtract(MeanSeaLevel2NN2000);
  
  return img.addBands(depth).addBands(depth_coastal_blue);//.set('m0', m0, 'm1', m1);
}

function compute_depth_landsat_img(img) {
  //calculate the big Rrs, rrs,and rrs*1000
  var bigrrs = img.divide(ee.Number(31415.926));//(img.subtract(img.select([4]))).divide(ee.Number(3.1415926));
  var rrsvec = bigrrs.divide((bigrrs.multiply(ee.Number(1.7))).add(ee.Number(0.52)));
  var rrsvec1k = rrsvec.multiply(ee.Number(1000));
  
  //calculate rrs vec
  var lnrrsvec = rrsvec1k.log();
  
  
  //Tests to set the Chla value in the area: 
  var w = bigrrs.select([2]).subtract(bigrrs.select([3]).multiply(0.46)).subtract(bigrrs.select([1]).multiply(0.54));
  var Chla_map = ee.Image(10).pow(w.multiply(191.659).subtract(0.4909));
  var m0_map =  Chla_map.multiply(0.9).exp().multiply(15);//Chla_map.multiply(0.957).exp().multiply(52.083);
  var m1_map =  Chla_map.multiply(0.9).exp().multiply(20);//Chla_map.multiply(0.957).exp().multiply(50.156);
  
  //set the chla value for depth processing
  //chla median value for our area: 0.3012
  var l8_comp_chla = 0.45;
  var l8_scale_factor_m1 = -2.2;
  
  var m0 = ee.Number(52.083 * Math.exp(0.957*l8_comp_chla));    //m0 = 69.40359032247949 while mean of our tests except 0: 69,306
  var m1 = ee.Number(50.156 * Math.exp(0.957*l8_comp_chla)).subtract(l8_scale_factor_m1);    //m1 = 66.83575209212759 while mean of our tests except 0: 66,758

  var depth = ((lnrrsvec.select([1]).divide(lnrrsvec.select([2]))).multiply(m0)).subtract(m1).subtract(10).multiply(-1);
  depth = depth.rename('depth');
  
  var depth_coastal_blue = ((lnrrsvec.select([0]).divide(lnrrsvec.select([2]))).multiply(m0)).subtract(m1);
  depth_coastal_blue = depth_coastal_blue.rename('depth_coastal_blue');
  //Correction of Chart Datum
//  var MeanSeaLevel2NN2000 = ee.Image(name_Mean_sea_correction).select(['b1']);
//  depth = depth.subtract(MeanSeaLevel2NN2000);
  
  return img.addBands(depth).addBands(depth_coastal_blue);//.set('m0', m0, 'm1', m1);
}

var sentinel_depth = sentinel_coastal_area.map(compute_depth_sentinel_img);
var landsat_depth = landsat_coastal_area.map(compute_depth_landsat_img);


//_______________________________________SMOOTHING (CONVOLVE)________________________________________

function correct_img(img) {
  img = img.where(img.lt(0), ee.Number(0));
  return img.where(img.gt(10), ee.Number(10));
}
//Map.addLayer(sentinel_depth.first().select(['depth']), {min: -5, max: 25, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'sentinel_depth_correct');
var sentinel_depth_correct = sentinel_depth.map(correct_img);
var landsat_depth_correct = landsat_depth.map(correct_img);

//Map.addLayer(sentinel_depth_correct.first().select(['depth']), {min: -5, max: 25, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'sentinel_depth_correct');
//sentinel_depth_correct = sentinel_depth_correct.map(function median(img) {return img.focalMedian(1)});
//Map.addLayer(sentinel_depth_correct.first().select(['depth']), {min: -5, max: 25, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'sentinel_depth_correct smooth');

//var Kernel = ee.Kernel.circle({radius: radius_kernel});   //circle-shaped boolean kernel
var Kernel =ee.Kernel.gaussian({radius: radius_kernel, sigma:3});

function smooth(img) {
  return img.convolve(Kernel);
}

var sentinel_depth_smooth = sentinel_depth_correct.map(smooth);
var landsat_depth_smooth = landsat_depth_correct.map(smooth);

//_______________________________________MEDIAN MAP________________________________________

var sentinel_depth_median = sentinel_depth_smooth.reduce(ee.Reducer.median());
//Map.addLayer(sentinel_depth_median.select(['depth_median']), {min: 0, max: 10, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'sentinel_depth_correct smooth');

var landsat_depth_median = landsat_depth_smooth.reduce(ee.Reducer.median());
//Map.addLayer(landsat_depth_median.select(['depth_median']), {min: 0, max: 10, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'landsat_depth_correct smooth');


//_______________________________COMPARISON GROUD TRUTH____________________________________

var sentinel_diff_gt = ground_truth_raster.add(sentinel_depth_median);
var landsat_diff_gt = ground_truth_raster.add(landsat_depth_median);

//Map.addLayer(sentinel_diff_gt.select(['depth_median']), {"min":-10,"max":10,"palette":["0000ff","ffffff","ff0000"]}, 'sentinel_diff_gt');
//Map.addLayer(landsat_diff_gt.select(['depth_median']), {"min":-10,"max":10,"palette":["0000ff","ffffff","ff0000"]}, 'landsat_diff_gt');


//___________________________________________________________________________________________________
//_____________________________________METHODS COMBINATION___________________________________________

/*___Available layers:___
s2_comp_depth_smooth
s2_comp_diff_gt
l8_comp_depth_smooth
l8_comp_diff_gt

sentinel_depth_median
sentinel_diff_gt
landsat_depth_median
landsat_diff_gt
*/

var depth_collection = ee.ImageCollection([ee.Image(s2_comp_depth_smooth).select(['B2_median']).rename('depth_median'), ee.Image(l8_comp_depth_smooth).select(['B2_median']).rename('depth_median'), ee.Image(sentinel_depth_median).select(['depth_median']), ee.Image(landsat_depth_median).select(['depth_median'])]);

var depth_combi_mean = depth_collection.reduce(ee.Reducer.mean());
var diff_combi_mean = ground_truth_raster.add(depth_combi_mean).round().int().clip(coastal_area);

var band_name = 'first';

var chartdiff_combi_mean = ui.Chart.image.histogram({image: diff_combi_mean.select([band_name]), region: roi, scale: scale_export})
      .setOptions({
          title: 'Combination of methods: Mean',
          series: { 0: {visibleInLegend: false}},
          hAxis: {title: 'Difference (meters)', titleTextStyle: {italic: false, bold: true}},
          vAxis: {
            //viewWindow: {min: 0, max: 0.7},
            title: 'Number of pixels',
            titleTextStyle: {italic: false, bold: true}
          },
          //lineWidth: 5,
          //colors: ['a50f15', 'fcae91'],
          //trendlines: { 0: {visibleInLegend: true}}
        });
print(chartdiff_combi_mean);


//___________________________________________________________________________________________________
//__________________________________________LAYERS DISPLAY___________________________________________
/*___Available layers:___
s2_comp_depth_smooth
s2_comp_diff_gt
l8_comp_depth_smooth
l8_comp_diff_gt

sentinel_depth_median
sentinel_diff_gt
landsat_depth_median
landsat_diff_gt

depth_combi_mean
diff_combi_mean
*/
Map.centerObject(roi, 14);

//Map.addLayer(depth_combi_median, {min: 0, max: 10, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'depth_combi_median');
Map.addLayer(depth_combi_mean, {min: 0, max: 10, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'depth_combi_mean');
Map.addLayer(sentinel_depth_median, {min: 0, max: 10, bands: 'depth_median', palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'sentinel_depth_median');
Map.addLayer(landsat_depth_median, {min: 0, max: 10, bands: 'depth_median', palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'landsat_depth_median');
Map.addLayer(s2_comp_depth_smooth, {min: 0, max: 10, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 's2_comp_depth_smooth');
Map.addLayer(l8_comp_depth_smooth, {min: 0, max: 10, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'l8_comp_depth_smooth');

var imageVisParam2 = {"opacity":1,"bands":["first"],"min":-10,"max":10,"palette":["0000ff","ffffff","ff0000"]};
var imageVisParam3 = {"opacity":1,"bands":["depth_median"],"min":-10,"max":10,"palette":["0000ff","ffffff","ff0000"]};
//Map.addLayer(diff_combi_median, imageVisParam2, 'diff_combi_median');
Map.addLayer(diff_combi_mean, imageVisParam2, 'diff_combi_mean');
Map.addLayer(sentinel_diff_gt, imageVisParam3, 'sentinel_diff_gt');
Map.addLayer(landsat_diff_gt, imageVisParam3, 'landsat_diff_gt');
Map.addLayer(s2_comp_diff_gt, imageVisParam2, 's2_comp_diff_gt');
Map.addLayer(l8_comp_diff_gt, imageVisParam2, 'l8_comp_diff_gt');

//Map.addLayer(depth_smooth, {min: 0, max: 10, palette: ['00FF00', '00FFFF', '00bfff','0080ff','0000FF', 'FF0000']}, 'Depth map smooth');

//__Difference with ground truth rounded
var imageVisParam2 = {"opacity":1,"bands":["first"],"min":-10,"max":10,"palette":["0000ff","ffffff","ff0000"]};
//Map.addLayer(difference_rd, imageVisParam2, 'Difference rounded');

//var imageVisParam = {"opacity":1,"bands":["first"],"min":-20,"max":0,"palette":["000093","00ffff"]};
//Map.addLayer(ground_truth_raster, imageVisParam, 'Ground truth');


/*
//__Depth map without smoothing:
Map.addLayer(depth_output, {min: 0, max: 15, palette: ['00FFFF', '0000FF']}, 'Depth map without smoothing');

//__Depth map with smoothing:
Map.addLayer(depth_smooth, {min: 0, max: 15, palette: ['00FFFF', '0000FF']}, 'Depth map smooth');

//__Depth map smooth rounded:
Map.addLayer(depth_smooth_rd, {min: 0, max: 15, palette: ['00FFFF', '0000FF']}, 'Depth map smooth rounded');

//__Ground truth raster:
var imageVisParam = {"opacity":1,"bands":["first"],"min":-20,"max":0,"palette":["000093","00ffff"]};
Map.addLayer(ground_truth_raster, imageVisParam, 'Ground truth');

//__Difference with ground truth rounded
var imageVisParam2 = {"opacity":1,"bands":["first"],"min":-20,"max":20,"palette":["0000ff","ffffff","ff0000"]};
Map.addLayer(difference_rd, imageVisParam2, 'Difference rounded');
*/

//___________________________________________________________________________________________________
//______________________________________________EXPORTS______________________________________________

var projection = roi.projection().getInfo(); 

/*
//__________________________PREPARING OUTPUT________________________________

//Raster
var raster_name_export = name_export + '_smooth_raster';

//Rounded raster
var raster_rd_name_export = name_export + '_smooth_rd_raster';
var depth_smooth_rd = depth_smooth.round().int();

//Vector
var vector_name_export = name_export + '_smooth_rd_vector';
var depth_smooth_rd_vector = depth_smooth_rd.addBands(depth_smooth_rd).reduceToVectors({
  geometry: coastal_area,
  crs: roi.projection(),
  scale: scale_export,
  geometryType: 'polygon',
  eightConnected: true,
  labelProperty: 'B2_median',
  reducer: ee.Reducer.min()
});

//Exporting
var difference_name_export = name_export + '_comparison_truth';

var difference_rd_vector = difference_rd.addBands(difference_rd).reduceToVectors({
  geometry: coastal_area,
  crs: roi.projection(),
  scale: scale_export,
  geometryType: 'polygon',
  eightConnected: true,
  labelProperty: 'B2_median',
  reducer: ee.Reducer.min()
});

//Depth raster basis
Export.image.toDrive({
  image: depth_output,
  description: name_export,
  folder: export_folder,
  scale: scale_export,
  crs: projection.crs,
  region: coastal_area
});

//Smooth raster
Export.image.toDrive({
  image: depth_smooth,
  description: raster_name_export,
  folder: export_folder,
  scale: scale_export,
  crs: projection.crs,
  region: coastal_area
});

//Rounded raster
Export.image.toDrive({
  image: depth_smooth_rd,
  description: raster_rd_name_export,
  folder: export_folder,
  scale: scale_export,
  crs: projection.crs,
  region: coastal_area
});

//Rounded vector
Export.table.toDrive({
  collection: depth_smooth_rd_vector,
  description: vector_name_export,
  folder: export_folder,
  fileFormat: 'SHP'
});

//Difference raster
Export.image.toDrive({
  image: difference,
  description: difference_name_export,
  folder: export_folder,
  scale: scale_export,
  crs: projection.crs,
  region: coastal_area
});

//Difference raster
Export.image.toDrive({
  image: difference_rd,
  description: difference_name_export + '_rd',
  folder: export_folder,
  scale: scale_export,
  crs: projection.crs,
  region: coastal_area
});

//Difference rounded vector
Export.table.toDrive({
  collection: depth_smooth_rd_vector,
  description: difference_name_export + '_rd_vector',
  folder: export_folder,
  fileFormat: 'SHP'
});

*/
