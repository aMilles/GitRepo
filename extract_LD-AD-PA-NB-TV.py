var pa = ee.Image("users/alexbenemilles/PA"),
    tlu = ee.Image("users/alexbenemilles/TLU"),
    lc = ee.Image("ESA/GLOBCOVER_L4_200901_200912_V2_3"),
    seg = ee.FeatureCollection("users/alexbenemilles/segments_GEE"),
    srtm = ee.Image("CGIAR/SRTM90_V4");
    
//preprocess lc: convert globcover map to ratio of cropland
lc = lc.select(['landcover']);
var agriculture = lc.lt(0)
agriculture = agriculture.where(lc.lt(15), 1)
agriculture = agriculture.where(lc.eq(20), .65)
agriculture = agriculture.where(lc.eq(30), .35)

pa = pa.int().gt(0)
//create buffered centroids

var slope = ee.Terrain.slope(srtm)

var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var bufferPoly = function(feature) {
  return feature.buffer(1500);   // substitute in your value of Z here
};

var centroids = seg.map(getCentroids);

var bufferedcents = centroids.map(bufferPoly);

//extract values

var ext_tlu = tlu.reduceRegions({
  collection: bufferedcents,
  reducer: ee.Reducer.mean(),
  scale: 100,
  tileScale: 1
});

var ext_agri = agriculture.reduceRegions({
  collection: bufferedcents,
  reducer: ee.Reducer.mean(),
  scale: 100,
  tileScale: 1
});

var ext_pa = pa.reduceRegions({
  collection: seg,
  reducer: ee.Reducer.max(),
  scale: 100,
  tileScale: 1
});

var ext_nb = slope.reduceRegions({
  collection: seg,
  reducer: ee.Reducer.mean(),
  scale: 100,
  tileScale: 1
});

var ext_tv = srtm.reduceRegions({
  collection: seg,
  reducer: ee.Reducer.stdDev(),
  scale: 100,
  tileScale: 1
});

var ext_sl = srtm.reduceRegions({
  collection: seg,
  reducer: ee.Reducer.mean(),
  scale: 100,
  tileScale: 1
});


print(ext_pa.limit(10))
    
//export extracted values

Export.table.toDrive({
  collection: ext_tlu,
  description:'LD',
  fileFormat: '',
  selectors: ['ID', 'mean'],
  folder: 'GEC'
});

Export.table.toDrive({
  collection: ext_agri,
  description:'AI',
  fileFormat: '',
  selectors: ['ID' , 'mean'],
  folder: 'GEC'
});

Export.table.toDrive({
  collection: ext_pa,
  description:'PA',
  fileFormat: '',
  selectors: ['ID', 'max'],
  folder: 'GEC'
});

Export.table.toDrive({
  collection: ext_nb,
  description:'NB',
  fileFormat: '',
  selectors: ['ID', 'mean'],
  folder: 'GEC'
});

Export.table.toDrive({
  collection: ext_tv,
  description:'TV',
  fileFormat: '',
  selectors: ['ID', 'stdDev'],
  folder: 'GEC'
});

Export.table.toDrive({
  collection: ext_sl,
  description:'SL',
  fileFormat: '',
  selectors: ['ID', 'mean'],
  folder: 'GEC'
});