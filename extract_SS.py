var soil = ee.Image("users/alexbenemilles/exchangeable_NA"),
    seg = ee.FeatureCollection("users/alexbenemilles/segments_GEE");
    
var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var bufferPoly = function(feature) {
  return feature.buffer(1500);   // substitute in your value of Z here
};

var centroids = seg.map(getCentroids);

var bufferedcents = centroids.map(bufferPoly);

//extract values

var ext_NA = soil.reduceRegions({
  collection: bufferedcents,
  reducer: ee.Reducer.mean(),
  scale: 100,
  tileScale: 1
});

print(ext_NA.limit(10))
    
//export extracted values

Export.table.toDrive({
  collection: ext_NA,
  description:'NA',
  fileFormat: '',
  selectors: ['ID', 'mean'],
  folder: 'GEC'
});
