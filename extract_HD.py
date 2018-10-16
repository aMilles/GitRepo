var HD_img = ee.ImageCollection('WorldPop/POP').filterDate('2014-01-01', '2016-01-01');
var seg = ee.FeatureCollection('users/alexbenemilles/segments_GEE');
var africa = ee.FeatureCollection('users/alexbenemilles/africa_hull5000');
var HD = ee.ImageCollection('CIESIN/GPWv4/unwpp-adjusted-population-density').filterDate('2014-01-01', '2016-01-01');
//Create functions to buffer the segments and buffer them

Map.addLayer(HD)


var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var bufferPoly = function(feature) {
  return feature.buffer(5000);   // substitute in your value of Z here
};

var centroids = seg.map(getCentroids);

var bufferedcents = centroids.map(bufferPoly);



//convert worldpop into a multiband image, each band is named by the number of 
//var empty = ee.Image(HD_img.first());
//var multiband = HD_img.iterate(function(image, result){
//  var name = image.get('country');
//  var name2 = image.get('UNadj')
//  image = image.rename(ee.String(ee.String(name).cat(ee.String(name2))));
//  return ee.Image(result).addBands(image);
//}, empty);
//
//print(multiband)
//


var HD = ee.Image(HD.first()).reduceRegions({
  collection: bufferedcents,
  tileScale: 1,
  scale: 1000,
  reducer: ee.Reducer.mean()
});


print(ee.FeatureCollection(HD).limit(10))

Export.table.toDrive({
  collection: HD,
  folder: 'GEC',
  description: 'HD',
  selectors: ['ID', 'mean']
})


