var seg = ee.FeatureCollection("users/alexbenemilles/segments_GEE"),
    precip_daily = ee.ImageCollection("NCEP_RE/surface_temp");

//var seg = seg.limit(10);
var precip_daily = precip_daily.filter(ee.Filter.date('2013-01-01', '2016-12-31'));
//print(precip_daily);
var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var centroids = seg.map(getCentroids);

var centroids_buffer = centroids.map(function(feature) {
  return feature.buffer(5000);
});


// create dummies to initialize datasets 
var features = seg.filter(ee.Filter.rangeContains('time', 1, 1))
var filter_image = precip_daily.filterDate(ee.Date(86400*365*1000*45), ee.Date(ee.Number(86400*365*1000*45).add(ee.Number(86400*1000 + 1))));
  
var empty = ee.Image().select()  
var multiband = filter_image.iterate(function(image, result) {
   return ee.Image(result).addBands(image)
}, empty)


var prep = ee.Image(multiband).reduceRegions({
  collection: features,
  reducer:ee.Reducer.mean(),
  scale: 100,
  tileScale: 1
})
//

var start = ee.Date('2014-01-01').millis()
var end = ee.Date('2015-12-31').millis()
var iter_dates = ee.List.sequence(start, end, 86400*1000)

var out = iter_dates.iterate(function(date, result) {
  var start = ee.Number(date).divide(1000000000000)
  var lod = ee.Number(86400*1000).divide(1000000000000)
  var end = start.add(lod)
  var features = seg.filter(ee.Filter.rangeContains('time', start, end))
  var filter_image = precip_daily.filterDate(ee.Date(date), ee.Date(ee.Number(date).add(ee.Number(86400*1000 + 1))));

  var multiband = filter_image.max()

 // var multiband = filter_image.iterate(function(image, result) {
  // return ee.Image(result).addBands(image)
//}, empty)

  //multiband = ee.Image(multiband).max()
  
  var ext = ee.Image(multiband).reduceRegions({
    collection: features,
    reducer:ee.Reducer.mean(),
    scale: 1000,
    tileScale: 1
  })
  
  return ee.FeatureCollection(result).merge(ee.FeatureCollection(ext))
}, prep)

Export.table.toDrive({
    collection: out,
    description: 'TC',
    folder: 'GEC',
    fileFormat: '',
    selectors: ['ID', 'mean']
});

