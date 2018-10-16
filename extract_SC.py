var seg = ee.FeatureCollection("users/alexbenemilles/segments_GEE"),
    precip_daily = ee.ImageCollection("NASA/GLDAS/V021/NOAH/G025/T3H").select('Rainf_tavg');

var precip_daily = precip_daily.filter(ee.Filter.date('2013-07-01', '2015-12-31'));
//print(precip_daily);
var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var cbuffer = function(feature) {
  return feature.buffer(10000);
};

var empty2 = ee.FeatureCollection([
  ee.Feature(null, {
    ID: 'has0', 
    p : 300.1
  }),
  ee.Feature(null, {
    ID: 'has1', 
    p : 300.1
  })]);
  

var start = ee.Date('2014-01-01').millis()
var end = ee.Date('2014-01-01').millis()
var iter_dates = ee.List.sequence(start, end, 86400*1000)


var out = iter_dates.iterate(function(date, result) {
  var start = ee.Number(date).divide(1000000000000)
  var lod = ee.Number(86400*1000).divide(1000000000000)
  var end = start.add(lod)
  var features = seg.filter(ee.Filter.rangeContains('time', start, end)).map(getCentroids).map(cbuffer)
  var filter_image = precip_daily.filterDate(ee.Date(date).advance(-120, 'day'), ee.Date(date));
  var multiband = filter_image.sum()
  
  var ext = ee.Image(multiband).reduceRegions({
    collection: features,
    reducer:ee.Reducer.mean(),
    scale: 1000,
    tileScale: 1
  })
  
  return ext
}, empty2)

var prep = out

var start = ee.Date('2014-01-01').millis()
var end = ee.Date('2015-12-31').millis()
var iter_dates = ee.List.sequence(start, end, 86400*1000)

var out = iter_dates.iterate(function(date, result) {
  var start = ee.Number(date).divide(1000000000000)
  var lod = ee.Number(86400*1000).divide(1000000000000)
  var end = start.add(lod)
  var features = seg.filter(ee.Filter.rangeContains('time', start, end)).map(getCentroids).map(cbuffer)
  var filter_image = precip_daily.filterDate(ee.Date(date).advance(-120, 'day'), ee.Date(date));
  var multiband = filter_image.sum().multiply(3600*3)
  
  var ext = ee.Image(multiband).reduceRegions({
    collection: features,
    reducer:ee.Reducer.mean(),
    scale: 1000,
    tileScale: 16
  })
  
  return ee.FeatureCollection(result).merge(ee.FeatureCollection(ext))
}, prep)

//print(ee.FeatureCollection(out).limit(1))

Export.table.toDrive({
    collection: out,
    description: 'SC',
    folder: 'GEC',
    fileFormat: '',
    selectors: ['ID', 'mean']
});