var table = ee.FeatureCollection("users/alexbenemilles/segments_GEE"),
    africa = ee.FeatureCollection("users/alexbenemilles/rough_africa");
var pr = ee.ImageCollection('NASA/NEX-GDDP').filterMetadata('scenario', 'equals', 'rcp85').select('pr').filterBounds(africa);
var pr_hist = ee.ImageCollection('NASA/NEX-GDDP').filterMetadata('scenario', 'equals', 'historical').select('pr').filterDate("1950-01-01", "1953-01-01").filterBounds(africa);
var temp = ee.ImageCollection('NASA/NEX-GDDP').filterMetadata('scenario', 'equals', 'rcp85').select('tasmax').filterBounds(africa);
var temp_hist = ee.ImageCollection('NASA/NEX-GDDP').filterMetadata('scenario', 'equals', 'historical').select('tasmax').filterDate("1950-01-01", "1953-01-01").filterBounds(africa);
var segs = ee.FeatureCollection('users/alexbenemilles/segments_GEE');



var dates = [ "2015-01-01", "2015-02-01", "2015-03-01", "2015-04-01",
"2015-05-01", "2015-06-01", "2015-07-01", "2015-08-01", "2015-09-01", "2015-10-01", "2015-11-01", "2015-12-01",
"2099-01-01", "2099-02-01", "2099-03-01", "2099-04-01", "2099-05-01", "2099-06-01", "2099-07-01", "2099-08-01",
"2099-09-01", "2099-10-01", "2099-11-01", "2099-12-01"]

var seq_dates = ee.List.sequence(12,35,1)
print(seq_dates)
var dates_hist = ["1951-01-01", "1951-02-01", "1951-03-01", "1951-04-01", "1951-05-01", "1951-06-01", "1951-07-01", "1951-08-01",
"1951-09-01", "1951-10-01", "1951-11-01", "1951-12-01"]

var seq_dates_hist = ee.List.sequence(0,11,1)

var empty = ee.Image().select()

Map.addLayer(ee.Image(pr.filterDate('2015-01-01', '2015-04-01').sum()))

var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var bufferPoly = function(feature) {
  return feature.buffer(5000);   // substitute in your value of Z here
};

var centroids = segs.map(getCentroids);

var bufferedcents = centroids.map(bufferPoly);

var image = ee.Image(pr.first())

var dummy = image.reduceRegions({
    collection: bufferedcents.limit(1),
    reducer: ee.Reducer.mean(),
    scale: 500,
    tileScale: 1
  })


//precipitation


var out = ee.List(seq_dates).iterate(function(date, result){
  date = ee.Date(ee.List(dates).get(ee.Number(date).add(-12)))
  var sum = pr.filterDate(date.advance(-120, "day"), date).sum().multiply(86400)
  
  var TC = sum.reduceRegions({
    collection: bufferedcents,
    reducer: ee.Reducer.mean(),
    scale: 2000,
    tileScale: 16
  })
  
  return ee.FeatureCollection(result).merge(TC)
}, dummy)


var out_hist = ee.List(seq_dates_hist).iterate(function(date, result){
  date = ee.Date(ee.List(dates_hist).get(date))
  var sum = pr_hist.filterDate(date.advance(-120, "day"), date).sum().multiply(86400)

  var TC = sum.reduceRegions({
    collection: bufferedcents,
    reducer: ee.Reducer.mean(),
    scale: 2000,
    tileScale: 16
  })
  
  return ee.FeatureCollection(result).merge(TC)
}, dummy)

//var all_pr = ee.Image(out_hist).addBands(ee.Image(out))
var all_pr = ee.FeatureCollection(out_hist).merge(ee.FeatureCollection(out))
print(ee.FeatureCollection(all_pr).limit(1))


//temperature
var out = ee.List(seq_dates).iterate(function(date, result){
  //var name = ee.String('temp_').cat(ee.String(ee.Number(date).toInt()))
  date = ee.Date(ee.List(dates).get(ee.Number(date).add(-12)))
  var sum = temp.filterDate(date.advance(-120, "day"), date).mean()
  //sum = sum.rename(ee.String(name))
  var TC = sum.reduceRegions({
    collection: bufferedcents,
    reducer: ee.Reducer.mean(),
    scale: 2000,
    tileScale: 16
  })
  
  return ee.FeatureCollection(result).merge(TC)
}, dummy)

var out_hist = ee.List(seq_dates_hist).iterate(function(date, result){
  date = ee.Date(ee.List(dates_hist).get(date))
  var sum = temp_hist.filterDate(date.advance(-120, "day"), date).mean()
  var TC = sum.reduceRegions({
    collection: bufferedcents,
    reducer: ee.Reducer.mean(),
    scale: 2000,
    tileScale: 16
  })
  
  return ee.FeatureCollection(result).merge(TC)
}, dummy)

var all_temp = ee.FeatureCollection(out_hist).merge(ee.FeatureCollection(out))
print(ee.FeatureCollection(all_temp).limit(1))


//export the ID and the mean distance-column to google drive
Export.table.toDrive({
    collection: all_temp,
    description: 'TC_cc',
    folder: 'Climate_Change',
    fileFormat: '',
    selectors: ["ID","mean"]
})

Export.table.toDrive({
    collection: all_pr,
    description: 'SC_cc',
    folder: 'Climate_Change',
    fileFormat: '',
    selectors: ["ID","mean"]
})