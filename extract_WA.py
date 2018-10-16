var gsw = ee.ImageCollection('JRC/GSW1_0/MonthlyHistory').filterDate('2014-01-01', '2015-12-31');
var gsw_season = ee.ImageCollection('JRC/GSW1_0/MonthlyRecurrence');
var africa = ee.FeatureCollection('users/alexbenemilles/africa_1500_iwbuffer')
var seg = ee.FeatureCollection('users/alexbenemilles/segments_GEE')
var empty = ee.Image().select();
var gsw_old = ee.ImageCollection('JRC/GSW1_0/MonthlyHistory').filterDate('2013-01-01', '2013-06-30');

//Map.addLayer(africa)

//collapse the imagecollection of monthly occurene to a binary multiband image 
//1 = occurence rate > 50%
//0 = occurence rate < 50%
var mb_gsw_season = gsw_season.iterate(function(image, result){
  var name = ee.Number(image.get('system:index'))
  image = image.select('monthly_recurrence').rename(ee.String(name))
  image = image.gt(50).clip(africa)
  return ee.Image(result).addBands(ee.Image(image))
}, empty)

//Map.addLayer(ee.Image(mb_gsw_season).select('10'))
Map.addLayer(ee.ImageCollection(gsw).first())

//select the corresponding occurence bands and fill no data areas. 
//afterwards, compute the distance maps and return them as an image collection
var empty = ee.ImageCollection(gsw_old);
var filled_gsw = gsw.iterate(function(image, result){
  //save image for system:time_start extraction
  var input = image //0 = no_data, 1 = no_water, 2 = water
  //get month and the respective seasonal occurence
  var month = ee.Number(ee.Image(image).date().get('month')).subtract(1)
  var fill = ee.Image(mb_gsw_season).select(month).multiply(2) //0 = no_occurence, 2 = occurence
  
  var filled_image = image.eq(0).multiply(fill) //0 = no_occurence, 2 = occurence
  //add the new layer to the image where 
  image = image.add(filled_image)
  image = image.clip(africa)
  var cost = ee.Image(gsw.first()).gte(-1)
  var source = image.gte(2)
  var dist = cost.cumulativeCost(source, 1000000)
  dist = ee.Image(dist).set('system:time_start', ee.Image(input).get('system:time_start'))
  //return dist
  return ee.ImageCollection(result).merge(ee.ImageCollection(dist))
}, empty)


//remove dummy
filled_gsw = ee.ImageCollection(filled_gsw).filterDate('2014-01-01', '2015-12-31')

//print(ee.Image(ee.ImageCollection(filled_gsw).first()))
//print(filled_gsw)
//print(ee.ImageCollection(filled_gsw).filterDate('2014-01-01', '2015-01-01'))
//Map.addLayer(ee.Image(ee.ImageCollection(filled_gsw).first()))

var dummy = seg.filter(ee.Filter.lt('time', 0))
var dummy_ext = ee.Image(filled_gsw.first()).reduceRegions({
  collection: dummy,
  reducer: ee.Reducer.mean(),
  tileScale: 16,
  scale: 1000
})

var WA = filled_gsw.iterate(function(image, result){
  var date = image.date().getRange('month')
  var start = date.start().millis().divide(1000000000000)
  var end = date.end().millis().divide(1000000000000)
  var features = seg.filter(ee.Filter.rangeContains('time', start, end))
  var out = image.reduceRegions({
    collection: features,
    reducer: ee.Reducer.mean(),
    scale: 1000,
    tileScale: 16
  })
  return ee.FeatureCollection(result).merge(ee.FeatureCollection(out))
}, dummy_ext)


var start = ee.Date("2015-11-01").millis().divide(1000000000000)
var end = ee.Date("2015-12-31").millis().divide(1000000000000)
var seg_late = seg.filter(ee.Filter.rangeContains('time', start, end))
var last_gsw = ee.Image(filled_gsw.filterDate('2015-10-01', '2015-10-31').first())


var WA_add = last_gsw.reduceRegions({
  collection: seg_late,
  scale: 1000,
  tileScale: 16,
  reducer: ee.Reducer.mean()
})



WA = ee.FeatureCollection(WA).merge(ee.FeatureCollection(WA_add))


//print(ee.FeatureCollection(WA).limit(10))

Export.table.toDrive({
  collection: WA,
  selectors: ['ID', 'mean'],
  description: 'WA',
  folder: 'GEC'
})

//Map.addLayer(ee.Image(filled_gsw))
//print(ee.Image(filled_gsw))
Map.addLayer(ee.Image(ee.ImageCollection(filled_gsw).first()))