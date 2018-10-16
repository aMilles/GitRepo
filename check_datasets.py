var seg = ee.FeatureCollection('users/alexbenemilles/segments_withinttime');
var VI = ee.ImageCollection('MODIS/006/MOD13Q1').select('SummaryQA');


var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var centroids = seg.map(getCentroids);

seg = centroids.map(function(feature) {
  return feature.buffer(5000);
});

var VI_filtered = VI.filterDate('2014-01-01', '2015-12-31')

var dummy_filter = VI.filter(ee.Filter.date('2000-01-01', '2000-12-31'))

//create an empty feature collection with properties for the output
var prep = dummy_filter.iterate(function(image, result){
  var start = ee.Number(image.get('system:time_start')).divide(1000000000000)
  var end = ee.Number(image.get('system:time_end')).divide(1000000000000)
  var features = seg.filter(ee.Filter.rangeContains('inttime', start, end))

  var out = image.reduceRegions({
    collection: features,
    reducer: ee.Reducer.max(),
    scale: 100,
    tileScale: 1
  })
  return out

})

//run the same function as for the dummy but us the filter from 2014-2015 and stack the consecutive results
var VD = VI_filtered.iterate(function(image, result){
  var start = ee.Number(image.get('system:time_start')).divide(1000000000000)
  var end = ee.Number(image.get('system:time_end')).divide(1000000000000)
  var features = seg.filter(ee.Filter.rangeContains('inttime', start, end))
  
  var out = image.reduceRegions({
    collection: features,
    reducer: ee.Reducer.min(),
    scale: 100,
    tileScale: 1
  })

  return ee.FeatureCollection(result).merge(out)
}, prep)


//export the ID and the mean distance-column to google drive
Export.table.toDrive({
    collection: VD,
    description: 'VD',
    folder: 'GEC_QC',
    fileFormat: '',
    selectors: ['ID', 'min']
})