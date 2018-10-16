var td = ee.ImageCollection('MODIS/051/MOD44B').select('Percent_Tree_Cover').filterDate('2014-01-01', '2016-01-01');
var seg = ee.FeatureCollection('users/alexbenemilles/segments_GEE');

var date = ee.Date('2015-01-01').millis().divide(1000000000000);

var seg_2015 = seg.filter(ee.Filter.gt('time', date));
var seg_2014 = seg.filter(ee.Filter.lte('time', date));
var td_2014 = td.filterDate('2014-01-01', '2014-12-01')
var td_2015 = td.filterDate('2015-01-01', '2015-12-01')
print(td_2014)

var first = ee.Image(td_2014.first()).reduceRegions({
  collection: seg_2014,
  tileScale: 1,
  reducer: ee.Reducer.mean(),
  scale: 100
})

var second = ee.Image(td_2015.first()).reduceRegions({
  collection: seg_2015,
  tileScale: 1,
  reducer: ee.Reducer.mean(),
  scale: 100
})
print(first)
var td_all = first.merge(second)
print(td_all)
Export.table.toDrive({
  collection: td_all,
  folder: 'GEC',
  description: 'TD',
  selectors: ['ID', 'mean'],
  fileFormat: ''
})
