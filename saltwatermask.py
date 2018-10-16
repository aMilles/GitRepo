var hull = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-18.984375, 12.574585793638613],
          [-11.07421875, 3.007412867685379],
          [6.15234375, 1.4266161024061292],
          [9.66796875, -5.245580594057638],
          [9.4921875, -15.941605424270408],
          [18.10546875, -35.587038207116386],
          [29.70703125, -34.14484242914207],
          [42.36328125, -15.603280230491912],
          [40.95703125, -5.070511240121141],
          [51.328125, 4.482969197365171],
          [52.0751953125, 13.030962296093984],
          [47.9443359375, 12.001354739547983],
          [43.59375, 12.40296388105832],
          [38.49609375, 19.164546002664657],
          [35.33203125, 30.010767189066062],
          [11.953125, 37.73481526855993],
          [-6.328125, 36.04791923135713],
          [-18.984375, 23.260193648691715]]]);
var gsw = ee.ImageCollection('JRC/GSW1_0/MonthlyHistory').filterDate('2014-01-01', '2015-12-31');
var gsw_season = ee.ImageCollection('JRC/GSW1_0/MonthlyRecurrence');
var countries = ee.FeatureCollection('USDOS/LSIB/2013');


var africa = countries.filterMetadata('region', 'equals', 'AFRICA').union().first().geometry()

var africa = ee.Geometry(africa)

var hull = hull.difference(africa)
var africa = hull
var hull = hull.buffer(1500).simplify(1000)

Export.table.toAsset({
  collection: ee.FeatureCollection(africa),
  description: 'rough_africa'
})

Export.table.toAsset({
  collection: ee.FeatureCollection(hull),
  description: 'africa_hull5000'
})

var hull = ee.FeatureCollection('users/alexbenemilles/africa_hull5000')
var africa = ee.Feature(africa).difference(hull.first().geometry())


Map.addLayer(africa)

Export.table.toAsset({
  collection: ee.FeatureCollection(africa),
  description: 'africa_1500_iwbuffer'
})
