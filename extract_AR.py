var lc = ee.ImageCollection('MODIS/051/MCD12Q1').select('Land_Cover_Type_1');
var roads = ee.FeatureCollection('users/alexbenemilles/roads');
var seg = ee.FeatureCollection('users/alexbenemilles/segments_GEE')
var africa = ee.FeatureCollection('users/alexbenemilles/rough_africa')

var lc_2013 = ee.Image(lc.filterDate('2013-01-01', '2017-01-01').first());
var source = ee.Image().toByte().paint(roads, 1);

var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var bufferPoly = function(feature) {
  return feature.buffer(1500);   // substitute in your value of Z here
};

var centroids = seg.map(getCentroids);

var buffered_centroids = centroids.map(bufferPoly);


//how many minutes per kilometer of a certain surface?
var cost = 
  lc_2013.eq(0).multiply(1/0.10).add(
  lc_2013.eq(1).multiply(1/3.24).add(
  lc_2013.eq(2).multiply(1/1.62).add(
  lc_2013.eq(3).multiply(1/3.24).add(
  lc_2013.eq(4).multiply(1/4.00).add(
  lc_2013.eq(5).multiply(1/3.24).add(
  lc_2013.eq(6).multiply(1/3.00).add(
  lc_2013.eq(7).multiply(1/4.20).add(
  lc_2013.eq(8).multiply(1/4.86).add(  
  lc_2013.eq(9).multiply(1/4.86).add(
  lc_2013.eq(10).multiply(1/4.86).add(
  lc_2013.eq(11).multiply(1/1.00).add(
  lc_2013.eq(12).multiply(1/2.50).add(
  lc_2013.eq(13).multiply(1/5.00).add(
  lc_2013.eq(14).multiply(1/3.24).add(
  lc_2013.eq(15).multiply(1/1.62).add(
  lc_2013.eq(16).multiply(1/3.00).add(
  lc_2013.gt(16).multiply(1/3.00))))))))))))))))));

cost = cost.divide(1000) // as meters


var cost_map = cost.cumulativeCost({
  source: source,
  maxDistance: 250 * 1000
});

Map.addLayer(cost_map)

var PI = cost_map.reduceRegions({
  collection: buffered_centroids,
  reducer: ee.Reducer.mean(),
  scale: 500,
  tileScale: 16
});

print(PI.limit(10))

Export.table.toDrive({
  collection: PI,
  folder: 'GEC',
  description: 'PI',
  selectors: ['ID', 'mean']
});

//0	1c0dff	Water .10 ???
//1	05450a	Evergreen needleleaf forest 3.24
//2	086a10	Evergreen broadleaf forest 1.62
//3	54a708	Deciduous needleleaf forest 3.24
//4	78d203	Deciduous broadleaf forest 4.00
//5	009900	Mixed forest 3.24
//6	c6b044	Closed shrublands 3.00
//7	dcd159	Open shrublands 4.20
//8	dade48	Woody savannas 4.86
//9	fbff13	Savannas 4.86
//10	b6ff05	Grasslands 4.86
//11	27ff87	Permanent wetlands 1.00 ???
//12	c24f44	Croplands 2.5
//13	a5a5a5	Urban and built-up 5.00
//14	ff6d4c	Cropland/natural vegetation mosaic 3.24
//15	69fff8	Snow and ice 1.62
//16	f9ffa4	Barren or sparsely vegetated 3.00