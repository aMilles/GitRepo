var segs = ee.FeatureCollection("users/alexbenemilles/segments_GEE"),
    africa = ee.FeatureCollection("users/alexbenemilles/africa_1500_iwbuffer");
    
var getCentroids = function(feature) {
  return feature.set({polyCent: feature.centroid()});
};

var cbuffer = function(feature) {
  return feature.buffer(5000);
};

var buffer_segs = segs.map(getCentroids).map(cbuffer)


var EVI = ee.ImageCollection('Oxford/MAP/EVI_5km_Monthly').filterDate('2010-01-01', '2014-12-31');

var empty = ee.Image().select()
var empty_IC = EVI.filterDate('2000-01-01', '2000-02-28')


// Do the same for Rain

var SC = ee.ImageCollection("NASA/GLDAS/V021/NOAH/G025/T3H").select('Rainf_tavg').filterDate('2009-06-30', '2014-12-31');

//var dates = ['2014-01-01', '2014-02-01', '2014-03-01', '2014-04-01', '2014-05-01', '2014-06-01', '2014-07-01', '2014-08-01', '2014-09-01', '2014-10-01', '2014-11-01', '2014-12-01'];
var dates = ["2010-01-01","2010-02-01","2010-03-01","2010-04-01","2010-05-01","2010-06-01","2010-07-01"
,"2010-08-01","2010-09-01","2010-10-01","2010-11-01","2010-12-01","2011-01-01","2011-02-01"
,"2011-03-01","2011-04-01","2011-05-01","2011-06-01","2011-07-01","2011-08-01","2011-09-01"
,"2011-10-01","2011-11-01","2011-12-01","2012-01-01","2012-02-01","2012-03-01","2012-04-01"
,"2012-05-01","2012-06-01","2012-07-01","2012-08-01","2012-09-01","2012-10-01","2012-11-01"
,"2012-12-01","2013-01-01","2013-02-01","2013-03-01","2013-04-01","2013-05-01","2013-06-01"
,"2013-07-01","2013-08-01","2013-09-01","2013-10-01","2013-11-01","2013-12-01","2014-01-01"
,"2014-02-01","2014-03-01","2014-04-01","2014-05-01","2014-06-01","2014-07-01","2014-08-01"
,"2014-09-01","2014-10-01","2014-11-01","2014-12-01"];
var dates = ee.List(dates)
var i = ee.List.sequence(0, 59, 1)

var months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]



var date = "2014-10-01";
var multiband = SC.filterDate(ee.Date(date).advance(-120, 'day'), ee.Date(date)).sum();//create sum of last 120 days 
  
var empty_FC = ee.Image(multiband).reduceRegions({
  collection: segs,
  reducer:ee.Reducer.mean(),
  scale: 1000,
  tileScale: 16
})
  
empty_FC = empty_FC.limit(0)


var SC_median = ee.List.sequence(0, 11, 1).iterate(function(iterator, result){
  var IC = ee.List.sequence(iterator, 59, 12).iterate(function(i, result){
    var date = dates.get(ee.Number(i).int());
    var multiband = SC.filterDate(ee.Date(date).advance(-120, 'day'), ee.Date(date)).sum(); //create sum of last 120 days 
    multiband = ee.Image(multiband).clip(africa);
    return ee.ImageCollection(result).merge(ee.ImageCollection(multiband))
  }, empty_IC)
  
  IC = ee.ImageCollection(IC).median();
  var name = ee.List(months).get(iterator)
  IC = IC.rename(ee.String(name));
  
  var SC_2014 = ee.Image(IC).reduceRegions({
    collection: buffer_segs,
    reducer:ee.Reducer.mean(),
    scale: 1000,
    tileScale: 16
  })
  
  return ee.FeatureCollection(result).merge(SC_2014)
}, empty_FC)

//VD

var VD_median = ee.List.sequence(0, 11, 1).iterate(function(iterator, result){
  var IC = ee.List.sequence(iterator, 59, 12).iterate(function(i, result){
    var date = dates.get(ee.Number(i).int());
    var multiband = EVI.filterDate(ee.Date(date), ee.Date(date).advance(15, 'day')).select('Mean').first();  
    multiband = ee.Image(multiband).clip(africa);
    return ee.ImageCollection(result).merge(ee.ImageCollection(multiband))
  }, empty_IC)
  
  IC = ee.ImageCollection(IC).median();
  var name = ee.List(months).get(iterator)
  IC = IC.rename(ee.String(name));
  
  var VD_2014 = ee.Image(IC).reduceRegions({
    collection: buffer_segs,
    reducer:ee.Reducer.mean(),
    scale: 1000,
    tileScale: 16
  })
  
  return ee.FeatureCollection(result).merge(VD_2014)
}, empty_FC)

print(ee.FeatureCollection(VD_median).limit(10))

//Distance to water

var gsw = ee.ImageCollection('JRC/GSW1_0/MonthlyHistory').filterDate('2010-01-01', '2014-12-31');
var gsw_season = ee.ImageCollection('JRC/GSW1_0/MonthlyRecurrence');
var africa = ee.FeatureCollection('users/alexbenemilles/africa_1500_iwbuffer')
var seg = ee.FeatureCollection('users/alexbenemilles/segments_GEE')
var empty = ee.Image().select();
var gsw_old = ee.ImageCollection('JRC/GSW1_0/MonthlyHistory').filterDate('2009-01-01', '2009-06-30');

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

//select the corresponding occurence bands and fill no data areas. 
//afterwards, compute the distance maps and return them as an image collection
//var empty = ee.ImageCollection(gsw_old);
//var WA = gsw.iterate(function(image, result){
  //save image for system:time_start extraction
//  var input = image //0 = no_data, 1 = no_water, 2 = water
  //get month and the respective seasonal occurence
//  var month = ee.Number(ee.Image(image).date().get('month')).subtract(1)
//  var fill = ee.Image(mb_gsw_season).select(month).multiply(2) //0 = no_occurence, 2 = occurence
  
//  var filled_image = image.eq(0).multiply(fill) //0 = no_occurence, 2 = occurence
  //add the new layer to the image where 
//  image = image.add(filled_image)
//  image = image.clip(africa)
//  var cost = ee.Image(gsw.first()).gte(-1)
//  var source = image.gte(2)
//  var dist = cost.cumulativeCost(source, 300000)
//  dist = ee.Image(dist).set('system:time_start', ee.Image(input).get('system:time_start'))
//  //return dist
//  return ee.Image(result).addBands(ee.Image(dist))
//}, empty)


var WA_median = ee.List.sequence(0, 11, 1).iterate(function(iterator, result){
     var season = ee.Image(mb_gsw_season).select(ee.Number(iterator)).multiply(2)
      
  var IC = ee.List.sequence(iterator, 59, 12).iterate(function(i, result){
    var date = ee.Date(ee.List(dates).get(i))
    var WA = ee.ImageCollection(gsw).filterDate(date, date.advance(5, "day")).first()
   //save image for system:time_start extraction
    //var input = WA //0 = no_data, 1 = no_water, 2 = water
    var filled_image = WA.eq(0).multiply(ee.Image(season)) //0 = no_occurence, 2 = occurence
    //add the new layer to the image where 
    WA = WA.add(filled_image)
    WA = WA.clip(africa)
    var cost = ee.Image(gsw.first()).gte(-1)
    var source = WA.gte(2)
    var dist = cost.cumulativeCost(source, 300000)
    //dist = ee.Image(dist).set('system:time_start', ee.Image(input).get('system:time_start'))
    return ee.ImageCollection(result).merge(ee.ImageCollection(dist))
  }, empty_IC)

  IC = ee.ImageCollection(IC).median();
  var name = ee.List(months).get(iterator)
  IC = IC.rename(ee.String(name));
  
  var WA_2014 = ee.Image(IC).reduceRegions({
    collection: segs,
    reducer:ee.Reducer.mean(),
    scale: 1000,
    tileScale: 16
  })
  
  return ee.FeatureCollection(result).merge(WA_2014)
}, empty_FC)



print(ee.FeatureCollection(WA_median))

//temp
var temp = ee.ImageCollection("NCEP_RE/surface_temp").filterDate('2010-01-01', '2014-12-31')

var TC_median = ee.List.sequence(0, 11, 1).iterate(function(iterator, result){
  var IC = ee.List.sequence(iterator, 59, 12).iterate(function(i, result){
    var date = dates.get(ee.Number(i).int());
    var multiband = temp.filterDate(ee.Date(date), ee.Date(date).advance(1, "month")).median();  
    multiband = ee.Image(multiband).clip(africa);
    return ee.ImageCollection(result).merge(ee.ImageCollection(multiband))
  }, empty_IC)
  
  IC = ee.ImageCollection(IC).median();
  var name = ee.List(months).get(iterator)
  IC = IC.rename(ee.String(name));
  
  var TC_2014 = ee.Image(IC).reduceRegions({
    collection: buffer_segs,
    reducer:ee.Reducer.mean(),
    scale: 1000,
    tileScale: 16
  })
  
  return ee.FeatureCollection(result).merge(TC_2014)
}, empty_FC)

print(ee.FeatureCollection(TC_median).limit(10))

//



//print(SC_2014.limit(10))
Export.table.toDrive({
    collection: SC_median,
    description: 'SC',
    folder: 'Seasonal Shift',
    fileFormat: '',
    selectors: ["ID", "mean"]
});

Export.table.toDrive({
    collection: VD_median,
    description: 'VD',
    folder: 'Seasonal Shift',
    fileFormat: '',
    selectors: ['ID', 'mean']
});

Export.table.toDrive({
    collection: WA_median,
    description: 'WA',
    folder: 'Seasonal Shift',
    fileFormat: '',
    selectors: ['ID', 'mean']
});

Export.table.toDrive({
    collection: TC_median,
    description: 'TC',
    folder: 'Seasonal Shift',
    fileFormat: '',
    selectors: ['ID', 'mean']
});