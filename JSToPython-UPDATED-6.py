#!/usr/bin/env python
# coding: utf-8

# In[24]:


import ee
import geemap


# In[25]:


ee.Initialize()


# In[26]:


table = ee.FeatureCollection("users/mazarderakhsh/apaLandClass202111") 


# In[27]:


lake = ee.FeatureCollection('projects/ee-mazard/assets/adk-samplelakes')


# In[28]:


# Create a Map object
Map = geemap.Map()


# In[29]:


# Map.addLayer(adb);


# In[30]:


Map


# In[31]:


station = ee.Geometry.Point(-75.09767330651451, 43.702075280163015)


# In[32]:


Map.addLayer(station)


# In[33]:


def maskL457sr(image):
    qaMask = image.select('QA_PIXEL').bitwiseAnd(int('11111', 2)).eq(0)
    saturationMask = image.select('QA_RADSAT').eq(0)
    qa = image.select('QA_PIXEL')
    cloud = qa.bitwiseAnd(1 << 5).Or(qa.bitwiseAnd(1 << 7)).Or(qa.bitwiseAnd(1 << 3))
    masked = image.updateMask(cloud.Not())
    opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermalBand = image.select('ST_B6').multiply(0.00341802).add(-124.15)
    return image.addBands(opticalBands, None, True).addBands(thermalBand, None, True).updateMask(masked).updateMask(qaMask).updateMask(saturationMask)



# In[34]:


def maskl7toa(image):
    # Bit 0 - Fill
    # Bit 1 - Dilated Cloud
    # Bit 2 - Unused
    # Bit 3 - Cloud
    # Bit 4 - Cloud Shadow
    qaMask = image.select('QA_PIXEL').bitwiseAnd(int('11111', 2)).eq(0)
    saturationMask = image.select('QA_RADSAT').eq(0)
    qa = image.select('QA_PIXEL')
    
    # Identify cloud pixels using the QA_PIXEL band
    cloud = qa.bitwiseAnd(1 << 5) \
                .Or(qa.bitwiseAnd(1 << 7)) \
                .Or(qa.bitwiseAnd(1 << 3))
    
    # Mask out cloud pixels
    masked = image.updateMask(cloud.Not())
    
    # Apply the scaling factors to the appropriate bands.
    opticalBands = image.select('B.')  #multiply(0.0000275).add(-0.2)
    thermalBand = image.select('B7')  #multiply(0.00341802).add(-124.15)
    
    # Replace the original bands with the scaled ones and apply the masks.
    return image.addBands(opticalBands, None, True) \
        .addBands(thermalBand, None, True) \
        .updateMask(masked) \
        .updateMask(qaMask) \
        .updateMask(saturationMask)


# In[35]:


def maskL8sr(image):
    # Bit 0 - Fill
    # Bit 1 - Dilated Cloud
    # Bit 2 - Cirrus
    # Bit 3 - Cloud
    # Bit 4 - Cloud Shadow
    qaMask = image.select('QA_PIXEL').bitwiseAnd(int('11111', 2)).eq(0)
    saturationMask = image.select('QA_RADSAT').eq(0)

    # Apply the scaling factors to the appropriate bands.
    opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermalBands = image.select('ST_B.*').multiply(0.00341802).add(-124.15)

    # Replace the original bands with the scaled ones and apply the masks.
    return image.addBands(opticalBands, None, True) \
        .addBands(thermalBands, None, True) \
        .updateMask(qaMask) \
        .updateMask(saturationMask)


# In[36]:


def cdom(image):
    co = image.expression("(20.3 - 10. * (b2 / b3) - 2.4 * (b3 / b4))", {
        'b1': image.select('blue'),
        'b2': image.select('green'),
        'b3': image.select('red'),
        'b4': image.select('nir')
    }).rename("CO")

    bad1 = co.where(co.gte(-2).And(co.lt(10)), 1).rename("bad1")
    bad2 = bad1.where(co.gte(10).Or(co.lt(-2)), 0).rename("bad2")
    cdo = co.multiply(bad2).rename("cdom")
    mask = cdo.neq(0)

    return ee.Image(image.addBands([co, cdo]).clip(lake).updateMask(mask))


# In[43]:


# // Map the function over one year of data.
collection2 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filter(ee.Filter.calendarRange(9,11,'month')).filterBounds(lake).map(maskL8sr);
composite = collection2.median().clip(lake);
# // Display the results.
Map.centerObject(lake, 16);    

 
   

# Assign a common name to the sensor-specific bands.
LC8_BANDS = ['SR_B2','SR_B3','SR_B4',  'SR_B5','SR_B6','SR_B7','ST_B10','QA_PIXEL']; #//Landsat 8
LC7_BANDS = ['SR_B1',   'SR_B2',    'SR_B3',  'SR_B4',  'SR_B5',    'SR_B7',  'ST_B6','QA_PIXEL']; #//Landsat 7
LC7T_BANDS = ['B1', 'B2', 'B3', 'B4','B5','B7', 'B8','QA_PIXEL']; #//Landsat 7
LC5_BANDS = ['SR_B1',   'SR_B2',    'SR_B3',  'SR_B4',  'SR_B5',    'SR_B7',    'ST_B6','QA_PIXEL']; #//Llandsat 5
STD_NAMES = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'temp','QA'];


l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(lake) \
.filter(ee.Filter.lt('CLOUD_COVER', 25)) \
.map(maskL8sr).filter(ee.Filter.calendarRange(5,11,'month')) \
.select(LC8_BANDS, STD_NAMES).map(cdom); #// Landsat 8



l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(lake) \
.filter(ee.Filter.lt('CLOUD_COVER', 35)) \
.filter(ee.Filter.calendarRange(5,11,'month')).map(maskL457sr). \
select(LC7_BANDS, STD_NAMES).map(cdom); #// Landsat 7

l7toa = ee.ImageCollection('LANDSAT/LE07/C02/T1_TOA').filterBounds(lake) \
.filter(ee.Filter.lt('CLOUD_COVER', 25)).map(maskl7toa) \
.filter(ee.Filter.calendarRange(5,11,'month')) \
.select(LC7T_BANDS, STD_NAMES) \
.map(cdom);

l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(lake) \
.filter(ee.Filter.lt('CLOUD_COVER', 45)) \
.filter(ee.Filter.calendarRange(7,11,'month')).map(maskL457sr). \
select(LC5_BANDS, STD_NAMES).map(cdom); #// Landsat 5

print(l7toa,'l7t');
# //var lall = ee.ImageCollection(l5.merge(l7).merge(l8));

lnew7t = l7toa;
lnew7=l5;


print(lnew7,'lnew7' );

Map.addLayer(l5.filter(ee.Filter.calendarRange(1985,1987,'year')).filter(ee.Filter.calendarRange(6,7,'month')).mean(), 
             {'bands': ['red', 'green', 'blue'], min: 0, max: 0.1},'l7-2000');
Map.addLayer(l5.filter(ee.Filter.calendarRange(2009,2011,'year')).filter(ee.Filter.calendarRange(6,7,'month')).mean().clip(lake), 
             {'bands': ['red', 'green', 'blue'], min: 0, max: 0.1},'l7-2020');


# // Define a region to calculate histogram for.
histRegion = lake;

# // Define the chart and print it to the console.
data1 = l5.select(['cdom']).filter(ee.Filter.calendarRange(1985,1986,'year')).mean().clip(lake);
data2 = l5.select(['cdom']).filter(ee.Filter.calendarRange(2010,2011,'year')).mean().clip(lake);

chart2 = ui.chart.Image.histogram({image: [data1, data2], region: histRegion,minBucketWidth: 0.1,scale: 30})\
    .setSeriesNames(['CDOM-1985', 'CDOM-2010'])\
    .setOptions({
        title: 'Histogram of Landsat-5 estimated CDOM for sample lakes',
        hAxis: {
            title: 'CDOM (1/m)',
            titleTextStyle: {italic: False, bold: True},
        },
        vAxis: {
            title: 'Count',
            titleTextStyle: {italic: False, bold: True},
        },
        colors: ['1d6b99', 'magenta'],
    })

# Print the chart
print(chart2)


# In[45]:


colle = l5.select("cdom");
collecr = colle.reduce(ee.Reducer.median()).clip(lake);
#//Map.addLayer(collecr.select("cdom.*"),{palette: ['231fff','418dff','14a8ff','d9ff08','ffa808','ff3406'], min:0, max:10}, 'water color');

col1 = l5.filter(ee.Filter.calendarRange(1985,1987,'year')).filter(ee.Filter.calendarRange(6,8,'month')).qualityMosaic('cdom'); #//.reduce(ee.Reducer.median()).clip(lake);
print(col1,'col1');

col2 = l5.select("cdom").filter(ee.Filter.calendarRange(2009,2011,'year')).filter(ee.Filter.calendarRange(6,8,'month')).qualityMosaic('cdom').clip(lake);

#//Map.addLayer(col1.select("cdom.*"),{palette: ['231fff','418dff','14a8ff','d9ff08','ffa808','ff3406'], min:0, max:7}, 'water color-2000');

Map.addLayer(col1.select("cdom"),{'palette': ['0049ed','379ae0','b9fbf5','eaf71e','ffaf1d','916d5f'], min:0.2, max:3.0}, 'water color-2000');
Map.addLayer(col2.select("cdom"),{'palette': ['0049ed','379ae0','b9fbf5','eaf71e','ffaf1d','916d5f'], min:0.2, max:3.0}, 'water color-2020');

colordiff = col2.subtract(col1);


#// Define a boxcar or low-pass kernel.
boxcar = ee.Kernel.square({
  radius: 5, units: 'pixels', normalize: true
});

#// Smooth the image by convolving with the boxcar kernel.

smooth = colordiff.convolve(boxcar);

#// Get a palette: a list of hex strings

palettes = require('users/gena/packages:palettes');
palettem = palettes.colorbrewer.RdYlBu[11].reverse();

vizParams = {
  min: -.5,
  max: .5,
 #// gamma: [1, 0, -.5]
palette : ['313695','4575b4','74add1','abd9e9','e0f3f8','ffffff','fefefe','fdae61','f46d43','d73027','a50026']
};
#//Map.addLayer(colordiff.select("cdom.*"),{palette: ['231fff','418dff','c6fffc','fff9af','ffa808','ff3406'], min:-1, max:4}, 'water color difference');

Map.addLayer(smooth.select("cdom.*"),vizParams, 'water color difference');


# In[59]:


chart = ui.Chart.image.series({
    imageCollection: lnew.select(['cdom']),
    region: station,
    reducer: ee.Reducer.mean(),
    scale: 30
     }).setOptions({
    title: 'Mean Surface Reflectance Value by Date for Bear Pond',
    hAxis: {title: 'Date', titleTextStyle: {italic: false, bold: true}},
    vAxis: {title: 'Surface Reflectance',titleTextStyle: {italic: false, bold: true}},})
  
    
print(chart)

chartstd = ui.Chart.image.series({
    imageCollection: lnew.select('cdom'),
    region: station,
    reducer: ee.Reducer.stdDev(),
    scale: 30
     }).setOptions({
    title: 'Mean Surface Reflectance Value by Date for Little Clear Pond',
    hAxis: {title: 'Date', titleTextStyle: {italic: false, bold: true}},
    vAxis: {title: 'water CI',titleTextStyle: {italic: false, bold: true}},
  })
  
print(chartstd)




# ////////////////// Slope map 

imgcoll1 = l5.select('cdom');
sequence = ee.List.sequence(1985, 2011);
statsCollection = ee.ImageCollection(sequence.map(lambda year: imgcoll1.filter(ee.Filter.calendarRange(year, year, 'year')).first()))
 
mean = oneYear.mean().rename('mean');
yr = ee.Image(ee.Number(year)).float().rename('year')
return ee.Image.cat(mean,yr).set('year', year)


#// Linear fit

linFit = statsCollection.select(['year', 'mean']).reduce(ee.Reducer.linearFit()).rename(['slope', 'offset'])


smooth2 = linFit.select(['slope']).convolve(boxcar)

Map.addLayer(smooth2.clip(lake), {bands: ['slope'], palette: ['100cff','ccefff','e9fcff','ffffff','c6ce4f','a15e0a'],
                    min:-0.02,
                    max:[ 0.02,]}, "Slope of linear fit")

lnew=l5.select(['blue','green','red','nir','cdom']);
lnew2=l8.select(['blue','green','red','nir','cdom']);

data1 = lnew.filter(ee.Filter.calendarRange(1985,1986,'year')).mean().clip(station)
data2 = lnew.filter(ee.Filter.calendarRange(1990,1991,'year')).mean().clip(station)
data3 = lnew.filter(ee.Filter.calendarRange(2000,2001,'year')).mean().clip(station)
data4 = lnew.filter(ee.Filter.calendarRange(2005,2006,'year')).mean().clip(station)
data5 = lnew.filter(ee.Filter.calendarRange(2010,2011,'year')).mean().clip(station)
data=[data1,data2,data3,data4,data5]

print(data,'data')
 
#// ui.Chart.image.byRegion(image, regions, reducer, scale, xProperty)
chartx = ui.Chart.image.byRegion({
    image: data,
    regions: station,
    reducer: ee.Reducer.mean(),
    scale: 30
     }).setOptions({
          title: 'Mean Surface Reflectance Value by Date for Big Moose Lake',
          hAxis: {title: 'Bands', titleTextStyle: {italic: false, bold: true}},
          vAxis: {title: 'Surface Reflectance',titleTextStyle: {italic: false, bold: true}},
  });
print(chartx)


# In[ ]:





# In[ ]:




