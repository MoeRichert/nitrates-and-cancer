var map;

// Credits for map located on bottom right of map


// Basemap options located on top right of map
var grayscale   = L.tileLayer('https://tiles.stadiamaps.com/tiles/alidade_smooth/{z}/{x}/{y}{r}.png', {
	maxZoom: 20,
	attribution: '&copy; <a href="https://stadiamaps.com/">Stadia Maps</a>, &copy; <a href="https://openmaptiles.org/">OpenMapTiles</a> &copy; <a href="http://openstreetmap.org">OpenStreetMap</a> contributors'
}),
    dark  = L.tileLayer('https://tiles.stadiamaps.com/tiles/alidade_smooth_dark/{z}/{x}/{y}{r}.png', {
	maxZoom: 20,
	attribution: '&copy; <a href="https://stadiamaps.com/">Stadia Maps</a>, &copy; <a href="https://openmaptiles.org/">OpenMapTiles</a> &copy; <a href="http://openstreetmap.org">OpenStreetMap</a> contributors'
}),
    outdoors = L.tileLayer('https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png', {
	maxZoom: 17,
	attribution: 'Map data: &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, <a href="http://viewfinderpanoramas.org">SRTM</a> | Map style: &copy; <a href="https://opentopomap.org">OpenTopoMap</a> (<a href="https://creativecommons.org/licenses/by-sa/3.0/">CC-BY-SA</a>)'
});

//Legend div in bottom right of map
var legend = new (L.Control.extend({
     options: { position: 'bottomright' }
    }));

var legend2 = new (L.Control.extend({
     options: { position: 'bottomright' }
    }));


// Initialize global variables for data layers
var censusTracts,
    wellPoints,
    censusGrid,
    wellPointsTurf,
    hexGrid,
    cancertract_hexgrid,
    cancerTurf = [],
    cancertract_grid = [],
    regression_hexgrid,
    Weight = 6,
    binArea = 2;

//Layer options located in key
var wellLayer = L.layerGroup(),
    censusLayer = L.layerGroup(),
    nitrateIDWLayer = L.layerGroup(),
    regressionLayer = L.layerGroup(),
    idwLayer = L.layerGroup();

//style functions

//break Arrays
var cBreaks = [],
    wBreaks = [],
    ciBreaks = [],
    nBreaks = [],
    rBreaks = [];

//color scheme arrays
var orangeScheme = ['#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f', '#990000'],
    blueScheme = ['#eff4ff', '#c6dbef', '#9ecae1', '#6baed6', '#3182bd', '#08519c'],
    purpleScheme = ['#bfd3e6', '#9ebcda', '#8c96c6', '#8c6bb1', '#88419d', '#6e016b'],
    greenredScheme = ['#7fcdbb', '#c7e9b4','#ffffcc','#fdae61', '#ff3c3c', '#990000'],
    regressionColors = ['#6600a1', '#00a3fe', '#00eb10', '#00eb10', '#00a3fe', '#6600a1'];

//default style

function getcBreaks(data) {

    // Create an empty array to store breaks
    var values = [];

    // Loop through each feature to get its cancer rate
    data.eachLayer(function (layer) {
        var value = layer.feature.properties.canrate;

        // Push each cancer rate into the array
        values.push(value);
    });

    // Determine 5 clusters of statistically similar values, sorted in ascending order
    var clusters = ss.ckmeans(values, 5);

    // Create a 2-dimensional array of the break points (lowest and highest values) in each cluster. The lowest value in each cluster is cluster[0]; the highest value is cluster.pop().
    var breaks = clusters.map(function (cluster) {
        return [cluster.pop()];
    });
    
    var lowestvalCalc = clusters.map(function (cluster) {
        return [cluster[0]];
    });
    
    var lowestval = lowestvalCalc[0];
    
    breaks.unshift(lowestval);

    // Return the array of class breaks
    return breaks;

}; // end getcBreaks()
function getnBreaks(data) {

    // Create an empty array to store breaks
    var values = [];

    // Loop through each feature to get its cancer rate
    data.eachLayer(function (layer) {
        var value = layer.feature.properties.nitr_ran;

        // Push each cancer rate into the array
        values.push(value);
    });

    // Determine 5 clusters of statistically similar values, sorted in ascending order
    var clusters = ss.ckmeans(values, 5);

    // Create a 2-dimensional array of the break points (lowest and highest values) in each cluster. The lowest value in each cluster is cluster[0]; the highest value is cluster.pop().
    var breaks = clusters.map(function (cluster) {
        return [cluster.pop()];
    });
    
    var lowestvalCalc = clusters.map(function (cluster) {
        return [cluster[0]];
    });
    
    var lowestval = lowestvalCalc[0];
    
    breaks.unshift(lowestval);

    // Return the array of class breaks
    return breaks;

}; // end getnBreaks() 
function getrBreaks(data) {
    console.log(data);

    var values = [];
	data.eachLayer(function(layer) {
		var value = layer.feature.properties.residual;
		values.push(value);
	});
	var standardDeviation = ss.sampleStandardDeviation(values);
    
    // Create an array of the break points for -2, -1, 0, 1, and 2 standard deviations
    var breaks = [-2 * standardDeviation, -1 * standardDeviation, 0, standardDeviation, 2 * standardDeviation];
    
    console.log(standardDeviation);
    console.log(breaks);
	return breaks;

}; // end getrBreaks() 

function getColor(b, breakArray, colorArray) {
    return b >=breakArray[0]-1 & b < breakArray[1] ? colorArray[0]:
        b >= breakArray[1] & b < breakArray[2] ? colorArray[1]:
        b >= breakArray[2] & b < breakArray[3] ? colorArray[2]:
        b >= breakArray[3] & b < breakArray[4] ? colorArray[3]:
        b >= breakArray[4] & b < breakArray[5] ? colorArray[4]:
        b >= breakArray[5] ? colorArray[5]:
        '#ffffff';
    
};

//Slider UI Update Values

function updatehexUI(valArea) {
    var areaVal = document.getElementById("areaVal");
    areaVal.innerHTML = valArea + " KM";
    binArea = parseInt(valArea);
};
function updatedecayUI(valDecay) {
    var decayVal = document.getElementById("decayVal");
    decayVal.innerHTML = valDecay;
    Weight = parseInt(valDecay);
};
    

//createMap builds map, returns variable globalMap
function createMap(){
    
	//create map
    map = L.map('map', {
            center: [43.78, -88.78],
            zoom: 7,
            minZoom: 4,
            maxZoom: 14,
            layers:grayscale
        });
    
    // Add the Leaflet easyPrint plugin to the map
    var printer = L.easyPrint({
      		sizeModes: ['Current', 'A4Landscape', 'A4Portrait'],
      		filename: 'myMap',
      		exportOnly: true,
      		hideControlContainer: true
		}).addTo(map);
    
    var baseLayers = {
        "Grayscale": grayscale,
       "Topographic": outdoors,
       "Darkscale": dark,
    };
    
   var displayLayers = {
       "Cancer Rate": censusLayer,
       "Well Layer": wellLayer,
       "Nitrate Concentration": nitrateIDWLayer,
       "Interpolated Cancer Rate": idwLayer, 
       "Regression Residuals": regressionLayer,
   };
    
    
    //call getData
    getData(map);

    L.control.layers(baseLayers, displayLayers).addTo(map);
    
};


//getData loads geoJSON

function getData(map){
    
    // load the cancer tract data 
    var ctd = $.getJSON("data/cancer_tracts.geojson"),
        wn = $.getJSON("data/well_nitrate.geojson"),
        cd = $.getJSON("data/cancergrid5.geojson");
    
    //call function responses
    $.when(ctd, wn, cd).then(function (response1, response2, response3) {
     //You have both responses at this point.
       
        // create layer and add to the layer group
        censusTracts = L.geoJson(response1, {
            style: style,
            onEachFeature: onEachFeature
        }).addTo(censusLayer);

        cBreaks = getcBreaks(censusTracts);

        //default layer style
        function style(feature) {
            return {
                weight: 0.7,
                opacity: 0.7,
                color: 'white',
                fillOpacity: 0.7,
                fillColor: getColor(feature.properties.canrate, cBreaks, greenredScheme)
            };
        };    


        //iterate through each feature in the geoJSON
        function onEachFeature(feature, layer) {
                layer.on({
                    mouseover: highlightFeature,
                    mouseout: resetHighlight
                });
            };

        //set the highlight on the map
        function highlightFeature(e) {

            var layer = e.target;
            layer.setStyle({fillOpacity: 1.0,
                           color: 'red'});

            var popup = "<b>Cancer Rate: </b>" + (layer.feature.properties.canrate * 100).toLocaleString() + "%";
            layer.bindPopup(popup);
        };

        function resetHighlight(e) {
            censusTracts.setStyle(style);
        };
        
        var censusLegend,
            headingcensus = "Cancer Rate <br> % per Census Unit",
            unitcensus = "%";
        
        //draw legend for the census tracts
        censusLegend = percentLegends(cBreaks, greenredScheme, headingcensus, unitcensus);

        // create layer and add to layer group
        wellPoints = L.geoJson(response2, {
            pointToLayer: createCircleMarker
        }).addTo(wellLayer);
        
        wBreaks = getnBreaks(wellPoints);

        // Create a style for the well points
        function createCircleMarker (feature, latlng) {

            return L.circleMarker(latlng, {
                fillOpacity: 1,
                color: '#3d3d3d',
                weight: 0.25,
                opacity: 1,
                radius: 2.5
            });
        }
        
        wellPoints.eachLayer(function(layer) {
            layer.setStyle({
                fillColor: getColor(layer.feature.properties.nitr_ran, wBreaks, purpleScheme)
            });
            var popup = "<b>Nitrate Concentration: </b>" + layer.feature.properties.nitr_ran.toFixed(2) + " ppm";
            layer.bindPopup(popup);
        });
        
        var wellLegend,
            headingwell = "Well Nitrate Concentration",
            unitwell = "ppm";
        
        //draw legend for the well points
        wellLegend = createLegends(wBreaks, purpleScheme, headingwell, unitwell);
        
        // create grid layer for interpolation
        censusGrid = L.geoJson(response3, {
            style: style,
            onEachFeature: onEachFeature
        });
    
        
        censusLayer.addTo(map);
        wellLayer.addTo(map);
        
});
    
    
    
};

function interpolateData (ddC, hexA){
    
    //remove old data and clean up display
    
    nitrateIDWLayer.clearLayers();
    regressionLayer.clearLayers();
    idwLayer.clearLayers();
    nitrateIDWLayer.remove();
    censusLayer.remove();
    legend.remove();
    legend2.remove();
    wellLayer.remove();
    idwLayer.remove();
    regressionLayer.remove();
    
    //get q value
    var sliderValue = ddC,
        hexbinArea = hexA;
    
    var wellPointsArray = [];
    
    // Loop through each feature
    wellPoints.eachLayer(function (layer) {

        // Build a Turf feature collection from the well points

        // Create shorthand variables to access the layer properties and coordinates
        var prop = layer.feature.properties;
        var coords = layer.feature.geometry.coordinates;

        // Create a Turf point feature for the well point, with its coordinates and attributes
        wellPointsFeature = turf.point(coords, prop);

        // Push the current well point feature into an array
        wellPointsArray.push(wellPointsFeature);

    });

    // Create a Turf feature collection from the array of well point features
    wellPointsTurf = turf.featureCollection(wellPointsArray);

    // Set options for the well point interpolation
    var optionsWP = {
        gridType: 'hex', // use hexbins as the grid type
        property: 'nitr_ran', // interpolate values from the nitrate concentrations
        units: 'kilometers', // hexbin size units
        weight: sliderValue // distance decay coefficient, q
    };

    // Interpolate the well point features using the grid size from the hexbinArea variable, the submitted distance decay coefficient, and the options just specified
    nitrateRatesHexbinsTurf = turf.interpolate(wellPointsTurf, hexbinArea, optionsWP);
    var interpolatedNitrate = [];

    // Loop through each hexbin and get its interpolated nitrate concentration
    for (var hexbin in nitrateRatesHexbinsTurf.features) {
        var interpolatedNitrateRate = nitrateRatesHexbinsTurf.features[hexbin].properties.nitr_ran;
        interpolatedNitrate.push(interpolatedNitrateRate);
    }
    // Convert the hexbins to a Leaflet GeoJson layer and add it to the Nitrate Concentrations layer group
    nitrateRatesHexbins = L.geoJson(nitrateRatesHexbinsTurf, {

        // Style the nitrate concentration hexbins
        onEachFeature: function(feature, layer){
            layer.bindPopup("<b>Nitrate Concentration: </b>" + layer.feature.properties.nitr_ran.toFixed(2) + " ppm")
        }
    }).addTo(nitrateIDWLayer);

    // Get the class breaks based on the ckmeans classification method
    nBreaks = getnBreaks(nitrateRatesHexbins);
    
    function styleN(feature) {
            return {
                color: '#585858', 
                weight: 0.5, 
                fillOpacity: 0.6, 
                opacity: 0.5,
                fillColor: getColor(feature.properties.nitr_ran, nBreaks, blueScheme)
            };
    };
    
    nitrateRatesHexbins.setStyle(styleN);
    
    var wellLegend,
        headingwell = "Interpolated Nitrate Concentration",
        unitwell = "ppm";
        
    //draw legend for the nitrate layer
    wellLegend = createLegends(nBreaks, blueScheme, headingwell, unitwell);



    
    ////////////////////////census tracts interpolation
    
    var censusArray = [];
    
    // Loop through each census tract feature and build a Turf feature collection from its centroid
    censusGrid.eachLayer(function (layer) {

        // Create shorthand variables to access the layer properties and coordinates
        var props = layer.feature.properties;
        var coords = layer.feature.geometry.coordinates;

        // Get the centroid of the census tract
        var censusPts = turf.point(coords, props);

        // Push the current census tract centroid into an array
        censusArray.push(censusPts);

    });

    // Create a Turf feature collection from the array of census tract centroid features
    censusTurf = turf.featureCollection(censusArray);
    
    var options = {gridType: 'hex', property: 'canrate', units: 'kilometers'};
    cancerTurf = turf.interpolate(censusTurf, hexbinArea, options);
    
    var gridToPoint = [];
    
    hexGrid = new L.GeoJSON(cancerTurf, {
        onEachFeature: function(feature, layer) {
            layer.bindPopup("<b>Interpolated Cancer Rate: </b>" + (layer.feature.properties.canrate * 100).toFixed(1) + "%");
            gridToPoint.push(turf.centroid(feature, feature.properties));
        }
        }).addTo(idwLayer);
        
    ciBreaks = getcBreaks(hexGrid);
    
    //default layer style
    function style(feature) {
        return {
            weight: 0.7,
            opacity: 0.7,
            color: 'white',
            fillOpacity: 0.7,
            fillColor: getColor(feature.properties.canrate, ciBreaks, orangeScheme)
        };
    };
    
    
    hexGrid.setStyle(style);
    
    var cancerLegend,
        headingcan = "Interpolated Cancer Percent",
        unitcan = "%";
        
    //draw legend for the nitrate layer
    cancerLegend = percentLegends(ciBreaks, orangeScheme, headingcan, unitcan);
    
    
    nitrateIDWLayer.addTo(map);
    
    
};

function joinFeatures(distanceDecayCoefficient, hexbinArea) {

    // Set options for the cancer rate interpolation by grid points
    var gridOptions = {
        gridType: 'point', // use points as the grid type, required to use the collect function
        property: 'canrate', // interpolate values from the cancer rates
        units: 'kilometers', // hexbin size units
        weight: distanceDecayCoefficient // distance decay coefficient, q
    };

    // Interpolate the cancer rate centroids into a surface of grid points (http://turfjs.org/docs#interpolate)
    cancerRatesGridPointsTurf = turf.interpolate(censusTurf, hexbinArea, gridOptions);

    // Use the collect function to join the cancer rates from the cancer rate grid points to the nitrate concentration hexbins (http://turfjs.org/docs/#collect)
    collectedFeaturesHexbinsTurf = turf.collect(nitrateRatesHexbinsTurf, cancerRatesGridPointsTurf, 'canrate', 'values');

    // Loop through each of the collected hexbins
    for (var i in collectedFeaturesHexbinsTurf.features) {

        // The collect function builds an array of cancer rates for features intersecting the current hexbin
        // Get the array of cancer rates for the current hexbin
        var canrateArray = collectedFeaturesHexbinsTurf.features[i].properties.values;

        // Loop through each feature in the cancer rates array and sum them
        var canrateArraySum = 0;
        for (var j in canrateArray) {

            if (canrateArray.length > 0) {
                canrateArraySum += parseFloat(canrateArray[j]);
            }

        }

        // Get the average cancer rate (sum / number of features in the array)
        var canrateArrayAvg = canrateArraySum / canrateArray.length;

        // Add the average cancer rate to the canrate property of the current hexbin
        if (canrateArrayAvg !== undefined) {
            collectedFeaturesHexbinsTurf.features[i].properties.canrate = canrateArrayAvg;
        } else {
            collectedFeaturesHexbinsTurf.features[i].properties.canrate = "";
        }

    }

    // Call the function to calculate linear regression using the joined hexbins
    calculateLinearRegression(collectedFeaturesHexbinsTurf);

};

function calculateLinearRegression(dataArray){
    
    censusLayer.remove();
    wellLayer.remove();
    nitrateIDWLayer.remove();
    idwLayer.remove();
    legend.remove();
    legend2.remove();
    
    var interpArray = [],
        obsArray = [];
    
    
    // Loop through the hexbin layer with nitrate concentrations and cancer rates
    // Create a two-dimensional array of [x, y] pairs where x is the nitrate concentration and y is the cancer rate

    // Loop through each of the collected hexbins
    for (var i in dataArray.features) {

        // Create a shorthand variable to access the layer properties
        var props = dataArray.features[i].properties;

        // Create variables to store the interpolated nitrate concentration and cancer rate
        var interpolatedNitrateConcentration = props.nitr_ran;
        var interpolatedCancerRate = props.canrate;

        // Create an array for the current feature of [nitrate concentration, cancer rate]
        var currentNitrateAndCancerRates = [parseFloat(interpolatedNitrateConcentration), parseFloat(interpolatedCancerRate)];

        // Push the array of the current feature's nitrate concentration and cancer rate into an array
        interpArray.push(currentNitrateAndCancerRates);

    }

    // Run the linearRegression method from the Simple Statistics library to return an object containing the slope and intercept of the linear regression line
    // where nitrate concentration is the independent variable (x) and cancer rate is the dependent variable (y)
    // The object returns m (slope) and b (y-intercept) that can be used to predict cancer rates (y) using the equation, y = mx + b
    var regressionEquation = ss.linearRegression(interpArray);

    // Create variables for the slope and y-intercept
    var m = regressionEquation.m;
    var b = regressionEquation.b;
    
    var regressionEquationFormatted = "y = " + parseFloat(m).toFixed(5) + "x + " + parseFloat(b).toFixed(5);
    console.log("Regression Equation: " + regressionEquationFormatted);
    
    var divequation = document.getElementById("equation")
    divequation.innerHTML = "<b>Regression Equation: </b><br>" + regressionEquationFormatted;
    

    // Loop through each of the collected hexbins
    for (var j in dataArray.features) {

        // calculate the predicted cancer rate from the interpolated nitrate concentration
        var predictedCancerRate = m * (parseFloat(dataArray.features[j].properties.nitr_ran)) + b;

        // Calculate the residual
        var residual = predictedCancerRate - dataArray.features[j].properties.canrate;

        // Add the predicted cancer rate and residual to the hexbin
        dataArray.features[j].properties.predictedCancerRate = predictedCancerRate;
        dataArray.features[j].properties.residual = residual;
        
        // Build an array of the observed nitrate concentration and cancer rate for the current feature
        var observedNitrateAndCancerRatesPair = [dataArray.features[j].properties.nitr_ran, dataArray.features[j].properties.canrate];
        
        // Push the current nitrate concentration and cancer rate pair into an array
        obsArray.push(observedNitrateAndCancerRatesPair);

    }
    
    // Calculate the r-squared for the regression
    
    // Build the linear regression line
    var regressionLine = ss.linearRegressionLine(regressionEquation);
    
    // Calculate the r-squared
    var rSquared = parseFloat(ss.rSquared(obsArray, regressionLine)).toFixed(5); // 1 is a perfect fit, 0 indicates no correlation
    console.log("r-Squared: " + rSquared);
    
    // Convert the collected hexbins to a Leaflet GeoJson layer
    regression_hexbins = L.geoJson(dataArray, {

        // Set a default style for the collected hexbins
        style: function (feature) {
            return {
                color: '#999999', // Stroke Color
                weight: 0.5, // Stroke Weight
                fillOpacity: 0.5, // Override the default fill opacity
                opacity: 0.5 // Border opacity
            };
        }

    }).addTo(regressionLayer);

    // Get the class breaks based on the ckmeans classification method
    var breaks = getrBreaks(regression_hexbins);

    // Loop through each feature, set its symbology, and build and bind its popup
    regression_hexbins.eachLayer(function (layer) {

        // Set its color based on the residual between the predicted and observed cancer rate
        layer.setStyle({
            fillColor: getColor(layer.feature.properties.residual, breaks, regressionColors)
        });

        // Build the popup for the feature
        var popup = "<b>Nitrate Concentration: </b>" + layer.feature.properties.nitr_ran.toFixed(2) + " ppm" + "<br/>" +
            "<b>Observed Cancer Rate: </b>" + (layer.feature.properties.canrate * 100).toFixed(2).toLocaleString() + "% of census tract population" + "<br/>" +
            "<b>Predicted Cancer Rate: </b>" + (layer.feature.properties.predictedCancerRate * 100).toFixed(2).toLocaleString() + "% of census tract population";
        layer.bindPopup(popup, {
			className: 'my-popup',
			closeButton: false
		});
		//hovering options
		layer.on('mouseover', function(e) {
			this.openPopup();
		});
		layer.on('mouseout', function(e) {
			this.closePopup();
		});

    });
    
    var regresLegend,
        headingreg = "Regression Residuals",
        unitres = "Std. Dev.",
        regBreaks = [-2, -1, 0, 1, 2];
        
    //draw legend for the nitrate layer
    regBregresLegend = regLegends(regBreaks, regressionColors, headingreg, unitres);
      
    regressionLayer.addTo(map);
};

function createLegends(breaks, colors, header, unit){
    
    legend.onAdd = function (map) {      
      var div = L.DomUtil.create('div', 'info legend'),
          labels = ['<strong>' + header + '</strong>'];
        
        //iterate through grades and create a color field and label for each
        for (var i = 0; i < breaks.length; i++) {
            labels.push(
                (breaks[i + 1]) ? '<i style="background:' + getColor(breaks[i+1]-.01, breaks, colors) + '"></i> ' +
                (parseFloat(breaks[i]).toFixed(2).toLocaleString() + '&ndash;' + parseFloat(breaks[i + 1]).toFixed(2).toLocaleString()) + ' ' + unit : '');
        }
        
        
        div.innerHTML = labels.join('<br>');
        return div;
      };
    legend.addTo(map);
    
};

function percentLegends(breaks, colors, header, unit){
    
    legend2.onAdd = function (map) {      
      var div = L.DomUtil.create('div', 'info legend'),
          labels = ['<strong>' + header + '</strong>'];
        
        //iterate through grades and create a color field and label for each
        for (var i = 0; i < breaks.length; i++) {
            labels.push(
                (breaks[i + 1]) ? '<i style="background:' + getColor(breaks[i+1]-.01, breaks, colors) + '"></i> ' +
                (parseFloat(breaks[i] * 100).toFixed(0).toLocaleString() + '&ndash;' + parseFloat(breaks[i + 1] * 100).toFixed(0).toLocaleString()) + ' ' + unit : '');
        }
        
        
        div.innerHTML = labels.join('<br>');
        return div;
      };
    legend2.addTo(map);
    
};

function regLegends(breaks, colors, header, unit){
    
    legend.onAdd = function (map) {      
      var div = L.DomUtil.create('div', 'info legend'),
          labels = ['<strong>' + header + '</strong>'],
          break1 = '<i style="background:' + colors[0] + '"></i> ' +
                unit + ' < ' + breaks[0],
          break2 = '<i style="background:' + colors[1] + '"></i> ' +
                breaks[0] + ' < ' + unit + ' < ' + breaks[1],
          break3 = '<i style="background:' + colors[2] + '"></i> ' +
                breaks[1] + ' < ' + unit + ' < ' + breaks[3],
          break5 = '<i style="background:' + colors[4] + '"></i> ' +
                breaks[3] + ' < ' + unit + ' < ' + breaks[4],
          break6 = '<i style="background:' + colors[5] + '"></i> ' +
                unit + ' > ' + breaks[4];
        
        //push colors and labels to an array
        labels.push(
                break1, break2, break3, break5, break6
        );
        
        
        div.innerHTML = labels.join('<br>');
        return div;
      };
    legend.addTo(map);
    
};


$("#calculate").click(function(){
    interpolateData(Weight, binArea);
    document.getElementById('runregression').style.display = 'block';
    document.getElementById('refresh').style.display = 'block';

});

$("#runregression").click(function(){
    joinFeatures(Weight, binArea);
    
});

$("#refresh").click(init);

//when the page loads, AJAX & call createMap to render map tiles and data.
$(document).ready(init);
function init(){
    //globalOutput = document.querySelector("header output");
    document.getElementById('runregression').style.display = 'none';
    document.getElementById('refresh').style.display = 'none';
    
    //remove old data and clean up display
    
    nitrateIDWLayer.clearLayers();
    regressionLayer.clearLayers();
    idwLayer.clearLayers();
    censusLayer.remove();
    wellLayer.remove();
    idwLayer.remove();
    regressionLayer.remove();

    createMap();
  	//create map home button
  	$("header button").on("click", function(){
    	globalMap.flyTo([43.78, -88.78], 5.5); //[lat, lng], zoom
    });
};


