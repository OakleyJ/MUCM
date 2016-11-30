// hover_both.js : Enable paired hovering of parameters and objectives plots in optim_results_app.R
//


// NOTE: would be good to have paired selections too...

function handleHover(value, plot) {
    if (value === null) {
        Plotly.Fx.unhover(plot);
    } else {
        Plotly.Fx.hover(plot, [
            { curveNumber: 0, pointNumber: value[0].pointNumber },
        ]);
    }
}

$(document).on("shiny:inputchanged", function(event) {
    if (event.name == ".clientValue-plotly_hover-Approx.vs.Truth") {
        //console.log(event.name + " changed");
        var value = eval(event.value);
        handleHover(value, "qqplot");
        handleHover(value, "pcd.plot");
        
    } else if (event.name == ".clientValue-plotly_hover-qqplot") {
        //console.log(event.name + " changed");
        var value = eval(event.value);
        handleHover(value, "Approx.vs.Truth");
        handleHover(value, "pcd.plot");
        
    } else if (event.name == ".clientValue-plotly_hover-pcd.plot") {
        //console.log(event.name + " changed");
        var value = eval(event.value);
        handleHover(value, "Approx.vs.Truth");
        handleHover(value, "qqplot");
        
    }/* else
        console.log(event.name + " changed"); */
});

