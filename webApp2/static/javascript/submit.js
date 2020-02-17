// Javascript controlling submission of input data, display and downloading of results, and error handling

$("#submit-button").click(function() {
    // Clear any old error messages and results
    $("#error-section").hide();
    $("#error-messages").empty();
    $("#results-section").hide();
    $("#downloads").empty();
    $("#plots").empty();

    // POST the form contents
    let request = {
        method: "POST",
        body: new FormData($("#input-form")[0])
    };
    fetch("/submit", request)
    .then(function(response) {
        return response.json();
    })
    .then(function(jsonData) {
        if (jsonData.form_validated) {
            $("#downloads").append("<p>Downloads will be shown here</p>");

            $("#plots").append("<div id='hsqc_plot'></div>");
            Bokeh.embed.embed_item(jsonData.hsqc_plot, "hsqc_plot");

            $("#plots").append("<div id='strip_plot'></div>");
            Bokeh.embed.embed_item(jsonData.strip_plot, "strip_plot");

            $("#results-section").show();
        }
        else {
            // Display the validation errors
            jsonData.validation_errors.forEach(function(error) {
                $("#error-messages").append("<p>"+error+"</p>");
            });
            $("#error-section").show();
        }
    });

});
