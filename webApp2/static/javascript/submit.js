// Javascript controlling submission of input data, display and downloading of results, and error handling

$("#submit-button").click(function() {
    // Clear any old error messages and results
    $("#error-section").hide();
    $("#error-messages").empty();
    $("#results-section").hide();
    $("#downloads").empty();
    $("#plots").empty();

    // Display a loading symbol
    $("#input-section").after("<div class='loader'></div>")

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
        $(".loader").remove();

        if (jsonData.form_validated) {
            // Add the download buttons
            $("#downloads").append(downloadButton("Download results table", jsonData.main_results, "results.txt"));
            $("#downloads").append(downloadButton("Download chemical shift list", jsonData.shiftlist, "shifts.txt"));
            $("#downloads").append(downloadButton("Download HSQC plot", jsonData.hsqc_plot_html, "HSQC.html"));
            $("#downloads").append(downloadButton("Download strip plot", jsonData.strip_plot_html, "strip plot.html"));

            // Add the plots
            $("#plots").append("<div id='hsqc_plot'></div>");
            Bokeh.embed.embed_item(jsonData.hsqc_plot_json, "hsqc_plot");
            $("#plots").append("<div id='strip_plot'></div>");
            Bokeh.embed.embed_item(jsonData.strip_plot_json, "strip_plot");

            $("#results-section").show();
        }
        else {
            // Display the validation errors
            jsonData.validation_errors.forEach(function(error) {
                $("#error-messages").append("<p>"+error+"</p>");
            });
            $("#error-section").show();
        }
    })
    .catch(function(error) {
        $(".loader").remove();

        console.log(error);
        $("#error-messages").append("<p>"+error+"</p>");
        $("#error-section").show();
    });

});

function downloadButton(buttonName, downloadFile, defaultFilename) {
    // Create a download button and an event handler for when it's clicked
    btn = document.createElement("button");
    $(btn).text(buttonName);
    $(btn).click(function() {
        let file = new Blob([downloadFile], { type: "application / octet - stream" });
        saveAs(file, defaultFilename);
    });
    return btn;
};
