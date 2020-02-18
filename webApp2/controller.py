# Flask web app for SNAPS

from flask import Flask, render_template, jsonify
from bokeh.embed import json_item, file_html
from bokeh.resources import INLINE
import sys, os

# Need to add the SNAPS python directory to the path before importing
main_SNAPS_file_path = os.path.dirname(os.path.realpath(__file__)) + '/../python'
sys.path.append(main_SNAPS_file_path)
from SNAPS import SNAPS_compute

app = Flask(__name__)
app.config["SECRET_KEY"] = "test"

@app.route("/")
def index():
    "Render the front page, incuding the submission form."
    from SNAPS_input_form import SNAPS_input_form
    
    form = SNAPS_input_form()
    
    return render_template("index.html", form=form)

@app.route("/info")
def info():
    return render_template("info.html")

@app.route("/howto")
def howto():
    return render_template("howto.html")

@app.route("/submit", methods=["POST"])
def submit():
    """Receive submitted data, validate it, and return results.
    If submitted data doesn't validate, return relevant errors."""
    from SNAPS_input_form import SNAPS_input_form
    
    form = SNAPS_input_form()
    
    if form.validate():
        
        results_tables, plots = SNAPS_compute(form.shift_file.data,
                                              form.shift_type.data,
                                              form.pred_file.data,
                                              form.pred_type.data)
        
        hsqc_plot_json = json_item(plots["hsqc_plot"])
        strip_plot_json = json_item(plots["strip_plot"])
        hsqc_plot_html = file_html(plots["hsqc_plot"], resources=INLINE)
        strip_plot_html = file_html(plots["strip_plot"], resources=INLINE)
        
        return jsonify(form_validated=True, validation_errors=[],
                       main_results=results_tables["main_results"],
                       shiftlist=results_tables["shiftlist"],
                       hsqc_plot_json=hsqc_plot_json, 
                       hsqc_plot_html=hsqc_plot_html,
                       strip_plot_json=strip_plot_json, 
                       strip_plot_html=strip_plot_html)
    else:
        validation_errors = []
        for field in form:
            for error in field.errors:
                validation_errors += [f"{field.label}: {error}"]
        return jsonify(form_validated=False, validation_errors=validation_errors)

if __name__ == '__main__':
    app.run(debug=True)
