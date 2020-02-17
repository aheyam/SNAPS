# Flask web app for SNAPS

from flask import Flask, render_template, jsonify, request
#from SNAPS import SNAPS_compute
import sys, os

#main_SNAPS_file_path = os.path.dirname(os.path.realpath(__file__)) + '/../python'
#sys.path.append(main_SNAPS_file_path)
#os.chdir(main_SNAPS_file_path)

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
        
        return jsonify(form_validated=True, validation_errors=[])
    else:
        validation_errors = []
        for field in form:
            for error in field.errors:
                validation_errors += [f"{field.label}: {error}"]
        return jsonify(form_validated=False, validation_errors=validation_errors)

if __name__ == '__main__':
    app.run(debug=True)