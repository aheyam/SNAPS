# Definition and validation of the SNAPS input form

# Nb. Validation of FileField from the standard wtforms module doesn't seem to 
# work, so I've used the equivalent from flask_wtf.

from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired
from wtforms import SelectField
from wtforms.validators import InputRequired, ValidationError

class SNAPS_input_form(FlaskForm):
    shift_file = FileField(label="Observed chemical shifts file",
                               validators=[FileRequired()]) # InputRequired()?
    shift_type = SelectField(label="Observed shifts format",
                             choices=[("snaps", "SNAPS"), 
                                      ("ccpn", "CCPN"), 
                                      ("sparky", "Sparky"), 
                                      ("mars", "MARS"),
                                      ("xeasy", "Xeasy"), 
                                      ("nmrpipe", "NMRpipe"), 
                                      ("test", "Test")],
                             validators=[InputRequired()])
    pred_file = FileField(label="Predicted chemical shifts file",
                               validators=[FileRequired()]) # InputRequired()?
    pred_type = SelectField(label="Predicted shifts format",
                            choices=[("shiftx2","SHIFTX2"),
                                     ("sparta+","SPARTA+")],
                            validators=[InputRequired()])
    
    def validate_shift_file(form, field):
        """Custom validation of the shift_file field
        This is run automatically as part of FlaskForm.validate()"""
        
        # placeholder for input file validation
        pass
        
    def validate_pred_file(form, field):
        """Custom validation of the pred_file field 
        This is run automatically as part of FlaskForm.validate()"""
        
        # placeholder for input file validation
        pass
    




