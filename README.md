# utci repo  
Calculates thermal comfort index UTCI. Model for UTCI calculation stands as a separate package outside of Deode-Workflow, EDF and OPTI_THRED where it is used. UTCI model is written in two versions:  
<ul>
  <li> as a class in calculate_utci.py</li>
  <li> as a callable function utci_function.py </li>
</ul>  
First version is used in Deode-Workflow and second in EDF and OPTI-THRED.

## Use in Deode-Workflow
UTCI it is impotred in Deode-Workflow as a module through pyproject.toml and used in the gribmodify.py that reads GRIB2 files, takes input parameters needed for UTCI calculation, calls the model and modifies existing GRIB2 files with UTCI fields.  
Parameters needed for the function call:  
The parameters are prescribed in a gribmodify_rules.json file. The json file contains a (single) list of parameter descriptors:    
<ul>
  <li>input_grib_id: a list of key values that uniquely define the GRIB record in strict order</li>  
  <li>operator: name of the calculation method</li>  
  <li>output_name: name of the output variable</li>  
  <li>output_grib_id: key value that uniquely define the GRIB record of the output variable.</li>
</ul>  
Example:<br>
{ "input_grib_id" : [ {"shortName":"2t"}, {"shortName":"mrt"}, {"shortName":"2r"},{"shortName":"10u"},{"shortName":"10v"} ],<br>
  "operator" : "utci",<br>
  "output_name" : "utci",<br>
  "output_grib_id" : { "shortName":"utci" }<br>
}<br>

## Use in EDF and OPTI-THRED  
UTCI is imported in EDF and OPTI_THRED as a library and will call UTCI to determine when the thresholds have been crossed, by what fraction, determine severity levels etc.
