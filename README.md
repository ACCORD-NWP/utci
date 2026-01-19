utci
Calculates thermal comfort UTCI from grib files.

25-05-2025 ines.muic@cirus.dhz.hr

Parameter list
The parameters are prescribed in a .json file. There is a default file that shows the main structures available. The json file should contain a (single) list of parameter descriptors explained below.

Simple parameters
Basic variables that can be read from a single GRIB record, can be described by a list with the following 4 elements

input_grib_id: a list of key values that uniquely define the GRIB record in strict order
operator: name of the calculation method
output_name: name of the output variable
output_grib_id: key value that uniquely define the GRIB record of the output variable.

{ "input_grib_id" : [ {"shortName":"2t"}, {"shortName":"mrt"}, {"shortName":"2r"},{"shortName":"10u"},{"shortName":"10v"} ],
  "operator" : "utci",
  "output_name" : "utci",
  "output_grib_id" : { "shortName":"utci" }
}
