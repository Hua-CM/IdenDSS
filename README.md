# IdenDSS

This is a software for identifying taxon-specific signature sequence from a DNA database. Pleas run IDSS.py in [src](src) directory.

The meta file should contain three columns:  

| column | description |  
| --- | --- |  
| Group | Usually species name |  
| Sample | Sample id, e.g. sample1,sample2,... |  
| Sequence | Sequence header |

The result file contain five colums

| column | description |  
| --- | --- |  
| group | Usually species name |  
| assembly | The reference assembly |  
| seq | DSS sequence |
| position | DSS position on reference assembly |
| GC |GC content|

For more help see example files and help option.



