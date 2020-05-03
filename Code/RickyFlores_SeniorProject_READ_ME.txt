The following will get you to run the Network Generation R-Script to develop the files used in my BIOI senior project.

Data provided by:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77221

All the necessary files including data are provided in the following git hub repo:
https://github.com/r-flores/Analysis-of-Hormones-and-the-Gut-Microbiome

Download or Clone the Repository to the location of your choosing.
Excecute the Network_Generation.R file this will create a text file which can then be imported to cytoscape under /Code/Results

By default, the code is set to generate one network (Conventionally Raised Mice)
Line 95 Must be modified to reflect which network you wish to generate the choices are
Conventionally Rasied Male Mice = convR
Germ Free Rasied Male Mice = GF
Growth Hormone Treated Male Mice = GH
Conventionally Rasied Female Mice = Female_convR
Germ Free Rasied Female Mice = Female_GF
*Note If you chose to generate a different network you must specify a different cut hight on line 114 

For a detailed breakdown of the code please refer to the comments within the R-script and or follow the walkthrough Provided at:
https://github.com/r-flores/Analysis-of-Hormones-and-the-Gut-Microbiome