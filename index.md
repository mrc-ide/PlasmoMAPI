# PlasmoMAPI
### Version 0.1.0

--------------------------------------------------------------------------------------------------------------------------------

PlasmoMAPI is an R package for converting pairwise spatial data (i.e. statistical values calculated between different sampling locations) into maps that show areas of high or low connectivity*. These maps can be useful for 1) visualising data, as an alternative to simple networks of nodes and edges which get confusing when when the number of nodes is large, and 2) statistical testing, because it can be very difficult to detect subtle patterns when looking at raw data. PlasmoMAPI is a spinoff from the original MAPI method (REF) by Sylvian Piry and co-authors, differing in a number of key assumptions and also being tailored toward analysis of malaria (Plasmodium) data.

<img src="https://raw.githubusercontent.com/mrc-ide/PlasmoMAPI/master/R_ignore/web_images/main_page.png" height="800px" width="800px" />

If you'd like to see an example analysis from start to finish, check out the 

* for PlasmoMAPI to give meaningful results it must be that pairwise data are connected in some way to a spatial process. A good example would be pairwise genetic distances calculated between samples extracted at different points in space.
