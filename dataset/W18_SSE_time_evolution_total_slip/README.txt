MODL_time_JJJ.atr - time window slip distributions on fault

Files contain fault attributes and quadrilaterals that can be used to make color plots of slip distributions.

First line in file is header:

# Time_1   Time_2  List_of_events
The header line for each fault segment has multiple attributes following the -Z; use awk to select the one to plot.

awk '{ if ($1 ==">") print $1, $2, $4; else print $1, $2 }' MODL_time_001.atr  | psxy -Cpalette.cpt .... 
where $1 = '>', $2 = '-Z' and the attributes are (in order):

3 time window index (from DW: option) 
4 total slip amplitude (mm) 
5 strike-slip component (Okada U1; see Fig. 2)
6 dip-slip component (U2)
7 fault opening component (U3)
8 Rake angle on sub-segment
9 Sub-segment centroid longitude
10 Sub-segment centroid latitude
11 Sub-segment centroid depth (in km, positive downward)
12 Node X-index
13 Node Z-index
14 Sub-segment X-index
15 Sub-segment Z-index
16 Area of patch (km^2)
Followed by coordinates of the trapezoid:
Lon1 Lat1 Depth1   
.
.
Lon4 Lat4 Depth4   
 
> -Z   1     0.006     0.000     0.006     0.000    90.000   265.502    14.802     -5.10   2   2   1   1
    265.52165000     14.78570000     -5.00000000
    265.47878889     14.81039444     -5.00000000
    265.48159022     14.81915733     -5.20000000
    265.52448000     14.79450400     -5.20000000