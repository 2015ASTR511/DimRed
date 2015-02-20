/*This file contains a query to get the M/H for the targets in our learning data from the CAS. Our query to the SAS put out was sorted on RA, so RA is selected here to verify order.*/

SELECT a.param_m_h, b.ra
FROM aspcapStar as a,apogeeStar as b
WHERE a.apstar_id = b.apstar_id AND b.location_id = 4558 
ORDER BY b.ra
