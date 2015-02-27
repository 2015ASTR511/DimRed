SELECT a.param_m_h, a.teff, a.logg, a.param_alpha_m, b.ra
FROM aspcapStar as a,apogeeStar as b
WHERE a.apstar_id = b.apstar_id AND b.location_id = 4558 
ORDER BY b.ra

