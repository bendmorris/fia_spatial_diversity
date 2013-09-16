SELECT p.lat as lat, p.lon as lon, s.genus, s.species, COUNT(*) 
FROM FIA_TREE t
JOIN fia_species s, FIA_PLOT p
ON s.spcd = t.spcd AND p.cn = t.plt_cn
WHERE t.invyr = 2010
AND p.plot_status_cd = 1
AND p.kindcd > 0 AND p.kindcd < 4
AND (p.designcd = 1 OR p.designcd = 220 OR p.designcd = 240 OR 
     (p.designcd >= 311 AND p.designcd <= 314) OR
     p.designcd = 328 OR
     (p.designcd >= 502 AND p.designcd <= 505)
    )
AND p.qa_status = 1
AND p.manual >= 1
AND p.samp_method_cd = 1
AND t.statuscd = 1
GROUP BY t.spcd, t.plot
ORDER BY lat, lon