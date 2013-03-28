SELECT ROUND(p.lat,2) as lat, ROUND(p.lon,2) as lon, s.genus, s.species, COUNT(*) 
FROM FIA_TREE t
JOIN fia_species s, FIA_PLOT p
ON s.spcd = t.spcd AND p.cn = t.plt_cn
WHERE t.invyr <= 2010 AND t.invyr >= 2006
GROUP BY t.spcd, t.plot
ORDER BY lat, lon