SELECT p.lat as lat, p.lon as lon, s.genus, s.species, COUNT(*) 
FROM FIA_TREE t
JOIN fia_species s, FIA_PLOT p
ON s.spcd = t.spcd AND p.cn = t.plt_cn
WHERE t.invyr = 2010
GROUP BY t.spcd, t.plot
ORDER BY lat, lon