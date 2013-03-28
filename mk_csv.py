''' 
This script requires the fia.sqlite database, which is several gigabytes in 
size, so it isn't included in this repository.

After downloading REF_SPECIES.CSV from the FIA website
(http://apps.fs.fed.us/fiadb-downloads/datamart.html),
you can generate it using the EcoData Retriever with the following commands:

retriever install fia -e s -f fia.sqlite
retriever install fia_species.script -e s -f fia.sqlite
''' 
import sqlite3 as dbapi
import csv

con = dbapi.connect('fia.sqlite')
query = open('query.sql').read()
cur = con.cursor()
cur.execute(query)

with open('fia.csv', 'w') as output_file:
    writer = csv.writer(output_file)
    writer.writerow(('lat','lon','genus','species','count'))
    for row in cur:
        writer.writerow(row)