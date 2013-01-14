#!/usr/bin/env python
#
# usage: general_query -h
#        general_query.py -tablelist
#        general_query.py [-table <tablename>] -fieldlist
#        general_query.py [arglist] (see -h)
#
#        e.g. general_query.py -table cresis_gr \
#                   -fields "wgs84surf,wgs84bed" -epsg 32622 \
#                   -and_clause "wgs84surf>-500" \
#                   -box -80.0 10.0 59.0 83.0 \
#                   -year_range 2010 2011 -mod_field gid \
#                   -mod_val
#
# input file is lon lat pairs one per line, first will be appended to
# end to close polygon
# 
# The idea is that we construct a polygon using a lat,lon box then
# extract all the points from the table within that box.  We might 
# also want to sort by track number and time ....
#
import os 
import sys
import psycopg2
import pyproj
import argparse
import re

dbname='icedb'
dbhost='localhost'

# set up command line arguments
parser = argparse.ArgumentParser( \
    description='Make database query for points in polygon.',
    epilog='>> either -in_poly_file OR -box MUST be specified <<',
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-tablelist', 
                    action='store_true', 
                    help='list the tables available in the database and exit')
parser.add_argument('-table', 
                    action='store', 
                    type=str, 
                    default='atm_icessn', 
                    help='the table in the database to query')
parser.add_argument('-fieldlist', 
                    action='store_true', 
                    help='list the fields available in the specified table and exit')
parser.add_argument('-fields', 
                    action='store', 
                    type=str, 
                    metavar='field1,field2...(no spaces)',
                    default='gid,track,filename,fltdate,timetag,height,height_str,' + \
                        'sn_slope,we_slope,rms_fit,nused,nedit,cbcl_dist,trkid', 
                    help='the fields to return')
parser.add_argument('-and_clause', 
                    action='store', 
                    type=str, 
                    default='', 
                    help='add an AND clause to the query:  trkid=0 (no spaces)')
parser.add_argument('-epsg_out', 
                    action='store', 
                    type=int, 
                    default=4326, 
                    help='the EPSG projecion number for output x y ' + \
                        'coordinates - default 4326 (lon lat)')
parser.add_argument('-year', 
                    action='store', 
                    type=str, 
                    default='None', 
                    help='the year of acquisition - default all years')
parser.add_argument('-year_range', 
                    action='store', 
                    type=str, 
                    metavar=('start_year','stop_year'), 
                    default=('None','None'), 
                    nargs=2, 
                    help='year range')
parser.add_argument('-all_fields', 
                    action='store_true', 
                    help='display all fields in record - default is display' + \
                        ' lon lat x y bed_elev surf_elev')
parser.add_argument('-in_poly_file', 
                    action='store', 
                    metavar='input_poly_file', 
                    type=str, default='None', 
                    help='input text file with query polygon (lon ' + \
                        'lat pairs on separate lines)')
parser.add_argument('-box', 
                    action='store', 
                    type=float, 
                    metavar=('min_lon','max_lon','min_lat','max_lat'), 
                    default=[0.0, 0.0, 0.0, 0.0], 
                    nargs=4, 
                    help='longitude latitude box for query included in command line')
parser.add_argument('-mod_val', 
                    action='store', 
                    type=int, 
                    metavar='decimation_factor', 
                    default=0, 
                    help='factor by which to decimate results (keep ' + \
                        'record_num mod(mod_val) == 0)')
parser.add_argument('-mod_field', 
                    action='store', 
                    type=str, 
                    default='gid', 
                    help='the field to use for modulus decimation')
parser.add_argument('--port',
                    type=int,
                    default=5433, 
                    help='ssh port')

args = parser.parse_args()
port = args.port

# set up connection to database
conn = psycopg2.connect("dbname=%s host=%s user=nobody port=%s" % (dbname, dbhost, str(port)))
curs=conn.cursor()

# list tables
if args.tablelist :
    
    # query table names
    qs2="SELECT tablename FROM pg_tables WHERE tablename !~* 'pg_*' AND tablename !~* 'sql_';"
    curs.execute(qs2)

    # extract results
    row_of_tables = curs.fetchall()
    print "tables in database %s" % (dbname,)
    print "==============================="
    for tname in row_of_tables:
        print tname[0]
    exit()

elif args.fieldlist :
    # fmtstr="SELECT a.attname AS field FROM pg_class c, pg_attribute a, " + \
    #     "pg_type t WHERE c.relname = '%s' AND a.attnum > 0 AND " + \
    #     "a.attrelid = c.oid AND a.atttypid = t.oid ORDER BY a.attnum;" 
    # qs1 = fmtstr % (args.table,)
    fmtstr="select ordinal_position, column_name, data_type " + \
        "from information_schema.columns where table_schema='public' " + \
        "and table_name='%s' order by ordinal_position;" 
    qs1 = fmtstr % (args.table,)
    curs.execute(qs1)

    # extract results
    rows = curs.fetchall()
    # print "# Number of rows fetched: ",len(rows)
    # print "#>>>>>>>>",rows[0],"<<<<<<<<<<<"
    # print rows
    print "fields in table %s " % (args.table,)
    print "===================================="
    print "name                | type"
    print "--------------------+---------------"
    for field in rows:
        fieldname=field[1];
        fieldtype=field[2];
        print "%-20s| %-20s" % (fieldname, fieldtype,)
    exit()


epsg_str=str(args.epsg_out)


if  args.year == 'None':
  if args.year_range[0] == 'None':
      year_str='' # all years uses a blank year string in the sql call
      years='all'
  else:
      year_start=args.year_range[0]
      year_stop=args.year_range[1]
      years='range'
else:
  if args.year_range[0] != 'None':
        print "Error: cannot define both -year and -year_range; exiting"
        sys.exit(0)
  else:
        year_str=args.year
        years='single'

#
# if -all_fields flag thrown, output all fields in record; else just
# location and surface and bed elevation
#
out_full_record=args.all_fields

# put the command line into the output file as a comment, so we know
# what was requested
print ('# ' + ' '.join(sys.argv[0:]))
sys.stderr.write(' '.join(sys.argv[0:]) + '\n')


if args.in_poly_file == 'None' :
    if (args.box[0]==0) & (args.box[1]==0) :
        print "Error: either -in_poly_file or -box must be specified; exiting"
        sys.exit(0)
    else:
        min_lon=args.box[0]
        max_lon=args.box[1]
        min_lat=args.box[2]
        max_lat=args.box[3]
        str1= \
             '%f %f\n %f %f\n %f %f\n %f %f' % \
             (min_lon,min_lat,min_lon,max_lat,max_lon,max_lat,max_lon,min_lat)  # box will be closed below
else:
    if args.box[0]!=args.box[1] :
        print "Error: both -in_poly_file and -box cannot be specified; exiting"
        sys.exit(0)
    else:
        in_poly_file=args.in_poly_file
        # read in the lon lat\n poly file - try to deal with linux/mac and pc line endings
        #
        f=open(in_poly_file,'r')
        str1=f.read().strip()                           # read it in, dump empty lines, etc
#
# now get rid of blanks
#
splitstr='\n'                                   # get ready for a mac file, but
if '\r\n' in str1: 
    splitstr='\r\n'                             # if input file came from a pc

# 
# The input polygon is assumed not to be closed, so we repeat the
# first coordinate pair at end to close it (one could test for this
# and do the right thing in both cases, but not implimented yet
#
# add first pair to end to close poly
out_pts_str=str1.replace(splitstr,', ') + \
    ', ' + str1[:str1.find(splitstr)]  

poly_str="PolyFromText(\'POLYGON ((" + out_pts_str + "))\',4326)"

print ('# ' + poly_str)

#
# set up connection to database
conn = psycopg2.connect("dbname=%s host=%s user=nobody port=%s" % (dbname, dbhost, str(port)))
# 
# setup query in parts
#
#  fields_str is the fields the query will return from the table
#
# fields_str = "lon, lat, thick, surface, bottom, elevation, quality, " + \
#     "frame, flight_date, gid, stime" 
#
# fields_str = "gid, track, filename, " + \
#     "fltdate, timetag, lon, lat, height, height_str, sn_slope, we_slope, " + \
#     "rms_fit, nused, nedit, cbcl_dist, trkid"

#
# argparse does not like spaces in strings (like the fields specifier)
# - actually, shell qoutes for this are problematic too - so don't
# use spaces...  but we need to put the spaces back in here to get
# things to format well for output later.
#
regpattern=re.compile(',')
args_fields_with_spaces=regpattern.sub(', ',args.fields)

fields_str = "lon, lat, " + args_fields_with_spaces

#
#
# table_name='atm_icessn'
table_name=args.table

if years=='all':
  query_str = "select %s from %s where Within(the_geom, %s);" % \
    (fields_str, table_name, poly_str)
elif years=='single':
  query_str = \
      "select %s from %s where Within(the_geom, %s) AND (extract ( YEAR from fltdate ) = %s );" % \
      (fields_str, table_name, poly_str, year_str)
elif years=='range':
    query_str = \
        "select %s from %s where Within(the_geom, %s) AND (extract ( YEAR from fltdate ) >= %s ) AND (extract ( YEAR from fltdate ) <= %s );" % \
    (fields_str, table_name, poly_str, year_start, year_stop)

if args.mod_val>0:
   query_str=query_str[:-1] + "AND ((%s %% %s) = 0);" % (args.mod_field,str(args.mod_val))

if args.and_clause <> '':
   query_str=query_str[:-1] + "AND (%s);" % (args.and_clause)
#
# open cursor for select
curs=conn.cursor()
curs.execute(query_str)
sys.stderr.write("# select string: " + query_str + '\n')

#
# extract results
#
rows = curs.fetchall()
print "# Number of rows fetched: ",len(rows)
if len(rows)>0:
	print "#>>>>>>>>",rows[0],"<<<<<<<<<<<"
sys.stderr.write("# Number of rows fetched: " + str(len(rows)) + '\n')

#
# print header in output as a comment
#
if out_full_record == True:
    print  "lon, lat, x, y, gid, track, filename, fltdate, timetag, height, height_str, sn_slope, we_slope, rms_fit, nused, nedit, cbcl_dist, trkid"
else:
    print "x, y, " + fields_str

# # clean up database connection
curs.close()
conn.close()
#
# add Projection X, Y (in meters?)
p1=pyproj.Proj(init="epsg:4326")

p2=pyproj.Proj(init="epsg:%s" % (epsg_str))
#
#
# output gives you: fields_str = "gid, track, filename, fltdate, timetag, lon, lat, height, height_str, sn_slope, we_slope, rms_fit, nused, nedit, cbcl_dist, trkid"
#
for irow in range(len(rows)):

#    gid=int(rows[irow][0])
#    track=rows[irow][1]
#    filename=rows[irow][2]
#    fltdate=rows[irow][3]
#    timetag=float(rows[irow][4])
#    height=float(rows[irow][7])
#    height_str=rows[irow][8]
#    sn_slope=float(rows[irow][9])
#    we_slope=float(rows[irow][10])
#    rms_fit=float(rows[irow][11])
#    nused=int(rows[irow][12])
#    nedit=int(rows[irow][13])
#    cbcl_dist=float(rows[irow][14])
#    trkid=int(rows[irow][15])
    x1=float(rows[irow][0])
    y1=float(rows[irow][1])
    x2, y2 = pyproj.transform(p1,p2,x1,y1)
#    if out_full_record == True:
#        outstr = "%12.6f, %12.6f, %12.2f, %12.2f, %d, %s, %s, %s, %12.8f, %12.6f, %s, %12.8f, %12.8f, %12.8f, %d, %d, %f, %d" % \
#            (x1, y1, x2, y2, gid, track, filename, fltdate, timetag, height, height_str, sn_slope, we_slope, rms_fit, nused, nedit, cbcl_dist, trkid)
#    else:
#        outstr = "%12.6f, %12.6f, %12.2f, %12.2f, %12.4f, %d" % \
#            (x1, y1, x2, y2, height, trkid)
#    outstr = "%12.2f, %12.2f, %s" % \
#            (x2, y2, rows[irow][2:])
#    print outstr
#    print rows[irow]
    print x2,",",y2,",",
    for iout in range(len(rows[irow][:-1])):
        print rows[irow][iout],",",
    print rows[irow][-1]

